configfile: "config.yaml"
import os
from snakemake.utils import makedirs

pipeline = "mapping" 

include: "rules/create_file_log.smk"

if "OUTDIR" in config:
    print("\nSaving to " + config["OUTDIR"] + "\n")
    workdir: config["OUTDIR"]

makedirs("logs_slurm")

ASSEMBLY = config["ASSEMBLY"]

def define_to_map_together(given_path):
    TO_MAP_TOGETHER = []
    FULLPATHS=[]
    SAMPLES=[]
    for root, dirs, files in os.walk(given_path):
        for name in files:
            if name.endswith("fq.gz"):
                fullpath = os.path.join(root, name)

                sample = fullpath.rsplit("/", 2)[1]
                name_sub = name.rsplit("_", 2)[0]

                if name_sub not in TO_MAP_TOGETHER:
                    TO_MAP_TOGETHER.append(name_sub)
                FULLPATHS.append(fullpath)
                if sample not in SAMPLES:
                    SAMPLES.append(sample)
    return(TO_MAP_TOGETHER, FULLPATHS, SAMPLES)

MAP = {} # to be mapped together
def create_bwa_map_input(given_path, map_dict):
    maptogether, fullpaths, samples=define_to_map_together(given_path)

    for var in maptogether:
        temp = []
        for file in fullpaths:  
            if var in file:
                temp.append(file)
        map_dict[var] = temp
    return(map_dict, samples)

ALL_SAMPLES = []
for file in config["PATHS_WITH_FILES"].values():
    MAP, samples = create_bwa_map_input(file, MAP)
    ALL_SAMPLES = ALL_SAMPLES + samples


localrules: create_file_log
rule all:
    input:
        files_log,
        expand("mapping_stats/qualimap/{sample_merged}/genome_results.txt", sample_merged=ALL_SAMPLES),
        "mapping_stats/sample_quality_summary.tsv"

def create_input_names(wildcards):
    return(MAP[wildcards.sample])

def create_names_to_merge(wildcards):
    input_files = []
    mapped_dir = expand("mapped_reads/{sample}.bam", sample=MAP.keys())
    for file in mapped_dir:
        file_nodir = file.rsplit("/")[1]
        if file_nodir.startswith(wildcards.sample_merged):
            input_files.append(file)
    return(input_files)

rule bwa_map:
    input:
        assembly = ASSEMBLY,
        reads = create_input_names
    output:
        temp("mapped_reads/{sample}.bam")
    message:
        "Rule {rule} processing"
    group:
        'group'
    shell:
        """
module load bwa samtools

rg="$(basename $(dirname "{input.reads[0]}"))"
echo $rg
bwa mem -t 16 -R "@RG\\tID:$rg\\tSM:$rg" {input} | samblaster -r | samtools view -b - > {output}
        """

rule merge_mapped:
    input:  
        create_names_to_merge
    output:
        temp("merged_reads/{sample_merged}.bam")
    message:
        "Rule {rule} processing"
    group:
        'group'
    run:
        if len(input) >1:
            shell("module load samtools && samtools merge -@ 16 {output} {input}")
        else:
            shell("mv {input} {output}")


rule samtools_sort:
    input: 
        rules.merge_mapped.output
    output: 
        "processed_reads/{sample_merged}.sorted.bam"
    message:
        "Rule {rule} processing"
    group:
        'group'
    shell: 
        "module load samtools && samtools sort -m 2G -@ 7 -O bam {input} > {output}"


rule samtools_index:
    input:
        rules.samtools_sort.output
    output:
        "processed_reads/{sample_merged}.sorted.bam.bai"
    message:
        "Rule {rule} processing"
    group:
        "group"
    shell:
        "module load samtools && samtools index -@ 16 {input}"

rule qualimap_report:
    input: 
        check=rules.samtools_index.output, 
        bam=rules.samtools_sort.output
    output: 
        "mapping_stats/qualimap/{sample_merged}/genome_results.txt"
    params:
        outdir = "mapping_stats/qualimap/{sample_merged}/"
    message:
        "Rule {rule} processing"
    group:
        "group"
    shell: 
        "unset DISPLAY && qualimap bamqc -bam {input.bam} --java-mem-size=16G -nt 8 -outformat PDF -outdir {params.outdir}"

rule qualimap_summary:
    input:
        expand("mapping_stats/qualimap/{sample_merged}/genome_results.txt", sample_merged = ALL_SAMPLES)
    output:
        "mapping_stats/sample_quality_summary.tsv"
    message:
        'Rule {rule} processing'
    group:
        "group"
    params:
        scripts_dir = os.path.join(workflow.basedir, "scripts/")
    shell:
        'sh {params.scripts_dir}create_qualimap_summary.sh'