configfile: "config.yaml"
import os

pipeline = "mapping" 

include: "rules/create_file_log.smk"

workdir: config["OUTDIR"]

ASSEMBLY = config["ASSEMBLY"]
MULTIPATH = config["MULTIPATH"] 
OTHER_PATH = config["OTHER_PATH"]


FQ_TO_MAP_dict={}
ALL_SAMPLES_dict={}
SINGLE_SAMPLES = []

for multipath_dir in os.listdir(MULTIPATH):
    if os.path.isdir(os.path.join(MULTIPATH , multipath_dir)):
        for intermediate_dir in os.listdir(MULTIPATH + multipath_dir):
            sample_dir = os.path.join(MULTIPATH, multipath_dir, intermediate_dir)
            if os.path.isdir(sample_dir) and not sample_dir.endswith("report"):
                for sample in os.listdir(sample_dir):
                    new_path = os.path.join(sample_dir,sample)
                    if os.path.isdir(new_path):
                        files_in_dir = [os.path.join(new_path,f) for f in os.listdir(new_path) if f.endswith('.fq.gz')]
                        files_in_dir = sorted(files_in_dir)
                        ALL_SAMPLES_dict[sample] = files_in_dir
                        if len(files_in_dir) >2:
                            FQ_TO_MAP_dict[sample+"_1"] =  files_in_dir[0:2]
                            FQ_TO_MAP_dict[sample+"_2"] =  files_in_dir[2:4]
                        elif len(files_in_dir) <=2:
                            FQ_TO_MAP_dict[sample] =  files_in_dir
                            SINGLE_SAMPLES.append(sample)


### for reads in another directory
samp_name_otherpath = OTHER_PATH.rsplit('/',1)[1]
otherpath_fq = [os.path.join(OTHER_PATH,f) for f in os.listdir(OTHER_PATH) if f.endswith('.fq.gz')]
otherpath_fq = sorted(otherpath_fq)
if len(otherpath_fq) >2:
    FQ_TO_MAP_dict[samp_name_otherpath+"_1"] =  otherpath_fq[0:2]
    FQ_TO_MAP_dict[samp_name_otherpath+"_2"] =  otherpath_fq[2:4]
elif len(otherpath_fq) <=2:
    FQ_TO_MAP_dict[samp_name_otherpath] =  otherpath_fq
    SINGLE_SAMPLES.append(samp_name_otherpath)

ALL_SAMPLES_dict[samp_name_otherpath] = otherpath_fq


# important variables :
# - ALL_SAMPLES_dict - dictionary with sample IDs as keys and list with path to the fastq as value
# - SINGLE_SAMPLES - list with name of samples that only have two fq files
# - FQ_TO_MAP_dict - dictionary where samples that have two pairs of fq files are split into SAMPLE_1 and SAMPLE_2
#                   keys are the sample ids (with _1 or _2 if they need to be run separately) and values are fq files

localrules: create_file_log
rule all:
    input:
        files_log,
        expand("processed_reads/{sample_merged}.sorted.bam.bai", sample_merged=ALL_SAMPLES_dict.keys()),
        expand("mapping_stats/qualimap/{sample_merged}/genome_results.txt", sample_merged=ALL_SAMPLES_dict.keys()),

def create_input_names(wildcards):
    return(FQ_TO_MAP_dict[wildcards.sample])

def create_names_to_merge(wildcards):
    input_files = []
    if wildcards.sample_merged in SINGLE_SAMPLES:
        input_files.append("mapped_reads/"+wildcards.sample_merged+".bam")
    else:
        input_files.append("mapped_reads/"+wildcards.sample_merged+"_1.bam")
        input_files.append("mapped_reads/"+wildcards.sample_merged+"_2.bam")
    return(input_files)

rule bwa_map:
    input:
        ASSEMBLY,
        create_input_names
    output:
        temp("mapped_reads/{sample}.bam")
    message:
        "Rule {rule} processing"
    params:
        readgroup = lambda wildcards: wildcards.sample.split("_")[0]
    shell:
        """
        module load bwa
        bwa mem -t 16 -R "@RG\\tID:{params.readgroup}\\tSM:{params.readgroup}" {input} | samblaster -r | samtools view -b - > {output}
        """


rule merge_mapped:
    input:  
        create_names_to_merge
    output:
        "merged_reads/{sample_merged}.bam"
    message:
        "Rule {rule} processing"
    run:
        if len(input) >1:
            if input[0].endswith("_1.bam") and input[1].endswith("_2.bam"):
                shell("samtools merge -@ 16 {output} {input}")
        else:
            shell("mv {input} {output}")


rule samtools_sort:
    input: 
        rules.merge_mapped.output
    output: 
        "processed_reads/{sample_merged}.sorted.bam"
    message:
        "Rule {rule} processing"
    # group:
    #     "end"
    shell: 
        "samtools sort -m 2G -@ 7 -O bam {input} > {output}"


rule samtools_index:
    input:
        rules.samtools_sort.output
    output:
        "processed_reads/{sample_merged}.sorted.bam.bai"
    message:
        "Rule {rule} processing"
    # group:
    #     "end"
    shell:
        "samtools index -@ 16 {input}"

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
    resources:
        time="2:0:0"
    shell: 
        "unset DISPLAY && qualimap bamqc -bam {input.bam} --java-mem-size=16G -nt 8 -outformat PDF -outdir {params.outdir}"