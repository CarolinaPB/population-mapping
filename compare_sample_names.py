import os

PATH="/lustre/nobackup/WUR/ABGC/shared/Chicken/Africa/"

dirs = ["X201SC20031230-Z01-F001_multipath/X201SC20031230-Z01-F001_1/raw_data",
        "X201SC20031230-Z01-F001_multipath/X201SC20031230-Z01-F001_2/raw_data",
        "X201SC20031230-Z01-F001_multipath/X201SC20031230-Z01-F001_3/raw_data",
        "X201SC20031230-Z01-F001_multipath/X201SC20031230-Z01-F001_4/raw_data",
        "X201SC20031230-Z01-F004/raw_data"]


dirs_path = [os.path.join(PATH,p) for p in dirs]

samples_list=[]
for sample_dir in dirs_path:
    samples = [samp for samp in os.listdir(sample_dir) if os.path.isdir(os.path.join(sample_dir,samp))]
    samples_list.append(samples)


for i in range(len(samples_list)):
    for n in range(len(samples_list)):
        if i != n:
            print(f"set{i} -- set{n}")
            print(set(samples_list[i]).intersection(samples_list[n]))
            print([x for x in samples_list[i] if x in samples_list[n]])