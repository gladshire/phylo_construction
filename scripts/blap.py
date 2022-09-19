import os


for file1 in os.listdir("/scratch/miles/allele_phasing/phylo_construction/scripts/05-filter_over_represented"):
    file_comps = file1.split(".")
    print(file_comps[0].split("_")[2:4])
