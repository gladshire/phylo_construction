import os
import sys
import subprocess

# Define post-assembly function

#   1. Convert fastq / fq to fasta (fq_to_fa.py)
#   2. Run chimera detection (chimera_detect.py)
#   3. Run corset, extract representative (corset_wrapper.py)
#   4. Filter corset output (filter_corset_output.py)
#   5. Run translation with transrate (transrate_wrapper.py)

def post_assembly(out_dir):
    
    for filename in os.listdir("05-filter_over_represented"):
        file_path, file_name = os.path.split(filename)
        fq_to_fa_cmd = ["python3", "fq_to_fa.py", "05-filter_over_represented/" + filename, out_dir]
        subprocess.run(" ".join(fq_to_fa_cmd), shell = True)




if __name__ == "__main__":
    post_assembly("/scratch/miles/allele_phasing/phylo_construction/scripts")
