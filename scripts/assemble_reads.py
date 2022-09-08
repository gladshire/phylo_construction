import os
import sys
import subprocess
import shutil
import process_reads
import trinity_wrapper
import concat_fasta
from Bio import Entrez

Entrez.email = "mwoodc2@uic.edu"
Entrez.api_key = "6512380b5cb6a15be6ef4919564e2e186508"


def assemble_trinity(processed_dir, threads, max_memory_gb, mult_samples = False):
    curr_dir = os.getcwd()
    if os.path.isabs(curr_dir) == False: curr_dir = os.path.abspath(curr_dir)
    if curr_dir[-1] != "/": curr_dir += "/"

    if os.path.exists(curr_dir + "06-trinity_assembly/") == False:
        os.mkdir(curr_dir + "06-trinity_assembly/")

    files_skip = []
    for filename in os.listdir(processed_dir):
        if filename in files_skip or "comb" in filename:
            continue
        file_comps = filename.split(".")
        file_sra = file_comps[0].split("_")[0]
        layout = process_reads.get_layout(file_sra)
        if mult_samples == True:
            org_name = "_".join(file_comps[0].split("_")[2:-1:])
            files_in = [processed_dir + fq for fq in os.listdir(processed_dir) if org_name in fq]
            
            concat_fasta.assemble_file(fasta_files = files_in,
                                       out_file = org_name + "_comb.fastq",
                                       out_dir = curr_dir + "05-filter_over_represented")
            cmd_trin = ["python3", "trinity_wrapper.py", processed_dir + org_name + "_comb.fastq",
                        str(threads), str(max_memory_gb), "non-stranded", curr_dir + "06-trinity_assembly/"]
            subprocess.Popen(" ".join(cmd_trin), shell = True).wait()
        elif mult_samples == False:
            if layout == "SINGLE":
                cmd_trin = ["python3", "trinity_wrapper.py", processed_dir + filename, str(threads),
                            str(max_memory_gb), "non-stranded", curr_dir + "06-trinity_assembly/"]
                subprocess.Popen(" ".join(cmd_trin), shell = True).wait()
            elif layout == "PAIRED":
                pe_file_1 = file_comps[0][0:-1:] + "1." + ".".join(file_comps[1::])
                pe_file_2 = file_comps[0][0:-1:] + "2." + ".".join(file_comps[1::])

                cmd_trin = ["python3", "trinity_wrapper.py", processed_dir + pe_file_1, processed_dir + pe_file_2,
                            str(threads), str(max_memory_gb), "non-stranded", curr_dir + "06-trinity_assembly/"]
                subprocess.Popen(" ".join(cmd_trin), shell = True).wait()


if __name__ == "__main__":
    curr_dir = os.getcwd()
    if os.path.isabs(curr_dir) == False: curr_dir = os.path.abspath(curr_dir)
    if curr_dir[-1] != "/": curr_dir += "/"

    if len(sys.argv) == 3:
        assemble_trinity(curr_dir + "05-filter_over_represented/", sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4 and "mult-samples" in sys.argv:
        assemble_trinity(curr_dir + "05-filter_over_represented/", sys.argv[1], sys.argv[2], mult_samples = True)
    else:
        print("Usage:")
        print("python3 assemble_reads.py threads max_memory_GB [optional:mult-samples]")
