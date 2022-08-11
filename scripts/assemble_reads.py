import os
import sys
import subprocess
import shutil
import process_reads
import trinity_wrapper
from Bio import Entrez

Entrez.email = "mwoodc2@uic.edu"
Entrez.api_key = "6512380b5cb6a15be6ef4919564e2e186508"


def assemble_trinity(processed_dir, threads, max_memory_gb):
    curr_dir = os.getcwd()
    if os.path.isabs(curr_dir) == False: curr_dir = os.path.abspath(curr_dir)
    if curr_dir[-1] != "/": curr_dir += "/"

    if os.path.exists(curr_dir + "06-trinity_assembly/") == False:
        os.mkdir(curr_dir + "06-trinity_assembly/")

    files_skip = []
    for filename in os.listdir(processed_dir):
        if filename in files_skip:
            continue
        file_comps = filename.split(".")
        file_sra = filename.split("_")[0]
        layout = process_reads.get_layout(file_sra)
        if layout == "SINGLE":
            cmd_trin = ["python3", "trinity_wrapper.py", processed_dir + filename, str(threads),
                        str(max_memory_gb), "non-stranded", curr_dir + "06-trinity_assembly/"]
            subprocess.run(" ".join(cmd_trin), shell = True)
        elif layout == "PAIRED":
            pe_file_1 = file_comps[0][0:-1:] + "1." + ".".join(file_comps[1::])
            pe_file_2 = file_comps[0][0:-1:] + "2." + ".".join(file_comps[1::])

            cmd_trin = ["python3", "trinity_wrapper.py", processed_dir + pe_file_1, processed_dir + pe_file_2,
                        str(threads), str(max_memory_gb), "non-stranded", curr_dir + "06-trinity_assembly/"]
            subprocess.run(" ".join(cmd_trin), shell = True)


if __name__ == "__main__":
    curr_dir = os.getcwd()
    if os.path.isabs(curr_dir) == False: curr_dir = os.path.abspath(curr_dir)
    if curr_dir[-1] != "/": curr_dir += "/"

    if len(sys.argv) == 3:
        assemble_trinity(curr_dir + "05-filter_over_represented/", sys.argv[1], sys.argv[2])
    else:
        print("Usage:")
        print("python3 assemble_reads.py threads max_memory_GB")
