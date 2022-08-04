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
            if os.path.exists(curr_dir + "06-trinity_assembly/" + file_comps[0] + ".Trinity.fasta"):
                print("Trinity assembly found for: " + filename)
                continue 
            cmd_trin = ["python3", "trinity_wrapper.py", processed_dir + filename, str(threads),
                        str(max_memory_gb), "non-stranded", curr_dir + "06-trinity_assembly/"]
            subprocess.run(" ".join(cmd_trin), shell = True)
        elif layout == "PAIRED":
            pe_file_1 = file_comps[0][0:-1:] + "1" + ".".join(file_comps[1::])
            pe_file_2 = file_comps[0][0:-1:] + "2" + ".".join(file_comps[1::])
            if os.path.exists(curr_dir + "06-trinity_assembly/" + pe_file_1) and\
               os.path.exists(curr_dir + "06-trinity_assembly/" + pe_file_2):
                print("Trinity assemblies found for: ")
                print(pe_file_1)
                print(pe_file_2)
                continue
            files_skip.append(pe_file_1)
            files_skip.append(pe_file_2)

            cmd_trin = ["python3", "trinity_wrapper.py", processed_dir + pe_file_1, processed_dir + pe_file_2, str(threads), str(max_memory_gb), "non-stranded", curr_dir + "06-trinity_assembly/"]
            subprocess.run(" ".join(cmd_trin), shell = True)
