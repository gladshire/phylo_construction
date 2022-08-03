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

    for filename in os.listdir(processed_dir):
        file_sra = filename.split("_")[0]
        layout = process_reads.get_layout(file_sra)
        if layout == "SINGLE":
            cmd_trin = ["python3", "trinity_wrapper.py", processed_dir + filename, str(threads),
                        str(max_memory_gb), "non-stranded", curr_dir + "06-trinity_assembly/"]
            subprocess.run(" ".join(cmd_trin), shell = True)
        elif layout == "PAIRED":
