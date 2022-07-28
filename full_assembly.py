import sys
import os
import subprocess

sys.path.append(os.getcwd() + "/data_retrieval/")
sys.path.append(os.getcwd() + "/read_processing/scripts/")
sys.path.append(os.getcwd() + "/assembly/")
sys.path.append(os.getcwd() + "/filtering_translation/")

import dnld_sra
import read_process
import trinity_wrapper

curr_wd = os.getcwd() + "/"

if __name__ == "__main__":
   subprocess.run(" ".join(["python3", curr_wd + "data_retrieval/dnld_sra.py",
                  "8", curr_wd + "data_retrieval/", curr_wd + "read_processing/00-raw_data/"]), shell=True)

   subprocess.run(" ".join(["python3", curr_wd + "read_processing/scripts/read_process.py",
                  curr_wd + "read_processing/00-raw_data/SRR364263.fastq", "8"]), shell=True)

  
