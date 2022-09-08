import os
import sys
import shutil
import re
import subprocess
import read_process
from Bio import Entrez

Entrez.email = "mwoodc2@uic.edu"
Entrez.api_key = "6512380b5cb6a15be6ef4919564e2e186508"


def dnld_sra_data(threads, out_dir):

    if os.path.exists(out_dir + "00-raw_seq_data/") == False:
        os.mkdir(out_dir + "00-raw_seq_data/")
    if os.path.exists(out_dir + "00-raw_seq_data/prefetch_repo/") == False:
        os.mkdir(out_dir + "00-raw_seq_data/prefetch_repo/")
    if os.path.exists(out_dir + "00-raw_seq_data/fasterq_dump/") == False:
        os.mkdir(out_dir + "00-raw_seq_data/fasterq_dump/")

    cmd_retrieve = ["python3", "retrieve_sra_data.py", str(threads),
                    out_dir + "00-raw_seq_data/prefetch_repo/",
                    out_dir + "00-raw_seq_data/fasterq_dump/"]

    subprocess.run(" ".join(cmd_retrieve), shell = True)


def get_layout(file_sra):
    handle = Entrez.efetch(db = "sra", id = file_sra)
    handle_str = str(handle.read())
    layout_fld_pat = re.compile("<LIBRARY_LAYOUT><.*/></LIBRARY_LAYOUT>")
    layout_fld_mat = layout_fld_pat.search(handle_str)
    layout_fld_sub = handle_str[layout_fld_mat.start():layout_fld_mat.end()]
    layout_pat = re.compile("><.*/><")
    layout_mat = layout_pat.search(layout_fld_sub)
    layout = layout_fld_sub[layout_mat.start():layout_mat.end() - 1]
    if "SINGLE" in layout:
        return "SINGLE"
    elif "PAIRED" in layout:
        return "PAIRED"
    else:
        return None


def process_fastq_dumps(fastq_dir, out_dir, threads, remove_inter):
    se_files = []
    pe_files_1 = []
    pe_files_2 = []

    for filename in os.listdir(fastq_dir):
        file_sra = filename.split("_")[0]
        layout = get_layout(file_sra)
        if layout == "SINGLE":
            se_files.append(filename)
        elif layout == "PAIRED":
            file_base = filename.split(".")[0]
            if file_base[-1] == "1":
                pe_files_1.append(filename)
            elif file_base[-1] == "2":
                pe_files_2.append(filename)
            else:
                print("Error: Unusual naming scheme for paired reads")
                sys.exit()

    print("\nNow initiating processing of single-ended runs ...\n")
    read_process.read_process_se(se_fq_files = [fastq_dir + se for se in se_files],
                                 threads = str(threads),
                                 remove_inter = remove_inter)
    print("\nNow initiating processing of paired-ended runs ...\n")
    read_process.read_process_pe(pe_fq1_files = [fastq_dir + pe_1 for pe_1 in pe_files_1],
                                 pe_fq2_files = [fastq_dir + pe_2 for pe_2 in pe_files_2],
                                 threads = str(threads),
                                 remove_inter = remove_inter)

if __name__ == "__main__":
    if len(sys.argv) == 3:
        threads = sys.argv[1]
        out_dir = os.getcwd()
        rm_inter = sys.argv[2]
        if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
        if out_dir[-1] != "/": out_dir += "/"

        # Download all SRA fastq files
        # dnld_sra_data(threads, out_dir)

        # Perform read processing on all SRA fastq files
        if rm_inter == "rem-inter":
            process_fastq_dumps(out_dir + "00-raw_seq_data/fasterq_dump/", out_dir, threads, remove_inter = True)
        elif rm_inter == "keep-inter":
            process_fastq_dumps(out_dir + "00-raw_seq_data/fasterq_dump/", out_dir, threads, remove_inter = False)
 
    else:
        print("Usage:")
        print("python3 process_reads.py threads remove_intermediates[rem-inter/keep-inter]")
        sys.exit()  

