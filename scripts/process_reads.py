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

    if os.path.exists(out_dir + "00-raw_seq_data/"):
        shutil.rmtree(out_dir + "00-raw_seq_data/")
    os.mkdir(out_dir + "00-raw_seq_data/")
    os.mkdir(out_dir + "00-raw_seq_data/prefetch_repo/")
    os.mkdir(out_dir + "00-raw_seq_data/fasterq_dump/")

    cmd_retrieve = ["python3", "retrieve_sra_data.py", str(threads),
                    out_dir + "00-raw_seq_data/prefetch_repo/",
                    out_dir + "00-raw_seq_data/fasterq_dump/"]

    subprocess.run(" ".join(cmd_retrieve), shell = True)


def get_layout(file_sra):
    handle = Entrez.efetch(db = "sra", id = file_sra)
    handle_str = str(handle.read())
    layout_fld_pat = re.compile("<LIBRARY_LAYOUT><[a-zA-Z]+\/>")
    layout_fld_mat = layout_fld_pat.search(handle_str)
    layout_fld_sub = handle_str[layout_fld_mat.start():layout_fld_mat.end()]
    layout_pat = re.compile("[a-zA-Z]+\/")
    layout_mat = layout_pat.search(layout_fld_sub)
    layout = layout_fld_sub[layout_mat.start():layout_mat.end() - 1]
    return layout


def process_fastq_dumps(fastq_dir, threads):
    for filename in os.listdir(fastq_dir):
        file_sra = filename.split("_")[0]
        layout = get_layout(file_sra)
        if layout == "SINGLE":
            read_process.read_process_se(fastq_dir + filename, str(threads))
        elif layout == "PAIRED":
            filename1 = fastq_dir + filename[:-1:] + "1"
            filename2 = fastq_dir + filename[:-1:] + "2"
            read_process.read_process_pe(filename1, filename2, str(threads))
        else:
            print("This should not happen")
            sys.exit()


if __name__ == "__main__":
    if len(sys.argv) == 2:
        threads = sys.argv[1]
        out_dir = os.getcwd()
        if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
        if out_dir[-1] != "/": out_dir += "/"

        # Download all SRA fastq files
        dnld_sra_data(threads, out_dir)

        # Perform read processing on all SRA fastq files
        process_fastq_dumps(out_dir + "00-raw_seq_data/fasterq_dump/", threads)
 
    else:
        print("Usage:")
        print("python3 process_reads.py threads")
        sys.exit()  

