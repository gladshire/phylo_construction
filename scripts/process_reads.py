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


def process_fastq_dumps(fastq_dir, out_dir, threads):
    file_skip = []
    out_dir_inter = [None] * 5 
    out_dir_inter[0] = "01-error_correction"
    out_dir_inter[1] = "02-filter_adapter_seq"
    out_dir_inter[2] = "03-filter_foreign_dna"
    out_dir_inter[3] = "04-quality_control"
    out_dir_inter[4] = "05-filter_over_represented"
    for d in out_dir_inter:
        if os.path.exists(d):
            shutil.rmtree(d)
        os.mkdir(out_dir + d)

    for filename in os.listdir(fastq_dir):
        if filename in file_skip:
            continue
        file_sra = filename.split("_")[0]
        layout = get_layout(file_sra)
        if layout == "SINGLE":
            continue
            read_process.read_process_se(fastq_dir + filename, str(threads))
        elif layout == "PAIRED":
            file_base = filename.split(".")[0]
            filename1 = fastq_dir + file_base[:-1:] + "1" + ".fastq"
            filename2 = fastq_dir + file_base[:-1:] + "2" + ".fastq"
            if file_base[-1] == "1":
                file_skip.append(file_base[:-1:] + "2" + ".fastq")
            elif file_base[-1] == "2":
                file_skip.append(file_base[:-1:] + "1" + ".fastq")
            else:
                print("This should not happen")
                sys.exit()
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
        process_fastq_dumps(out_dir + "00-raw_seq_data/fasterq_dump/", out_dir, threads)
 
    else:
        print("Usage:")
        print("python3 process_reads.py threads")
        sys.exit()  

