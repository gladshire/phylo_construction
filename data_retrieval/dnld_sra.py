import os
import sys
import subprocess

SRA_TOOLKIT_BIN_PATH = "~/miles/packages/sratoolkit.3.0.0-ubuntu64/bin/"
SRA_LIST_PATH = "~/miles/allele_phasing/phylo_construction/data_retrieval/sras.txt"


def get_sra_numbers(sra_file_path):
    sra_list_file = open(os.getcwd() + "/sras.txt", 'r')
    ln = sra_list_file.readline()
    ln = sra_list_file.readline()
    sra_list = []
    while ln:
        sra_list.append(ln[:-1:])
        ln = sra_list_file.readline()
    return sra_list


def get_sra_data(sra_list, threads, out_dir_prefetch, out_dir_fastq):
    if out_dir_prefetch == ".": out_dir_prefetch = os.getcwd()
    if os.path.isabs(out_dir_prefetch) == False:
        out_dir_prefetch = os.path.abspath(out_dir_prefetch)
    if out_dir_prefetch[-1] != "/": out_dir_prefetch += "/"

    if out_dir_fastq == ".": out_dir_fastq = os.getcwd()
    if os.path.isabs(out_dir_fastq) == False:
        out_dir_fastq = os.path.abspath(out_dir_fastq)
    if out_dir_fastq[-1] != "/": out_dir_fastq += "/"

    for sra in sra_list:
        print("Downloading: " + sra + "to " + out_dir_prefetch)
        prefetch_cmd = [SRA_TOOLKIT_BIN_PATH + "prefetch", sra,
                        "--max-size", "100g", "--progress",
                        "--output-directory", out_dir_prefetch]
        subprocess.run(" ".join(prefetch_cmd), shell=True)
        print("Dumping " + sra + ".sra to fastq ...")
        fasterq_cmd = [SRA_TOOLKIT_BIN_PATH + "fasterq-dump",
                       out_dir_prefetch + sra, "--threads",
                       str(threads), "--outdir", out_dir_fastq]
        subprocess.run(" ".join(fasterq_cmd), shell=True)
    print("Done.")    


if __name__ == "__main__":
    if len(sys.argv) == 4:
        sra_list = get_sra_numbers(SRA_LIST_PATH)
        get_sra_data(sra_list, threads = sys.argv[1], out_dir_prefetch = sys.argv[2], out_dir_fastq = sys.argv[3])
    else:
        print("Usage:")
        print("Enter SRA numbers in sras.txt prior to running this script")
        print("To download SRA data: python3 dnld_sra.py threads output_directory_prefetch output_directory_fastq")
