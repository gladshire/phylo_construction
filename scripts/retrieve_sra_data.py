import os
import sys
import re
import subprocess
from Bio import Entrez

Entrez.email = "mwoodc2@uic.edu"
Entrez.api_key = "6512380b5cb6a15be6ef4919564e2e186508"

SRA_TOOLKIT_BIN_PATH = "~/miles/packages/sratoolkit.3.0.0-ubuntu64/bin/"
SRA_LIST_PATH = "~/miles/allele_phasing/phylo_construction/data_retrieval/sras.txt"


def get_sra_numbers(sra_file_path):
    sra_list_file = open("sras.txt", 'r')
    ln = sra_list_file.readline()
    ln = sra_list_file.readline()
    sra_list = []
    while ln:
        sra_list.append(ln[:-1:])
        ln = sra_list_file.readline()
    return sra_list


def get_tax_id(sra_number):
    handle = Entrez.efetch(db = "sra", id = sra_number)
    handle_str = str(handle.read())
    tax_id_fld_pat = re.compile('tax_id="([^"]*)"')
    tax_id_fld_mat = tax_id_fld_pat.search(handle_str)
    tax_id_fld_sub = handle_str[tax_id_fld_mat.start():tax_id_fld_mat.end()]
    tax_id_pat = re.compile('"([^"]*)"')
    tax_id_mat = tax_id_pat.search(tax_id_fld_sub)
    taxon_id = tax_id_fld_sub[tax_id_mat.start():tax_id_mat.end()][1:-1:]
    return taxon_id


def get_org(sra_number):
    handle = Entrez.efetch(db = "sra", id = sra_number)
    handle_str = str(handle.read())
    org_fld_pat = re.compile('organism="([^"]*)"')
    org_fld_mat = org_fld_pat.search(handle_str)
    org_fld_sub = handle_str[org_fld_mat.start():org_fld_mat.end()]
    org_pat = re.compile('"([^"]*)"')
    org_mat = org_pat.search(org_fld_sub)
    org = org_fld_sub[org_mat.start():org_mat.end()][1:-1:]
    return org


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
        taxon_id = get_tax_id(sra)
        organism = get_org(sra)
        out_file = taxon_id + "_" + organism.lower().replace(" ", "_") + ".fastq"

        print("Downloading: " + organism + " (" + taxon_id + ")" + " ...")
        prefetch_cmd = [SRA_TOOLKIT_BIN_PATH + "prefetch", sra,
                        "--max-size", "100g", "--progress",
                        "--output-directory", out_dir_prefetch]
        subprocess.run(" ".join(prefetch_cmd), shell=True)

        print("Dumping " + sra + ".sra to fastq ...")
        fasterq_cmd = [SRA_TOOLKIT_BIN_PATH + "fasterq-dump",
                       out_dir_prefetch + sra, "--threads",
                       str(threads), "--outfile",
                       out_dir_fastq + out_file]
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
