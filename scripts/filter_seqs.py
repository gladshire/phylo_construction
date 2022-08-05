#!/bin/env python3


import sys
import os
import psutil
import shutil
import gzip
from time import sleep
import subprocess

from itertools import zip_longest
from os.path import exists

KRAKEN2_LOC = "~/miles/packages/kraken2/"
FILT_LIST = ["bacteria", "archaea", "fungi", "viral", "mitochondrion_and_plastid"]


def filter_summary(filt_list, in_file, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    print("\nKraken2 filtering complete!")
    print("SUMMARY:\n")
    total_reads_removed = 0
    total_percent_removed = 0
    for i in filt_list:
        curr_rep = open(out_dir + in_file.split("/")[-1].split(".")[0] + "." + i + ".report", 'r')
        lines = curr_rep.readlines()
        reads_removed = int(lines[1].strip().split()[1])
        percent_removed = float(lines[1].strip().split()[0])
        total_reads_removed += reads_removed
        total_percent_removed += percent_removed
        if i == filt_list[0]:
            total_reads = int(lines[0].strip().split()[1]) + int(lines[1].strip().split()[1])
        for line in lines:
            if "Bacteria" in line:
                print("  Bacterial reads: " + str(reads_removed) + " (" + str(percent_removed) + "%)")
                break
            elif "Fungi" in line:
                print("  Fungal reads: " + str(reads_removed) + " (" + str(percent_removed) + "%)")
                break
            elif "Archaea" in line:
                print("  Archaeal reads: " + str(reads_removed) + " (" + str(percent_removed) + "%)")
                break
            elif "Viruses" in line:
                print("  Viral reads removed: " + str(reads_removed) + " (" + str(percent_removed) + "%)")
                break
            elif "plant" in line:
                print("  Mitochondrial/Plastid reads: " + str(reads_removed) + " (" + str(percent_removed) + "%)")
                break
            else:
                continue

    print("\n  Total reads removed: " + str(total_reads_removed) + " (" + str(round(total_percent_removed, 2)) + "%)")
    print("  Leaving: " + str(total_reads - total_reads_removed) + " (" + str(round(100.0 - total_percent_removed, 2)) + "%) reads\n")
     
        
def gzip_filtered(out_dir, filtered_file, threads):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    file_name = out_dir + filtered_file.split("/")[-1].split(".")[0] + ".filt.fq"
    print("Compressing filtered file...")
    subprocess.run(["pigz", "-p" + str(threads), file_name], shell=True)
        

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(fillvalue = fillvalue, *args)


def kraken_filter_se(db, in_file, threads, out_dir, confidence=0.2):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    if exists(out_dir + in_file.split("/")[-1].split(".")[0] + ".filt.fq"):
        in_file = out_dir + in_file.split("/")[-1].split(".")[0] + ".filt.fq"
    out_file = "temp.fq"
    cmd = [KRAKEN2_LOC + "kraken2", "--threads", str(threads), "--unclassified-out",
           out_dir + out_file, "--db", db, in_file, "--output", "-", "--report",
           out_dir + in_file.split("/")[-1].split(".")[0] + "." + db.split("/")[-1] + ".report"]
    mem_available = psutil.virtual_memory()[1]
    
    print("Now filtering: " + db.split("/")[-1].replace("_", " "))
    if confidence != None:
        cmd.append("--confidence " + str(confidence))
    if os.path.getsize(db) < mem_available / 2:
        print("Not enough RAM for Kraken2 database {}".format(db))
        print("Switching to slower memory-mapping mode...")
        cmd.append("--memory-mapping")
    
    subprocess.run(" ".join(cmd), shell=True)
    os.rename(out_dir + "temp.fq", out_dir + in_file.split("/")[-1].split(".")[0] + ".filt.fq")
    

def kraken_filter_pe(db, in_file1, in_file2, threads, out_dir, confidence=0.2):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    if exists(out_dir + in_file1.split("/")[-1].split(".")[0] + ".filt.fq"):
        in_file = out_dir + in_file1.split("/")[-1].split(".")[0] + ".filt.fq"
    out_file = "temp.fq"
    cmd = [KRAKEN2_LOC + "kraken2", "--threads", str(threads), "--unclassified-out",
           out_dir + out_file, "--db", db, "--paired", in_file1, in_file2,
           "--output" "-", "--report", out_file.split("/")[-1].split(".")[0] + "." + db.split("/")[-1] + ".report"]
    mem_available = psutil.virtual_memory()[1]
    
    print("Now filtering: " + db.split("/")[-1].replace("_", " "))
    if confidence != None:
        cmd.append("--confidence " + str(confidence))
    if os.path.getsize(db) < memory_available / 2:
        print("Not enough RAM for Kraken2 database {}.".format(db))
        print("Switching to slower memory-mapping mode...")
        cmd.append("--memory-mapping")

    subprocess.run(" ".join(cmd), shell=True)
    os.rename(out_dir + "temp.fq", out_dir + in_file.split("/")[-1].split(".")[0] + ".filt.fq")




if __name__ == "__main__":
    if len(sys.argv) == 4:
        file_path, file_name = os.path.split(sys.argv[1])
        file_base = str(file_name).split(".")
        if os.path.exists(sys.argv[3] + "/" + file_base[0] + ".filt.fq"):
            print("Filtered file found for: " + sys.argv[1])
            sys.exit()
        for i in FILT_LIST:
            kraken_filter_se(db = "/scratch/kraken2_dbs/" + i,
                             in_file = sys.argv[1],
                             threads = sys.argv[2],
                             out_dir = sys.argv[3])
        filter_summary(FILT_LIST, in_file = sys.argv[1], out_dir = sys.argv[3])
    elif len(sys.argv) == 5:
        file_path_1, file_name_1 = os.path.split(sys.argv[1])
        file_path_2, file_path_2 = os.path.split(sys.argv[2])
        file_base_1 = str(file_name_1).split(".")
        file_base_2 = str(file_name_2).split(".")
        for i in FILT_LIST:
            if os.path.exists(sys.argv[4] + "/" + file_base_1[0] + ".filt.fq") and\
               os.path.exists(sys.argv[4] + "/" + file_base_2[0] + ".filt.fq"):
                print("Filtered files found for: ")
                print(sys.argv[1])
                print(sys.argv[2])
                sys.exit()
            kraken_filter_pe(db = "/scratch/kraken2_dbs/" + i,
                             in_file1 = sys.argv[1],
                             in_file2 = sys.argv[2],
                             threads = sys.argv[3],
                             out_dir = sys.argv[4])
        filter_summary(FILT_LIST, in_file = sys.argv[1], out_dir = sys.argv[3])
    else:
        print("Usage:")
        print("For single-end reads: python3 filter_seqs.py fq_file threads output_directory")
        print("For paired-end reads: python3 filter_seqs.py fq_file1 fq_file2 threads output_directory")
        print("python3 filter_seqs.py in_file threads out_dir")
