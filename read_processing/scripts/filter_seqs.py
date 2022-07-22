#!/bin/env python3


import sys
import os
import psutil
import shutil
import gzip
from time import sleep

from itertools import zip_longest
from os.path import exists

KRAKEN2_LOC = "~/miles/packages/kraken2/"


def filter_summary(filt_list, in_file, out_dir):
    if out_dir == ".": out_dir = os.getcwd() + "/"
    print("\nKraken2 filtering complete!")
    print("SUMMARY:")
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
                print("    Bacterial reads: " + str(reads_removed) + " (" + str(percent_removed) + "%)")
                break
            elif "Fungi" in line:
                print("    Fungal reads: " + str(reads_removed) + " (" + str(percent_removed) + "%)")
                break
            elif "Archaea" in line:
                print("    Archaeal reads: " + str(reads_removed) + " (" + str(percent_removed) + "%)")
                break
            elif "Viruses" in line:
                print("    Viral reads removed: " + str(reads_removed) + " (" + str(percent_removed) + "%)")
                break
            elif "plant" in line:
                print("    Mitochondrial/Plastid reads: " + str(reads_removed) + " (" + str(percent_removed) + "%)")
                break
            else:
                continue

    print("\n    Total reads removed: " + str(total_reads_removed) + " (" + str(round(total_percent_removed, 2)) + "%)")
    print("    Leaving " + str(total_reads) + " (" + str(round(100.0 - total_percent_removed, 2)) + ") reads")
     
        
def gzip_filtered(out_dir, filtered_file, threads):
    if out_dir == ".": out_dir = os.getcwd() + "/"
    file_name = out_dir + filtered_file.split("/")[-1].split(".")[0] + ".filt.fq"
    print("Compressing filtered file...")
    os.system("pigz -p" + str(threads) + " " + file_name)
        

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(fillvalue = fillvalue, *args)


def kraken_filter_se(db, in_file, threads, out_dir, confidence=None):
    if out_dir == ".": out_dir = os.getcwd() + "/" 
    
    if exists(out_dir + in_file.split("/")[-1].split(".")[0] + ".filt.fq"):
        in_file = out_dir + in_file.split("/")[-1].split(".")[0] + ".filt.fq"
    #out_file = in_file.split("/")[-1].split(".")[0] + ".filt.fq"
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
    
    os.system(" ".join(cmd))
    os.rename(out_dir + "temp.fq", out_dir + in_file.split("/")[-1].split(".")[0] + ".filt.fq")
    

def kraken_filter_pe(db, in_file1, in_file2, threads, out_dir, confidence=None):
    if out_dir == ".": out_dir = os.getcwd() + "/"
    
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

    os.system(" ".join(cmd))

"""
def kraken_filter_se(filter_file, seq_file, out_dir):
    
    out_filt_name = seq_file.split(".")[0] + ".filt.fq"
    
    seq_in = open(seq_file, 'rb')
    seq_out = open(out_dir + out_filt_name, 'wb')
    
    filt_file = open(filter_file)
    iter_seqs = grouper(seq_in, 4)
    filt_ct = 0
    for line_num, line in enumerate(filt_file):
        
        print(iter_seqs[i])
        head,seq,placeholder,qual = [i.strip() for i in iter_seqs[line_num]]
        if line.split()[0] == "C":
            filt_ct += 1
        else:
            seq_out.write(b'%s\n' % b'\n'.join([head,seq,placeholder,qual]))

    seq_in.close()
    seq_out.close()
    filt_file.close()

    if seq_file.split(".")[-1] == "gz":
        os.system("gzip " + out_dir + out_filt_name)
    
    print("Kraken2 filter complete!")
    print("\nSummary:\n")
    print("Reads processed: " + str(line))
    print("Reads removed:   " + str(filt_ct))
"""

if __name__ == "__main__":
    filt_list = ["bacteria", "archaea", "fungi", "viral", "mitochondrion_and_plastid"]
    if len(sys.argv) == 4:
        for i in filt_list:
            kraken_filter_se(db = "/scratch/kraken2_dbs/" + i,
                             in_file = sys.argv[1],
                             threads = sys.argv[2],
                             out_dir = sys.argv[3])
        gzip_filtered(out_dir = sys.argv[3], filtered_file = sys.argv[1], threads = sys.argv[2])
        filter_summary(filt_list, in_file = sys.argv[1], out_dir = sys.argv[3])
    elif len(sys.argv) == 5:
        for i in filt_list:
            kraken_filter_pe(db = "/scratch/kraken2_dbs/" + i,
                             in_file = sys.argv[1],
                             threads = sys.argv[2],
                             out_dir = sys.argv[3])
        filter_summary(filt_list, in_file = sys.argv[1], out_dir = sys.argv[3])
    else:
        print("Usage:")
        print("python3 filter_seqs.py in_file threads out_dir")
