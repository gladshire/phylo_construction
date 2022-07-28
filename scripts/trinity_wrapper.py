import os
import sys
import subprocess

TRINITY_LOC = "~/miles/packages/trinityrnaseq-v2.14.0/"
TRINITY_CMD = TRINITY_LOC + "Trinity"

def trinity_pe(pe_fq1, pe_fq2, taxon_id, threads, max_memory_gb, strand, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    assembly_file = taxon_id + ".Trinity.fasta"
    
    if strand == "stranded":
        stranded = "--SS_lib_type RF"
    else:
        stranded = ""
   
    #cmd = "ulimit -s unlimited\n"
    cmd = [TRINITY_CMD, "--seqType fq", "--left", pe_fq1, "--right", pe_fq2,
           "--max_memory", str(max_memory_gb) + "G", "--CPU", str(threads),
           "--bflyCalculateCPU", "--full_cleanup", "--no_normalize_reads",
           stranded, "--output", out_dir + taxon_id + ".trinity"]
    
    subprocess.run(" ".join(cmd), shell=True)
    os.rename(out_dir + taxon_id + ".trinity.Trinity.fasta", out_dir + assembly_file)

def trinity_se(se_fq, taxon_id, threads, max_memory_gb, strand, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != ".": out_dir += "/"

    assembly_file = taxon_id + ".Trinity.fasta"

    if strand == "stranded":
        stranded = "--SS_lib_type R"
    else:
        stranded = ""

    cmd = [TRINITY_CMD, "--seqType", "fq", "--single", se_fq, "--max_memory",
           str(max_memory_gb) + "G", "--CPU", str(threads), "--bflyCalculateCPU",
           "--full_cleanup", "--no_normalize_reads", stranded, "--output",
           out_dir + str(taxon_id) + ".trinity"]
    subprocess.run(" ".join(cmd), shell=True)
    #os.rename(out_dir + taxon_id + ".trinity.Trinity.fasta", out_dir + assembly_file)


if __name__ == "__main__":
    if len(sys.argv) == 8:
        trinity_pe(pe_fq1 = sys.argv[1], pe_fq2 = sys.argv[2],
                   taxon_id = sys.argv[3], threads = sys.argv[4],
                   max_memory_gb = sys.argv[5], strand = sys.argv[6],
                   out_dir = sys.argv[7])
    elif len(sys.argv) == 7:
        trinity_se(se_fq = sys.argv[1], taxon_id = sys.argv[2],
                   threads = sys.argv[3], max_memory_gb = sys.argv[4],
                   strand = sys.argv[5], out_dir = sys.argv[6])
    else:
        print("Usage:")
        print("For single-end reads: python3 trinity_wrapper.py fastq_file taxonID threads max_memory_GB strand(stranded / non-stranded) output_directory")
        print("For paired-end reads: python3 trinity_wrapper.py fastq_file1 fastq_file2 taxonID threads max_memory_GB strand(stranded / non-stranded) output_directory")
        sys.exit()
