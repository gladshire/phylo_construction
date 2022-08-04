import os
import sys
import subprocess

RCORRECTOR_LOC = "~/miles/packages/rcorrector"
RCORRECTOR_CMD = "perl " + RCORRECTOR_LOC + "/run_rcorrector.pl "

def rcorrector_se(se_fq, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_se, file_se = os.path.split(se_fq)
    se_file = str(file_se)
    corr_name = se_file.split(".")[0] + ".cor.fq"
   
    print(out_dir + corr_name) 
    if os.path.exists(out_dir + corr_name):
        print("Corrected file found for: " + se_fq)
        return
    cmd = [RCORRECTOR_CMD, "-s", se_fq, "-t", str(threads), "-od", out_dir]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)

def rcorrector_pe(pe_fq1, pe_fq2, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_pe_1, file_pe_1 = os.path.split(pe_fq1)
    path_pe_2, file_pe_2 = os.path.split(pe_fq2)

    pe_file_1 = str(file_pe_1)
    corr_name_1 = pe_file_1.split(".")[0] + ".cor.fq"

    pe_file_2 = str(file_pe_2)
    corr_name_2 = pe_file_2.split(".")[0] + ".cor.fq"

    if os.path.exists(out_dir + corr_name_1) and \
       os.path.exists(out_dir + corr_name_2):
        print("Corrected files found for: ")
        print(pe_fq1)
        print(pe_fq2)
        return
    cmd = [RCORRECTOR_CMD, "-1", pe_fq1, "-2", pe_fq2, "-t", str(threads), "-od", out_dir]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)


if __name__ == "__main__":
    if len(sys.argv) == 4:
        rcorrector_se(se_fq = sys.argv[1], threads = int(sys.argv[2]), out_dir = sys.argv[3])
    elif len(sys.argv) == 5:
        rcorrector_pe(pe_fq1 = sys.argv[1], pe_fq2 = sys.argv[2], threads = int(sys.argv[3]),
                      out_dir = sys.argv[4])
    else:
        print("Usage:")
        print("For single end reads: python3 rcorrector.py fastq_se_reads threads output_directory")
        print("For paired end reads: python3 rcorrector.py fastq_pe_reads1 fastq_pe_reads2 threads output_directory")
