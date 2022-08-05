import sys
import os
import subprocess

def fastQC_pe(pe_fq1, pe_fq2, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_pe_1, file_pe_1 = os.path.split(pe_fq1)
    pe_fq1_name = str(file_pe_1)
    base_name_pe_1 = pe_fq1_name.split(".")[0] + ".org_filtered_fastqc.html"
   
    path_pe_2, file_pe_2 = os.path.split(pe_fq2)
    pe_fq2_name = str(file_pe_2)
    base_name_pe_2 = pe_fq2_name.split(".")[0] + ".org_filtered_fastqc.html"
    
    file_base = pe_fq1.split(".")
    if os.path.exists(out_dir + file_base[0][:-2:] + ".filt_fastqc"):
        print("FastQC file found for: ")
        print(pe_fq1)
        print(pe_fq2)
        return
    cmd = ["fastqc", pe_fq1, pe_fq2, "-t", str(threads), "-o", out_dir, "--extract"]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)

def fastQC_se(se_fq, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"
    
    path_se, file_se = os.path.split(se_fq)
    se_fq_name = str(file_se)
    base_name_se = se_fq_name.split(".")[0] + ".org_filtered_fastqc.html"
    base_file = se_fq_name.split(".")    

    if os.path.exists(out_dir + base_file[0] + ".filt_fastqc"):
        print("FastQC file found for: " + se_fq)
        return

    cmd = ["fastqc", se_fq, "-t", str(threads), "-o", out_dir, "--extract"]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)


if __name__ == "__main__":
    if len(sys.argv) == 4:
        fastQC_se(se_fq = sys.argv[1], threads = int(sys.argv[2]), out_dir = sys.argv[3])
    elif len(sys.argv) == 5:
        fastQC_pe(pe_fq1 = sys.argv[1], pe_fq2 = sys.argv[2], threads = int(sys.argv[3]),
                  out_dir = sys.argv[4])
    else:
        print("Usage:")
        print("For single-end reads: python3 fastqc_wrapper.py fastq_file threads output_directory")
        print("For paired-end reads: python3 fastqc_wrapper.py fastq_file1 fastq_file2 threads output_directory")
        sys.exit()
