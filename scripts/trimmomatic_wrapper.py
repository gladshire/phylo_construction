import sys
import os
import subprocess

TRIMMOMATIC_LOC = "~/miles/packages/Trimmomatic-0.39"
TRIMMOMATIC_CMD = "java -jar " + TRIMMOMATIC_LOC + "/dist/jar/trimmomatic-0.39.jar"
TRUSEQ_ADAPTER = os.path.expanduser("~/miles/packages/phylogenomic_dataset_construction/databases/trimmomatic_adapters.fasta")

def trimmomatic_pe(pe_fq1, pe_fq2, threads, out_dir):
    trim_setting = "ILLUMINACLIP:" + TRUSEQ_ADAPTER + ":2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"

    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_pe_1, file_pe_1 = os.path.split(pe_fq1)
    pe_fq1_name = str(file_pe_1)
    base_name_pe_1 = pe_fq1_name.split(".")[0] + ".paired.trim.fq"
    base_name_up_1 = pe_fq1_name.split(".")[0] + ".unpaired.trim.fq"

    path_pe_2, file_pe_2 = os.path.split(pe_fq2)
    pe_fq2_name = str(file_pe_2)
    base_name_pe_2 = pe_fq2_name.split(".")[0] + ".paired.trim.fq"
    base_name_up_2 = pe_fq1_name.split(".")[0] + ".unpaired.trim.fq"

    if os.path.exists(out_dir + base_name_pe_1) and\
       os.path.exists(out_dir + base_name_up_1) and\
       os.path.exists(out_dir + base_name_pe_2) and\
       os.path.exists(out_dir + base_name_up_2):
        print("Trimmed files found for: ")
        print(os.path.split(pe_fq1)[-1])
        print(os.path.split(pe_fq2)[-1])
        return

    cmd = [TRIMMOMATIC_CMD, "PE", "-threads", str(threads),
           pe_fq1, pe_fq2,
           out_dir + base_name_pe_1,
           out_dir + base_name_up_1,
           out_dir + base_name_pe_2,
           out_dir + base_name_up_2,
           trim_setting]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)

def trimmomatic_se(se_fq, threads, out_dir):
    #trim_setting = "ILLUMINACLIP:" + TRUSEQ_ADAPTER + ":2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"


    trim_setting = "ILLUMINACLIP:" + TRUSEQ_ADAPTER + ":2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:25"
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_se, file_se = os.path.split(se_fq)
    se_fq_name = str(file_se)
    base_name_se = se_fq_name.split(".")[0] + ".trim.fq"

    if os.path.exists(out_dir + base_name_se):
        print("Trimmed file found for: " + os.path.split(se_fq)[-1])
        return

    cmd = [TRIMMOMATIC_CMD, "SE", "-threads", str(threads),
           se_fq,
           out_dir + base_name_se,
           trim_setting]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)


if __name__ == "__main__":
    if len(sys.argv) == 4:
        trimmomatic_se(se_fq = sys.argv[1], threads = int(sys.argv[2]), out_dir = sys.argv[3])
    elif len(sys.argv) == 5:
        trimmomatic_pe(pe_fq1 = sys.argv[1], pe_fq2 = sys.argv[2], threads = int(sys.argv[3]),
                       out_dir = sys.argv[4])
    else:
        print("Usage:")
        print("For single-end reads: python3 trimmomatic_wrapper.py fastq_file threads output_directory")
        print("For paired-end reads: python3 trimmomatic_wrapper.py fastq_file1 fastq_file2 threads output_directory")
        sys.exit()
