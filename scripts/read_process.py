import os
import sys
import subprocess
import rcorrector_wrapper
import filter_unfixable
import trimmomatic_wrapper
import filter_seqs
import fastqc_wrapper
import filter_overrep

def read_process_se(se_fq, threads, out_dir=None, remove_inter=False):
    
    # Determine where to store output
    if out_dir != None:
        if out_dir == ".": out_dier = os.getcwd()
        if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
        if out_dir[-1] != "/": out_dir += "/"
    else:
        out_dir = os.path.abspath(os.path.dirname(__file__)) + "/"
   
    # If intermediate files not removed, store them in folders step-wise
    out_dir_inter = [None] * 5
    for d in out_dir_inter:
        d = ""
    if remove_inter == False:
        out_dir_inter[0] = "01-error_correction"
        out_dir_inter[1] = "02-filter_adapter_seq"
        out_dir_inter[2] = "03-filter_foreign_dna"
        out_dir_inter[3] = "04-quality_control"
        out_dir_inter[4] = "05-filter_over_represented"
        for d in out_dir_inter:
            os.mkdir(out_dir + d)
    os.mkdir(out_dir + "processed_reads")

    path_se, file_se = os.path.split(se_fq)
    se_fq_name = str(file_se)
    base_name_se = se_fq_name.split(".")[0]

    # Error correct with Rcorrector
    print("Running error correction...")
    rcorrector_wrapper.rcorrector_se(se_fq, threads, out_dir + out_dir_inter[0])
    corrected_se = out_dir + out_dir_inter[0] + "/" +  base_name_se + ".cor.fq"

    # Remove unfixable errors found by Rcorrector
    print("Removing unfixable error reads...")
    filter_unfixable.filter_unfix_se(corrected_se, out_dir + out_dir_inter[0])
    filtered_se = out_dir + out_dir_inter[0] + "/" + base_name_se + ".fix.fq"

    # Trim adapter sequences from reads with Trimmomatic
    print("Trimming adapter sequences...")
    trimmomatic_wrapper.trimmomatic_se(filtered_se, threads, out_dir + out_dir_inter[1])
    trimmed_se = out_dir + out_dir_inter[1] + "/" + base_name_se + ".trim.fq"

    # Filter undesired genome/transcriptome reads with Kraken2
    print("Filtering foreign reads...")
    subprocess.run(["python3", "filter_seqs.py", trimmed_se, threads, out_dir + "/" +  out_dir_inter[2]])
    foreign_filtered_se_reads = out_dir + out_dir_inter[2] + "/" + base_name_se + ".filt.fq"

    # Perform quality analysis with FastQC
    print("Running quality analysis...")
    fastqc_wrapper.fastQC_se(foreign_filtered_se_reads, threads, out_dir + out_dir_inter[3])
    se_fqc_path = out_dir + out_dir_inter[3] + "/" + base_name_se + ".filt_fastqc/fastqc_data.txt"

    # Remove over-represented reads found by FastQC
    print("Removing over-represented reads...")
    filter_overrep.filter_overrep_se(foreign_filtered_se_reads, se_fqc_path, out_dir + out_dir_inter[4])


def read_process_pe(pe_fq1, pe_fq2, threads, out_dir=None, remove_inter=False):
    
    # Determine where to store output
    if out_dir != None:
        if out_dir == ".": out_dier = os.getcwd()
        if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
        if out_dir[-1] != "/": out_dir += "/"
    else:
        out_dir = os.path.abspath(os.path.dirname(__file__)) + "/"


    # If intermediate files not removed, store them in folders step-wise
    out_dir_inter = [None] * 5
    for d in out_dir_inter:
        d = ""
    if remove_inter == False:
        out_dir_inter[0] = "01-error_correction"
        out_dir_inter[1] = "02-filter_adapter_seq"
        out_dir_inter[2] = "03-filter_foreign_dna"
        out_dir_inter[3] = "04-quality_control"
        out_dir_inter[4] = "05-filter_over_represented"
        for d in out_dir_inter:
            os.mkdir(out_dir + d)
    os.mkdir(out_dir + "processed_reads")

    path_pe_1, file_pe_1 = os.path.split(pe_fq1)
    pe_fq1_name = str(file_pe_1)
    base_name_pe_1 = pe_fq1_name.split(".")[0]

    path_pe_2, file_pe_2 = os.path.split(pe_fq2)
    pe_fq2_name = str(file_pe_2)
    base_name_pe_2 = pe_fq2_name.split(".")[0]

    # Error correct with Rcorrector
    print("Running error correction...")
    rcorrector_wrapper.rcorrector_pe(pe_fq1, pe_fq2, threads, out_dir + out_dir_inter[0])
    corrected_pe_1 = out_dir + out_dir_inter[0] + "/" + base_name_pe_1 + ".cor.fq"
    corrected_pe_2 = out_dir + out_dir_inter[0] + "/" + base_name_pe_2 + ".cor.fq"

    # Remove unfixable errors found by Rcorrector
    print("Removing unfixable error reads...")
    filter_unfixable.filter_unfix_pe(corrected_pe_1, corrected_pe_2, out_dir + out_dir_inter[0])
    filtered_pe_1 = out_dir + out_dir_inter[0] + "/" + base_name_pe_1 + ".fix.fq"
    filtered_pe_2 = out_dir + out_dir_inter[0] + "/" + base_name_pe_2 + ".fix.fq"

    # Trim adapter sequences from reads with Trimmomatic
    print("Trimming adapter sequences...")
    trimmomatic_wrapper.trimmomatic_pe(filtered_pe_1, filtered_pe_2, threads, out_dir + out_dir_inter[1])
    trimmed_pe_1 = out_dir + out_dir_inter[1] + "/" + base_name_pe_1 + ".trim.fq"
    trimmed_pe_2 = out_dir + out_dir_inter[1] + "/" + base_name_pe_2 + ".trim.fq"
    
    # Filter undesired genome/transcriptome reads with Kraken2
    print("Filtering foreign reads...")
    subprocess.run(["python3", "filter_seqs.py", trimmed_pe_1, trimmed_pe_2, threads, out_dir + out_dir_inter[2]])
    foreign_filtered_pe_1 = out_dir + out_dir_inter[2] + "/" + base_name_pe_1 + ".filt.fq"
    foreign_filtered_pe_2 = out_dir + out_dir_inter[2] + "/" + base_name_pe_2 + ".filt.fq"

    # Perform quality analysis with FastQC
    print("Running quality analysis...")
    fastqc_wrapper.fastQC_pe(foreign_filtered_pe_1, foreign_filtered_pe_2, threads, out_dir + out_dir_inter[3])
    pe_fqc1_path = out_dir + out_dir_inter[3] + "/" + base_name_pe_1 + ".filt_fastqc/fastqc_data.txt"
    pe_fqc2_path = out_dir + out_dir_inter[3] + "/" + base_name_pe_2 + ".filt_fastqc/fastqc_data.txt"

    # Remove over-represented reads found by FastQC
    print("Removing over-represented reads...")
    filter_overrep.filter_overrep_se(foreign_filtered_pe_1, foreign_filtered_pe_2, pe_fqc1_path, pe_fqc2_path, out_dir + out_dir_inter[4])



if __name__ == "__main__":
    if len(sys.argv) == 3:
        read_process_se(se_fq = sys.argv[1], threads = sys.argv[2])
    elif len(sys.argv) == 4:
        read_process_pe(pe_fq1 = sys.argv[1] , pe_fq2 = sys.argv[2] , threads = sys.argv[3])
    else:
        print("Usage:")
        print("For single-end reads: python3 read_process.py fastq_file threads [optional: out_directory, remove_inter=True]")
        print("For paired-end reads: python3 read_process.py fastq_file1 fastq_file2 threads [optional: out_directory, remove_inter=True]")

