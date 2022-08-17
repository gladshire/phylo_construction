import os
import sys
import subprocess
import shutil
import rcorrector_wrapper
import filter_unfixable
import trimmomatic_wrapper
import filter_seqs
import fastqc_wrapper
import filter_overrep

def read_process_se(se_fq_files, threads, out_dir=None, remove_inter=False):
    # Determine where to store output
    if out_dir != None:
        if out_dir == ".": out_dir = os.getcwd()
        if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
        if out_dir[-1] != "/": out_dir += "/"
    else:
        out_dir = os.path.abspath(os.path.dirname(__file__)) + "/"

    out_dir_inter = [None] * 5
    for d in out_dir_inter:
        d = ""
    if remove_inter == False:
        out_dir_inter[0] = "01-error_correction"
        out_dir_inter[1] = "02-filter_adapter_seq"
        out_dir_inter[2] = "03-filter_foreign_dna"
        out_dir_inter[3] = "04-quality_control"
        out_dir_inter[4] = "05-filter_over_represented"
    
    paths_se = []
    files_se = []
    se_fq_names = []
    base_names_se = []

    for filename in se_fq_files:
        path_file_tup = os.path.split(filename)
        paths_se.append(path_file_tup[0])
        files_se.append(path_file_tup[1])
        se_fq_names.append(str(path_file_tup[1]))
        base_names_se.append(se_fq_names[-1].split(".")[0])

    # Error correct with Rcorrector
    print("\nRunning error correction ...\n")
    if os.path.exists(out_dir + out_dir_inter[0]) == False:
        os.mkdir(out_dir + out_dir_inter[0])
        
    corrected_ses = []
    for file_num, file_name in enumerate(se_fq_files):
        print("Correcting: " + file_name)
        rcorrector_wrapper.rcorrector_se(file_name, threads, out_dir + out_dir_inter[0])
        corrected_ses.append(out_dir + out_dir_inter[0] + "/" +  base_names_se[file_num] + ".cor.fq")

    # Remove unfixable errors found by Rcorrector
    print("\nRemoving unfixable error reads ...\n")
    filtered_ses = []
    for file_num, file_cor in enumerate(corrected_ses):
        print("Filtering: " + file_cor)
        filter_unfixable.filter_unfix_se(file_cor, out_dir + out_dir_inter[0])
        filtered_ses.append(out_dir + out_dir_inter[0] + "/" + base_names_se[file_num] + ".fix.fq")

    # Trim adapter sequences from reads with Trimmomatic
    print("\nTrimming adapter sequences ...\n")
    if os.path.exists(out_dir + out_dir_inter[1]) == False:
        os.mkdir(out_dir + out_dir_inter[1])

    trimmed_ses = []
    for file_num, file_filt in enumerate(filtered_ses):
        print("Trimming: " + file_filt)
        trimmomatic_wrapper.trimmomatic_se(file_filt, threads, out_dir + out_dir_inter[1])
        trimmed_ses.append(out_dir + out_dir_inter[1] + "/" + base_names_se[file_num] + ".trim.fq")

    # Filter undesired genome/transcriptome reads with Kraken2
    print("\nFiltering foreign reads ...\n")
    if os.path.exists(out_dir + out_dir_inter[2]) == False:
        os.mkdir(out_dir + out_dir_inter[2])

    foreign_filt_se_reads = []
    for file_num, file_trim in enumerate(trimmed_ses):
        print("Filtering: " + file_trim)
        subprocess.run(["python3", "filter_seqs.py", file_trim, threads, out_dir + "/" +  out_dir_inter[2]])
        foreign_filt_se_reads.append(out_dir + out_dir_inter[2] + "/" + base_names_se[file_num] + ".filt.fq")

    # Perform quality analysis with FastQC
    print("\nRunning quality analysis ...\n")
    if os.path.exists(out_dir + out_dir_inter[3]) == False:
        os.mkdir(out_dir + out_dir_inter[3])

    se_fqc_paths = []
    for file_num, file_for_filt in enumerate(foreign_filt_se_reads):
        print("Analyzing: " + file_for_filt)
        fastqc_wrapper.fastQC_se(file_for_filt, threads, out_dir + out_dir_inter[3])
        se_fqc_paths.append(out_dir + out_dir_inter[3] + "/" + base_names_se[file_num] + ".filt_fastqc/fastqc_data.txt")

    # Remove over-represented reads found by FastQC
    print("\nRemoving over-represented reads ...\n")
    if os.path.exists(out_dir + out_dir_inter[4]):
        os.mkdir(out_dir + out_dir_inter[4])

    for file_num, file_for_filt in enumerate(foreign_filt_se_reads):
        print("Filtering: " + file_for_filt)
        filter_overrep.filter_overrep_se(file_for_filt, se_fqc_paths[file_num], out_dir + out_dir_inter[4])

def read_process_pe(pe_fq1_files, pe_fq2_files, threads, out_dir=None, remove_inter=False):
    # Determine where to store output
    pe_fq1_files.sort()
    pe_fq2_files.sort()
    if out_dir != None:
        if out_dir == ".": out_dier = os.getcwd()
        if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
        if out_dir[-1] != "/": out_dir += "/"
    else:
        out_dir = os.path.abspath(os.path.dirname(__file__)) + "/"

    out_dir_inter = [None] * 5
    for d in out_dir_inter:
        d = ""
    if remove_inter == False:
        out_dir_inter[0] = "01-error_correction"
        out_dir_inter[1] = "02-filter_adapter_seq"
        out_dir_inter[2] = "03-filter_foreign_dna"
        out_dir_inter[3] = "04-quality_control"
        out_dir_inter[4] = "05-filter_over_represented"
    
    paths_pe_1 = []
    paths_pe_2 = []

    files_pe_1 = []
    files_pe_2 = []

    pe_fq_names_1 = []
    pe_fq_names_2 = []

    base_names_pe_1 = []
    base_names_pe_2 = []

    for filename in pe_fq1_files:
        path_file_tup = os.path.split(filename)
        paths_pe_1.append(path_file_tup[0])
        files_pe_1.append(path_file_tup[1])
        pe_fq_names_1.append(str(path_file_tup[1]))
        base_names_pe_1.append(pe_fq_names_1[-1].split(".")[0])

    for filename in pe_fq2_files:
        path_file_tup = os.path.split(filename)
        paths_pe_2.append(path_file_tup[0])
        files_pe_2.append(path_file_tup[1])
        pe_fq_names_2.append(str(path_file_tup[1]))
        base_names_pe_2.append(pe_fq_names_2[-1].split(".")[0])

    num_pes = len(pe_fq_names_1)

    # Error correct with Rcorrector
    print("\nRunning error correction ...\n")
    if os.path.exists(out_dir + out_dir_inter[0]):
        os.mkdir(out_dir + out_dir_inter[0])

    corrected_pes_1 = []
    corrected_pes_2 = []
    for i in range(0, num_pes):
        print("Correcting: ")
        print(pe_fq1_files[i])
        print(pe_fq2_files[i])
        rcorrector_wrapper.rcorrector_pe(pe_fq1_files[i], pe_fq2_files[i], threads, out_dir + out_dir_inter[0])
        corrected_pes_1.append(out_dir + out_dir_inter[0] + "/" + base_names_pe_1[i] + ".cor.fq")
        corrected_pes_2.append(out_dir + out_dir_inter[0] + "/" + base_names_pe_2[i] + ".cor.fq")

    # Remove unfixable errors found by Rcorrector
    print("\nRemoving unfixable error reads ...\n")
    filtered_pes_1 = []
    filtered_pes_2 = []
    for i in range(0, num_pes):
        print("Filtering: ")
        print(corrected_pes_1[i])
        print(corrected_pes_2[i])
        filter_unfixable.filter_unfix_pe(corrected_pes_1[i], corrected_pes_2[i], out_dir + out_dir_inter[0])
        filtered_pes_1.append(out_dir + out_dir_inter[0] + "/" + base_names_pe_1[i] + ".fix.fq")
        filtered_pes_2.append(out_dir + out_dir_inter[0] + "/" + base_names_pe_2[i] + ".fix.fq")

    # Trim adapter sequences from reads with Trimmomatic
    print("\nTrimming adapter sequences ...\n")
    if os.path.exists(out_dir + out_dir_inter[1]):
        os.mkdir(out_dir + out_dir_inter[1])

    trimmed_pes_1 = []
    trimmed_pes_2 = []
    for i in range(0, num_pes):
        print("Trimming:")
        print(filtered_pes_1[i])
        print(filtered_pes_2[i])
        trimmomatic_wrapper.trimmomatic_pe(filtered_pes_1[i], filtered_pes_2[i], threads, out_dir + out_dir_inter[1])
        trimmed_pes_1.append(out_dir + out_dir_inter[1] + "/" + base_names_pe_1[i] + ".paired.trim.fq")
        trimmed_pes_2.append(out_dir + out_dir_inter[1] + "/" + base_names_pe_2[i] + ".paired.trim.fq")
    
    # Filter undesired genome/transcriptome reads with Kraken2
    print("\nFiltering foreign reads ...\n")
    if os.path.exists(out_dir + out_dir_inter[2]):
        os.mkdir(out_dir + out_dir_inter[2])

    foreign_filt_pes_1 = []
    foreign_filt_pes_2 = []
    for i in range(0, num_pes):
        print("Filtering:")
        print(trimmed_pes_1[i])
        print(trimmed_pes_2[i])
        subprocess.run(["python3", "filter_seqs.py", trimmed_pes_1[i], trimmed_pes_2[i], threads, out_dir + out_dir_inter[2]])
        foreign_filt_pes_1.append(out_dir + out_dir_inter[2] + "/" + base_names_pe_1[i] + ".filt.fq")
        foreign_filt_pes_2.append(out_dir + out_dir_inter[2] + "/" + base_names_pe_2[i] + ".filt.fq")

    # Perform quality analysis with FastQC
    print("\nRunning quality analysis ...\n")
    if os.path.exists(out_dir + out_dir_inter[3]):
        os.mkdir(out_dir + out_dir_inter[3])

    pe_fqc1_paths = []
    pe_fqc2_paths = []
    for i in range(0, num_pes):
        print("Analyzing:")
        print(foreign_filt_pes_1[i])
        print(foreign_filt_pes_2[i])
        fastqc_wrapper.fastQC_pe(foreign_filt_pes_1[i], foreign_filt_pes_2[i], threads, out_dir + out_dir_inter[3])
        pe_fqc1_paths.append(out_dir + out_dir_inter[3] + "/" + base_names_pe_1[i] + ".filt_fastqc/fastqc_data.txt")
        pe_fqc2_paths.append(out_dir + out_dir_inter[3] + "/" + base_names_pe_2[i] + ".filt_fastqc/fastqc_data.txt")

    # Remove over-represented reads found by FastQC
    print("\nRemoving over-represented reads ...\n")
    if os.path.exists(out_dir + out_dir_inter[4]):
        os.mkdir(out_dir + out_dir_inter[4])

    for i in range(0, num_pes):
        print("Filtering:")
        print(foreign_filt_pes_1[i])
        print(foreign_filt_pes_2[i])
        filter_overrep.filter_overrep_pe(foreign_filt_pes_1[i], foreign_filt_pes_2[i], pe_fqc1_paths[i], pe_fqc2_paths[i], out_dir + out_dir_inter[4])



if __name__ == "__main__":
    if len(sys.argv) == 3:
        read_process_se(se_fq = sys.argv[1], threads = sys.argv[2])
    elif len(sys.argv) == 4:
        read_process_pe(pe_fq1 = sys.argv[1] , pe_fq2 = sys.argv[2] , threads = sys.argv[3])
    else:
        print("Usage:")
        print("For single-end reads: python3 read_process.py se_fastq_list threads [optional: out_directory, remove_inter=True]")
        print("For paired-end reads: python3 read_process.py pe_fastq_1_list pe_fastq_2_list threads [optional: out_directory, remove_inter=True]")

