import sys, os
from itertools import zip_longest


def filter_summary(tot_reads, corr_reads, unfx_reads):
    print("\nSUMMARY:\n")
    print("  Total reads processed: {}".format(tot_reads))
    print("  Reads retained: {}".format(corr_reads))
    print("  Reads removed: {}\n".format(unfx_reads))


def filter_unfix_se(se_fq, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_se, file_se = os.path.split(se_fq)
    file_name = str(file_se)

    filt_name = file_name.split(".")[0] + ".fix.fq"
    
    if os.path.exists(out_dir + filt_name):
        print("Filtered file found for: " + se_fq)
        return
 
    se_in = open(se_fq, 'r')
    se_out = open(out_dir + filt_name, 'w')
    
    se_corr_ct = 0
    se_unfx_ct = 0
    tot_rds_ct = 0
    line_ct = 0
    
    with se_in as se:
        line = se.readline()
        while line:
            if line_ct % 4 == 0:
                tot_rds_ct += 1
                print("{} reads processed".format(tot_rds_ct), end = "\r")
                if "unfixable" in line:
                    se_unfx_ct += 1
                    for i in range(4):
                        line = se.readline()
                else:
                    se_corr_ct += 1
                    for i in range(4):
                        se_out.write("{}".format(line))
                        line = se.readline()
            line_ct += 1
    filter_summary(tot_rds_ct, se_corr_ct, se_unfx_ct)

def filter_unfix_pe(pe_fq1, pe_fq2, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"
    
    path_pe_1, file_pe_1 = os.path.split(pe_fq1)
    file_name_1 = str(file_pe_1)
    filt_name_1 = file_name_1.split(".")[0] + ".fix.fq"

    path_pe_2, file_pe_2 = os.path.split(pe_fq2)
    file_name_2 = str(file_pe_2)
    filt_name_2 = file_name_2.split(".")[0] + ".fix.fq"

    if os.path.exists(out_dir + filt_name_1) and \
       os.path.exists(out_dir + filt_name_2):
        print("Filtered files found for: ")
        print(pe_fq1)
        print(pe_fq2)
        return

    pe_in_1 = open(pe_fq1, 'r')
    pe_in_2 = open(pe_fq2, 'r')
    pe_out_1 = open(out_dir + filt_name_1, 'w')
    pe_out_2 = open(out_dir + filt_name_2, 'w')

    rd_corr_ct_1 = 0
    rd_corr_ct_2 = 0
    pair_corr_ct = 0
    pe_unfx_ct = 0
    tot_rds_ct = 0
    line_ct = 0

    print("Initiating removal of unfixable reads...")
    with pe_in_1 as p1, pe_in_2 as p2:
        line_1 = p1.readline()
        line_2 = p2.readline()
        while line_1 and line_2:
            if line_ct % 4 == 0:
                tot_rds_ct += 1
                print("Reads processed: {}".format(tot_rds_ct), end = "\r")
                if "unfixable" in line_1 or "unfixable" in line_2:
                    pe_unfx_ct += 1
                    for i in range(4):
                        line_1 = p1.readline()
                        line_2 = p2.readline()
                else:
                    pair_corr_ct += 1
                    if "cor" in line_1:
                        rd_corr_ct_1 += 1
                    if "cor" in line_2:
                        rd_corr_ct_2 += 1
                    if "cor" in line_1 and "cor" in line_2:
                        pair_corr_ct += 1
                    for i in range(4):
                        pe_out_1.write("{}".format(line_1))
                        pe_out_2.write("{}".format(line_2))
                        line_1 = p1.readline()
                        line_2 = p2.readline()
            line_ct += 1
    filter_summary(tot_rds_ct, pair_corr_ct, pe_unfx_ct)



if __name__ == "__main__":
    if len(sys.argv) == 3:
        filter_unfix_se(se_fq = sys.argv[1], out_dir = sys.argv[2])
    elif len(sys.argv) == 4:
        filter_unfix_pe(pe_fq1 = sys.argv[1], pe_fq2 = sys.argv[2], out_dir = sys.argv[3])
    else:
        print("Usage:")
        print("For single-end reads: python3 filter_unfixable.py fastq_reads output_directory")
        print("For paired-end reads: python3 filter_unfixable.py fastq_reads1 fast1_reads2 output_directory")
        sys.exit()
