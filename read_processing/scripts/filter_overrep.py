import sys
import os
import re
import time


def parse_fastqc_log(fastqc_log):
    with open(fastqc_log) as fp:
        for result in re.findall("Overrepresented sequences(.*?)END_MODULE", fp.read(), re.S):
            seqs = ([i.split("\t")[0] for i in result.split("\n")[2:-1]])
    return seqs

def seqs_match(overrep_list, read):
    flag = False
    if overrep_list != []:
        for seq in overrep_list:
            if seq in read:
                flag = True
                break
    return flag

def overrep_summary(tot_rds_ct, tot_fail_ct):
    print("SUMMARY:")
    print("  Total reads processed: {}".format(tot_rds_ct))
    print("  Reads retained: {}".format(tot_rds_ct - tot_fail_ct))
    print("  Reads removed: {}".format(tot_fail_ct))

def filter_overrep_pe(pe_fq1, pe_fq2, pe_fqc1, pe_fqc2, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_pe_1, file_pe_1 = os.path.split(pe_fq1)
    pe_fq1_name = str(file_pe_1)
    base_name_pe_1 = pe_fq1_name.split(".")[0] + ".overrep_filtered.fq"

    path_pe_2, file_pe_2 = os.path.split(pe_fq2)
    pe_fq2_name = str(file_pe_2)
    base_name_pe_2 = pe_fq2_name.split(".")[0] + ".overrep_filtered.fq"

    fqin1 = open(pe_fq1, 'r')
    fqin2 = open(pe_fq2, 'r')
    fq1out = open(out_dir + base_name_pe_1, 'w')
    fq2out = open(out_dir + base_name_pe_2, 'w')

    leftseqs = parse_fastqc_log(pe_fqc1)
    rightseqs = parse_fastqc_log(pe_fqc2)

    tot_rds_ct = 0
    tot_fail_ct = 0
    line_num = 0

    print("Initiating removal of over-represented reads...")
    with fqin1 as f1, fqin2 as f2:
        head1 = f1.readline()
        head2 = f2.readline()
        while head1 and head2:
            if line_num % 4 == 0:
                tot_rds_ct += 1
                print("{} reads processed".format(tot_rds_ct), end = "\r")
                
                if line_num != 0:
                    head1 = f1.readline()
                    head2 = f2.readline()
         
                seq1 = f1.readline()
                seq2 = f2.readline()
                 
                tmp1 = f1.readline()
                tmp2 = f2.readline()

                qual1 = f1.readline()
                qual2 = f2.readline()

                flag_left = seqs_match(leftseqs, seq1)
                flag_right = seqs_match(rightseqs, seq2)

                if True not in (flag_left, flag_right):
                    fq1out.write("{}".format("".join([head1, seq1, tmp1, qual1])))
                    fq2out.write("{}".format("".join([head2, seq2, tmp2, qual2])))
                else:
                    tot_fail_ct += 1
            line_num += 1

        fq1out.close()
        fq2out.close()
        
        overrep_summary(tot_rds_ct, tot_fail_ct)


def filter_overrep_se(se_fq, se_fqc, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_se, file_se = os.path.split(se_fq)
    se_fq_name = str(file_se)
    base_name_se = se_fq_name.split(".")[0] + ".overrep_filtered.fq"

    sein = open(se_fq, 'r')
    seout = open(out_dir + base_name_se, 'w')
    
    se_seqs = parse_fastqc_log(se_fqc)

    tot_rds_ct = 0
    tot_fail_ct = 0
    line_num = 0

    print("Initiating removal of over-represented reads")
    with sein as f1:
        head = f1.readline()
        while head:
            if line_num % 4 == 0:
                tot_rds_ct += 1
                print("{} reads processed".format(tot_rds_ct), end = "\r")
                
                if line_num != 0:
                    head = f1.readline()
                seq = f1.readline()
                tmp = f1.readline()
                qual = f1.readline()
               
                """
                print("head " + head)
                print(" seq " + seq)
                print(" tmp " + tmp)
                print("qual " + qual)
                time.sleep(1)
                """
                se_flag = seqs_match(se_seqs, seq)
                
                if se_flag != True:
                    seout.write("{}".format("".join([head, seq, tmp, qual])))
                else:
                    tot_fail_ct += 1
            line_num += 1

        seout.close()

        overrep_summary(tot_rds_ct, tot_fail_ct)



if __name__ == "__main__":
    if len(sys.argv) == 4:
        filter_overrep_se(se_fq = sys.argv[1], se_fqc = sys.argv[2], out_dir = sys.argv[3])
    elif len(sys.argv) == 6:
        filter_overrep_pe(pe_fq1 = sys.argv[1], pe_fq2 = sys.argv[2],
                          pe_fqc1 = sys.argv[3], pe_fqc2 = sys.argv[4],
                          out_dir = sys.argv[5])
    else:
        print("Usage:")
        print("For single-end reads: python3 filter_overrep.py fastq_file fastQC_file output_directory")
        print("For paired-end reads: python3 filter_overrep.py fastq_file1 fastq_file2 fastQC_file1 fastQC_file2 output_directory")
        sys.exit()


