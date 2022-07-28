import sys
import os
import re
import pandas
from Bio import SeqIO


def filter_vad_transcripts(transcripts, transrate_contigs_csv, out_dir):

    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_transcript, files_transcript = os.path.split(transcripts)
    transcripts_name = str(files_trnascript)
    base_name_transcripts = transcript_name.split(".")[0]

    good_transcripts_file = base_name_transcripts + ".good_transcripts.fa"
    bad_transcripts_file = base_name_transcripts + ".bad_transcripts.fa"
    good_transcripts_short_name_file = base_name_transcripts + ".good_transcripts.short_name.fa"

    contigs = pandas.read_csv(transrate_contigs_csv)
    transcripts_original = SeqIO.intex(transcripts, "fasta")

    Cord = contigs["sCord"]
    Cord_cutoff = Cord <= 0.50

    Ccov = contigs["sCcov"]
    Ccov_cutoff = Ccov <= 0.25

    Cnuc = contigs["sCnuc"]
    Cnuc_cutoff = Cnuc <= 0.25

    misassembly = contigs[Cord_cutoff]
    uncovered = contigs[Ccov_cutoff]
    nonagreement = contigs[Cnuc_cutoff]

    componentd = [misassembly, uncovered, nonagreement]
    reads_to_filter = pandas.concat(components)
    undup_reads_to_filter = pandas.DataFrame.drop_duplicates(reads_to_filter, keep = 'first')
    bad_transcripts_names = undup_reads_to_filter["contig_name"]

    bad_transcripts = []
    if len(bad_transcripts_names) == 0:
        print("Bad transcripts not found")
    for i in bad_transcripts_names:
        bad_transcripts.append(transcripts_original[i])
    count = SeqIO.write(bad_transcripts, out_dir + bad_transcripts_file, "fasta")
    print("Removed {} bad transcripts".format(count))
    undup_reads_to_filter.to_csv(out_dir + str(base_name_transcripts) + ".bad_transcripts.csv", index = False)

    good_transcripts_to_keep = pandas.concat([contigs, undup_reads_to_filter]).drop_duplicates(keep = False)
    good_transcripts_names = good_transcripts_to_keep["contig_name"]
    good_transcripts = []

    if len(good_transcripts_names) == 0:
        print("Good transcripts not found")
    for i in good_transcripts_names:
        good_transcripts.append(transcripts_original[i])
    count = SeqIO.write(good_transcripts, out_dir + good_transcripts_file, "fasta")
    print("Kept {} good transcripts".format(count))
    good_transcripts_to_keep.to_csv(out_dir + str(base_name_transcripts) + ".good_transcripts.csv", index = False)

    searchstr = r'(>\w+)(\slen.*)'
    replacestr = r'\1'
    infile = open(out_dir + good_transcripts_short_name_file, 'w')

    with open(out_dir + good_transcripts_file, "rU") as fasta_file:
        reg = re.compile(searchstr)
        for line in fasta_file:
            line = line.strip('\n')
            if line.startswith('>'):
                fixline = reg.sub(replacestr, line)
                infile.write(fixline + '\n')
            else:
                infile.write(line + '\n')
    infile.close()

    

if __name__ == "__main__":
    if len(sys.argv) == 4:
        filter_bad_transcripts(transcripts = sys.argv[1], transrate_contigs_csv = sys.argv[2], out_dir = sys.argv[3])
    else:
        print("Usage:")
        print("python3 filter_transcripts_transrate.py transcripts_fasta transrate_csv output_directory")
        sys.exit()
