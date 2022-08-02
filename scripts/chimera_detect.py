import sys
import os
import subprocess
import pandas as pd
from Bio import SeqIO

SCRIPTS_HOME = os.path.expanduser("~/miles/allele_phasing/phylo_construction/scripts/")


def make_blast_db(proteome_ref, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_ref, file_ref = os.path.split(proteome_ref)
    ref_name = str(file_ref)
    ref_base_name = ref_name.split(".")[0]

    db_file_1 = ref_base_name + ".phr"
    db_file_2 = ref_base_name + ".pin"
    db_file_3 = ref_base_name + ".psq"

    cmd = ["makeblastdb", "-in", proteome_ref, "-dbtype", "prot", "-out",
           out_dir + ref_base_name]

    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell = True)


def run_blastx(transcripts, blast_db, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_transcript, file_transcript = os.path.split(transcripts)
    transcript_name = str(file_transcript)
    transcript_base_name = transcript_name.split(".")[0]

    blastx_output_file = transcript_base_name + ".blastx"

    cmd = ["blastx", "-db", blast_db, "-query", transcripts, "-evalue", "0.01", "-outfmt",
           '"6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore"', "-out", out_dir + blastx_output_file, "-num_threads",
           str(threads), "-max_target_seqs", "100"]

   
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell = True)


def run_chimera_detection(blastx_output_file, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_blastx, file_blastx = os.path.split(blastx_output_file)
    blastx_name = str(file_blastx)
    blastx_base_name = blastx_name.split(".")[0]

    cut_file = blastx_base_name + ".cut"
    info_file = blastx_base_name + ".info"

    cmd = ["python3", SCRIPTS_HOME + "detect_chimera_from_blastx.py",
           blastx_output_file, out_dir]

    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell = True)


def remove_chimeras_from_fasta(transcripts, info_file, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_transcript, file_transcript = os.path.split(transcripts)
    transcript_name = str(file_transcript)
    transcript_base_name = transcript_name.split(".")[0]

    filtered_transcripts = transcript_base_name + ".filtered_transcripts.fa"
    chimeric_transcripts = transcript_base_name + ".chimeric_transcripts.fa"

    transcripts_original = SeqIO.index(transcripts, "fasta")

    df = pd.read_table(info_file, header = None)
    blastx_columns = "qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore empty".strip().split(" ")
    df.columns = blastx_columns
    chimera_names = (df["qseqid"]).drop_duplicates()

    if len(chimera_names) == 0:
        print("No chimeras found")
    else:
        chimeras = []
        for i in chimera_names:
            chimeras.append(transcripts_original[i])
        count = SeqIO.write(chimeras, out_dir + chimeric_transcripts, "fasta")
        print("Removed {} chimeras".format(count))

        chimera_list = chimera_names.tolist()
        all_transcripts_ids = sorted(transcripts_original)
        non_chimeras_names = [x for x in all_transcripts_ids if x not in chimera_list]


    if len(chimera_names) == 0:
        print("No non-chimeras found")
    else:
        non_chimeras = []
        for i in non_chimeras_names:
            non_chimeras.append(transcripts_original[i])
        count = SeqIO.write(non_chimeras, out_dir + filtered_transcripts, "fasta")
        print("Retained {} transcripts".format(count))



if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage:")
        print("python3 chimera_detect.py transcripts_fasta ref_proteome_fasta threads output_directory")
        sys.exit()

    transcripts = os.path.abspath(sys.argv[1])
    proteome_ref = os.path.abspath(sys.argv[2])
    threads = int(sys.argv[3])
    out_dir = os.path.abspath(sys.argv[4])

    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_ref, file_ref = os.path.split(proteome_ref)
    ref_name = str(file_ref)
    ref_base_name = ref_name.split(".")[0]

    path_transcript, file_transcript = os.path.split(transcripts)
    transcript_name = str(file_transcript)
    transcript_base_name = transcript_name.split(".")[0]
    blastx_output_file = transcript_base_name + ".blastx"

    path_blastx, file_blastx = os.path.split(out_dir + blastx_output_file)
    blastx_name = str(file_blastx)
    blastx_base_name = blastx_name.split(".")[0]
    info_file = blastx_base_name + ".info"

    make_blast_db(proteome_ref, out_dir)
    run_blastx(transcripts, out_dir + ref_base_name, threads, out_dir)
    run_chimera_detection(out_dir + blastx_output_file, out_dir)
    remove_chimeras_from_fasta(transcripts, out_dir + info_file, out_dir)
