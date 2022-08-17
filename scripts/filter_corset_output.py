import sys
import os
import pandas as pd
from Bio import SeqIO


def filter_corset(transcripts, corset_cluster, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_transcript, file_transcript = os.path.split(transcripts)
    base_name_transcripts = file_transcript.split(".")

    largest_cluster_transcripts = base_name_transcripts[0] + ".largest_cluster_transcripts.fa"
    redundant_transcripts = base_name_transcripts[0] + ".redundant_cluster_transcripts.fa"

    if os.path.exists(out_dir + largest_cluster_transcripts) and\
       os.path.exists(out_dir + redundant_transcripts):
        print("Largest / redundant transcript files found for: " + base_name_transcripts[0])
        return

    clusters_df = pd.read_table(corset_cluster, header = None)
    cluster_columns = "seqid cluster".strip().split(' ')
    clusters_df.columns = cluster_columns

    seqid = []
    length = []

    for rec in SeqIO.parse(transcripts, "fasta"):
        name = rec.id
        seq = rec.seq
        seqLen = len(rec)
        seqid.append(name)
        length.append(seqLen)
    d = {"seqid": seqid, "length": length}
    seq_len_df = pd.DataFrame(d)

    seq_len_filtered_df = seq_len_df[seq_len_df["seqid"].isin(clusters_df["seqid"])]

    clusters_with_len_df = pd.merge(clusters_df, seq_len_filtered_df, on = "seqid", how = "left")

    largest_cluster_df = clusters_with_len_df.sort_values("length", ascending = False).drop_duplicates("cluster").sort_index()

    removed_cluster_df = clusters_with_len_df[~clusters_with_len_df["seqid"].isin(largest_cluster_df["seqid"])]

    transcripts = SeqIO.index(transcripts, "fasta")

    largest_cluster = []
    largest_cluster_name = largest_cluster_df["seqid"]

    for i in largest_cluster_name:
        largest_cluster.append(transcripts[i])
    
    count = SeqIO.write(largest_cluster, out_dir + largest_cluster_transcripts,
                        "fasta")

    print("Kept {} largest transcripts from corset clusters".format(count))

    largest_cluster_df.to_csv(out_dir + base_name_transcripts[0] + ".largest_cluster.csv", index = False)

    removed_cluster = []
    removed_cluster_name = removed_cluster_df["seqid"]

    for i in removed_cluster_name:
        removed_cluster.append(transcripts[i])

    count = SeqIO.write(removed_cluster, out_dir + redundant_transcripts,
                        "fasta")

    print("Removed {} redundant transcripts".format(count))

    removed_cluster_df.to_csv(out_dir + base_name_transcripts[0] + ".redundant_cluster.csv", index = False)



if __name__ == "__main__":
    if len(sys.argv) == 4:
        filter_corset(transcripts = sys.argv[1], corset_cluster = sys.argv[2],
                      out_dir = sys.argv[3])
    else:
        print("\nUsage:")
        print("python3 filter_corset_output.py transcripts_fasta corset_cluster_file output_directory")
        sys.exit()
