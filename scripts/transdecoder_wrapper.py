import os
import sys
import shutil
import seq
import gzip
import subprocess
from retrieve_sra_data import get_tax_id


MIN_FASTA = 1000
CDS_REFERENCE = "~/miles/allele_phasing/phylo_construction/databases/"
DB_PATH = "~/miles/allele_phasing/phylo_construction/scripts/07-filter_chimeras/proteome_ref"
TRANSDEC_LOC = "~/miles/packages/TransDecoder-TransDecoder-v5.3.0/"

def fasta_ok(fasta, min_count = MIN_FASTA):
    if not os.path.exists(fasta):
        return False
    fasta_count = 0
    for i in seq.read_fasta_file(fasta):
        if len(i.seq) > 0:
            fasta_count += 1
    print("{} contains {} non-empty sequences".format(fasta, fasta_count))
    if fasta_count >= min_count:
        return True
    else:
        return False


def blastpout_ok(blastpout, min_count = MIN_FASTA):
    if not os.path.exists(blastpout):
        return False
    with open(blastpout) as infile:
        count = len(set([line.split("\t")[0] for line in infile]))
    print("{} contains {} unique query ids".format(blastpout, count))
    if count >= min_count:
        return True
    else:
        return False


def shorten_fasta_names(inname, outname, taxon_id):
    infile = open(inname, 'r')
    outfile = open(outname, 'w')
    for line in infile:
        if line[0] == ">":
            newid = line.split(" ")[0].split(".")[0].replace(">TRINITY_", "")
            outfile.write(">" + taxon_id + "@" + newid + "\n")
        else:
            outfile.write(line)
    infile.close()
    outfile.close()


def run_transdecoder_blastp(transcripts, threads, strand, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"
    
    if os.path.isabs(transcripts) == False: transcripts = os.path.abspath(transcripts)

    path_transcript, file_transcript = os.path.split(transcripts)
    base_name_transcripts = file_transcript.split(".")

    blastpout = base_name_transcripts[0] + ".blastp.outfmt6"
    
    outpep = file_transcript + ".transdecoder.pep"
    outcds = file_transcript + ".transdecoder.cds"
   
    if fasta_ok(out_dir + outpep) and fasta_ok(out_dir + outcds):
        print("Skipping transdecoder")
    else:

        allpep = file_transcript + ".transdecoder_dir/longest_orfs.pep"
        if os.path.exists(out_dir + allpep):
            print("Skipping search for long orfs")
        else:
            if strand == "stranded":
                stranded = "-S"
            else:
                stranded = ""
            cmd = ["TransDecoder.LongOrfs", "-t", transcripts, stranded]
            print(" ".join(cmd))
            os.chdir(out_dir)
            subprocess.Popen(" ".join(cmd), shell = True).wait()
            
        if blastpout_ok(out_dir + blastpout):
            print("Skipping blastp")
        else:
            cmd = ["blastp", "-query", out_dir + allpep,
                   "-db", DB_PATH, "-max_target_seqs", "1", "-outfmt", "6",
                   "-evalue", "10", "-num_threads", str(threads), ">",
                   out_dir + blastpout]
            print(" ".join(cmd))
            subprocess.Popen(" ".join(cmd), shell = True).wait()

        if fasta_ok(out_dir + outpep) and fasta_ok(out_dir + outcds):
            print("Skip finding final CDS and PEP")
        else:
            cmd = ["TransDecoder.Predict", "-t", transcripts,
                   "--retain_blastp_hits", out_dir + blastpout, "--cpu", str(threads)]
            print(" ".join(cmd))
            os.chdir(out_dir)
            subprocess.Popen(" ".join(cmd), shell = True).wait()

    shorten_fasta_names(out_dir + outpep, out_dir + taxon_id + "_" + org_name + ".pep.fa", taxon_id)
    shorten_fasta_names(out_dir + outcds, out_dir + taxon_id + "_" + org_name + ".cds.fa", taxon_id)



if __name__ == "__main__":
    if len(sys.argv) == 5:
        transcripts = sys.argv[1]
        threads = sys.argv[2]
        strand = sys.argv[3]
        out_dir = sys.argv[4]
        path_transcript, file_transcript = os.path.split(transcripts)
        base_name_transcripts = file_transcript.split(".")
        taxon_id = base_name_transcripts[0].split("_")[0]
        org_name = "_".join(base_name_transcripts[0].split("_")[:-1:])

        blastpout = base_name_transcripts[0] + ".blastp.outfmt6"
        outcds = file_transcript + ".transdecoder.cds"

        run_transdecoder_blastp(transcripts, threads, strand, out_dir)
    else:
        print("\nUsage:\n")
        print("python3 transdecoder_wrapper.py transcripts threads strand[stranded or non-stranded] output_directory")
