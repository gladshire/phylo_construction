import os
import sys
import subprocess
import chimera_detect
from process_reads import get_layout
from Bio import Entrez

#   1. Run chimera detection (chimera_detect.py)
#   2. Run corset, extract representative (corset_wrapper.py)
#   3. Filter corset output (filter_corset_output.py)
#   4. Run translation with transdecoder (transdecoder_wrapper.py)

def post_assembly(prot_ref, threads, mult_samples = False):
   
    out_dir = os.getcwd() + "/"
    assembly_dir = out_dir + "06-trinity_assembly/"
    overrep_dir = out_dir + "05-filter_over_represented/"
    
    # Run chimera detection / removal
    chim_dir = "07-filter_chimera"
    if os.path.exists(out_dir + chim_dir) == False:
        os.mkdir(out_dir + chim_dir)
    
    print("\nRemoving chimera transcripts ...\n")
    for file_asmbly in os.listdir(assembly_dir):
        if file_asmbly.split(".")[-1] == "gene_trans_map":
            continue
        chim_cmd = ["python3", "chimera_detect.py", assembly_dir + file_asmbly,  prot_ref,
                    str(threads), out_dir + chim_dir]
        subprocess.Popen(" ".join(chim_cmd), shell = True).wait()
    
    # Run corset, extract representative clusters
    cors_dir = "08-clustering/"
    
    if os.path.exists(out_dir + cors_dir) == False:
        os.mkdir(out_dir + cors_dir)

    print("\nGenerating corset clusters ...\n")
    files_skip = []
    for file_overrep in os.listdir(overrep_dir):
        if file_overrep in files_skip or "comb." in file_overrep:
            continue
        curr_sra = file_overrep.split("_")[0]
        if get_layout(curr_sra) == "SINGLE":
            file_transcript = file_overrep.split(".")[0] + ".trinity.Trinity.fasta"
            
            cors_cmd = ["python3", "corset_wrapper.py", assembly_dir + file_transcript,
                        overrep_dir + file_overrep, threads, out_dir + cors_dir]
            subprocess.Popen(" ".join(cors_cmd), shell = True).wait()
        elif get_layout(curr_sra) == "PAIRED":
            if mult_samples:
                file_comps = file_overrep.split(".")
                file_transcript = "_".join(file_comps[0].split("_")[2:4:]) + "_comb.trinity.Trinity.fasta"
                org_name = "_".join(file_comps[0].split("_")[2:4])
                files_overrep_1 = [overrep_dir + file1 for file1 in os.listdir(overrep_dir) \
                                   if org_name + "_1" in file1]
                files_overrep_2 = [overrep_dir + file2 for file2 in os.listdir(overrep_dir) \
                                   if org_name + "_2" in file2]
                cors_cmd = ["python3", "corset_wrapper.py", assembly_dir + file_transcript,
                            ",".join(files_overrep_1), ",".join(files_overrep_2), threads,
                            out_dir + cors_dir]
                files_skip += files_overrep_1
                files_skip += files_overrep_2
                print(" ".join(cors_cmd))
                subprocess.Popen(" ".join(cors_cmd), shell = True).wait()
            else:
                file_comps = file_overrep.split(".")
                file_transcript = file_comps[0][:-2:] + ".trinity.Trinity.fasta"
                file_overrep_1 = file_comps[0][:-2:] + "_1" + "." + ".".join(file_comps[1::])
                file_overrep_2 = file_comps[0][:-2:] + "_2" + "." + ".".join(file_comps[1::])
                cors_cmd = ["python3", "corset_wrapper.py", assembly_dir + file_transcript,
                            overrep_dir + file_overrep_1, overrep_dir + file_overrep_2,
                            threads, out_dir + cors_dir]
                files_skip.append(file_overrep_1)
                files_skip.append(file_overrep_2)
                subprocess.Popen(" ".join(cors_cmd), shell = True).wait()
                            
    
    # Filter corset output
    print("\nFiltering corset output ...\n")
    for file_asmbly in os.listdir(assembly_dir):
        if file_asmbly.split(".")[-1] == "gene_trans_map":
            continue
        file_base = file_asmbly.split(".")[0]
        cors_filt_cmd = ["python3", "filter_corset_output.py", assembly_dir + file_asmbly,
                         out_dir + cors_dir + file_base + "_salmon-clusters.txt", out_dir + cors_dir]
        subprocess.Popen(" ".join(cors_filt_cmd), shell = True).wait()

    # Run transdecoder to translate
    trans_dir = "09-translate/"
    if os.path.exists(out_dir + trans_dir) == False:
        os.mkdir(out_dir + trans_dir)

    print("\nTranslating ...\n")
    for file_asmbly in os.listdir(out_dir + cors_dir):
        if "largest_cluster_transcripts.fa" not in file_asmbly:
            continue
        file_base = file_asmbly.split(".")[0]
        trans_cmd = ["python3", "transdecoder_wrapper.py", out_dir + cors_dir + file_asmbly,
                     str(threads), "non-stranded", out_dir + trans_dir]
        subprocess.Popen(" ".join(trans_cmd), shell = True).wait()


if __name__ == "__main__":
    
    if len(sys.argv) == 3:
        post_assembly(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4 and "mult-samples" in sys.argv:
        post_assembly(sys.argv[1], sys.argv[2], mult_samples = True)
    else:
        print("\nUsage:")
        print("python3 post_assembly.py proteome_ref_fasta threads [optional: mult-samples]")

    
