import os
import sys
import subprocess
import chimera_detect

# Define post-assembly function

#   1. Run chimera detection (chimera_detect.py)
#   2. Run corset, extract representative (corset_wrapper.py)
#   3. Filter corset output (filter_corset_output.py)
#   4. Run translation with transrate (transrate_wrapper.py)

def post_assembly(prot_ref, threads):
   
    out_dir = os.getcwd() + "/"
    assembly_dir = out_dir + "06-trinity_assembly/"

    # Run chimera detection / removal
    chim_dir = "07-chimera_filtered"
    if os.path.exists(out_dir + chim_dir) == False:
        os.mkdir(out_dir + chim_dir)
    
    print("\nNow initiating removal of chimeras ...\n")
    for file_asmbly in os.listdir(assembly_dir):
        if file_asmbly.split(".")[-1] == "gene_trans_map":
            continue
        chim_cmd = ["python3", "chimera_detect.py", assembly_dir + file_asmbly,  prot_ref,
                    str(threads), out_dir + chim_dir]
        subprocess.Popen(" ".join(chim_cmd), shell = True).wait()



if __name__ == "__main__":
    
    if len(sys.argv) == 3:
        post_assembly(sys.argv[1], sys.argv[2])
    else:
        print("\nUsage:")
        print("python3 post_assembly.py proteome_ref_fasta threads")

    
