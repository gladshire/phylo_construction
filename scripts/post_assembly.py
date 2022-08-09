import os
import sys
import subprocess

# Define post-assembly function

#   1. Convert fastq / fq to fasta (fq_to_fa.py)
#   2. Run chimera detection (chimera_detect.py)
#   3. Run corset, extract representative (corset_wrapper.py)
#   4. Filter corset output (filter_corset_output.py)
#   5. Run translation with transrate (transrate_wrapper.py)

def post_assembly():
    
    for filename in os.listdir("05-filter_over_represented"):
        fq_to_fa_cmd = ["python3", "fq_to_fa.py", filename, 
