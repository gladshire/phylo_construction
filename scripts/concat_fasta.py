import os
import sys
import subprocess


def assemble_file(fasta_files, out_file, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    output_file = out_dir + out_file
    with open(output_file, 'w') as fast_out:
        for fast_in in fasta_files:
            curr_in = open(fast_in, 'r')
            curr_in_str = curr_in.read() 
            fast_out.write(curr_in_str)
            curr_in.close()


if __name__ == "__main__":
    assemble_file(sys.argv[:-2:], sys.argv[-2], sys.argv[-1])
