import os
import sys
import subprocess

def assemble_file(files_concat, out_dir, out_file):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    output_file = out_dir + out_file
    cmd = ["cat"] + files_concat + [">", output_file]
    subprocess.Popen(cmd, shell = True).wait()


if __name__ == "__main__":
    assemble_file(sys.argv[1:-2:], sys.argv[-2], sys.argv[-1])
