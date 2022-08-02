import subprocess
import os
import sys


if __name__ == "__main__":
    if len(sys.argv) == 3:

        out_dir = sys.argv[2]
        in_file = os.path.abspath(sys.argv[1])

        if out_dir == ".": out_dir = os.getcwd()
        if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
        if out_dir[-1] != "/": out_dir += "/"

        path_fa, file_fa = os.path.split(sys.argv[1])
        name_fa = str(file_fa)
        fa_base_name = name_fa.split(".")[0]

        cmd_help = "'{if(NR%4==1) {printf(\">%s\\n\",substr($0,2));} else if(NR%4==2) print;}'"

        cmd = ["cat", in_file, "|", "awk", cmd_help, ">", out_dir + fa_base_name + ".fa"]

        print(" ".join(cmd))
        os.system(" ".join(cmd))
    else:
        print("Usage:")
        print("python3 fq_to_fa.py input_fastq output_directory")
        sys.exit()
