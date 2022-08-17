import sys
import os
import subprocess


SALMON_LOC = "~/miles/packages/Salmon-latest_linux_x86_64/bin/"
SALMON_CMD = SALMON_LOC + "salmon-0.9.1"

CORSET_LOC = "~/miles/packages/corset-1.09-linux64/"
CORSET_CMD = CORSET_LOC + "corset"

def salmon_index(transcript, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_transcript, file_transcript = os.path.split(transcript)
    salmon_base_name = file_transcript.split(".")

    index_name = salmon_base_name[0] + "_salmon_index"

    if os.path.exists(out_dir + index_name):
        print("Salmon index found for: " + salmon_base_name[0])
        return

    salm_cmd = [SALMON_CMD, "index", "-t", transcript, "-i", out_dir + index_name,
                "--type", "quasi", "-p", str(threads)]
    print(" ".join(salm_cmd))
    subprocess.Popen(" ".join(salm_cmd), shell = True).wait()


def salmon_quant_pe(transcript, index, pe_fq1, pe_fq2, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_transcript, file_transcript = os.path.split(transcript)
    salmon_base_name = file_transcript.split(".")

    index_name = salmon_base_name[0] + "_salmon_index"
    quant_name = salmon_base_name[0] + "_salmon_quant"

    x = pe_fq1.replace(",", " ")
    y = pe_fq2.replace(",", " ")

    if os.path.exists(out_dir + quant_name):
        print("Salmon quant found for: " + salmon_base_name[0])
        return
    
    salm_cmd = [SALMON_CMD, "quant", "-i", out_dir + index_name, "--dumpEq",
                "--libType", "A", "-p", str(threads), "-1", x, "-2", y, "-o",
                out_dir + quant_name]
    print(" ".join(salm_cmd))
    subprocess.Popen(" ".join(salm_cmd), shell = True).wait()


def salmon_quant_se(transcript, index, se_fq, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_transcript, file_transcript = os.path.split(transcript)
    salmon_base_name = file_transcript.split(".")
   
    index_name = salmon_base_name[0] + "_salmon_index"
    quant_name = salmon_base_name[0] + "_salmon_quant"

    x = se_fq.replace(",", " ")

    if os.path.exists(out_dir + quant_name):
        print("Salmon quant found for: " + salmon_base_name[0])
        return

    salm_cmd = [SALMON_CMD, "quant", "-i", out_dir + index_name, "--dumpEq",
                "--libType", "A", "-p", str(threads), "-r", x, "-o",
                out_dir + quant_name]
    print(" ".join(salm_cmd))
    subprocess.Popen(" ".join(salm_cmd), shell = True).wait()


def corset_salmon_eq_classes(transcript, eq_classes, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_transcript, file_transcript = os.path.split(transcript)
    salmon_base_name = file_transcript.split(".")

    clusters_name = salmon_base_name[0] + "_salmon" + "-clusters.txt"
    counts_name = salmon_base_name[0] + "_salmon" + "-counts.txt"

    if os.path.exists(out_dir + clusters_name) and os.path.exists(out_dir + counts_name):
        print("Corset-salmon files found for: " + salmon_base_name[0])
        return
    
    cors_cmd = [CORSET_CMD, "-i", "salmon_eq_classes", eq_classes, "-m", "5",
                "-p", out_dir + salmon_base_name[0] + "_salmon"]
    print(" ".join(cors_cmd))
    subprocess.Popen(" ".join(cors_cmd), shell = True).wait()


def run_pe_salmon(transcript, pe_fq1, pe_fq2, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_transcript, file_transcript = os.path.split(transcript)
    salmon_base_name = file_transcript.split(".")

    index_name = salmon_base_name[0] + "_salmon_index"
    quant_name = salmon_base_name[0] + "_salmon_quant"

    salmon_index(transcript, threads, out_dir)
    salmon_quant_pe(transcript, out_dir + index_name, pe_fq1, pe_fq2,
                       threads, out_dir)

    corset_salmon_eq_classes(transcript, out_dir + quant_name + "/aux_info/eq_classes.txt", out_dir)


def run_se_salmon(transcript, se_fq, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_transcript, file_transcript = os.path.split(transcript)
    salmon_base_name = file_transcript.split(".")

    index_name = salmon_base_name[0] + "_salmon_index"
    quant_name = salmon_base_name[0] + "_salmon_quant"

    salmon_index(transcript, threads, out_dir)
    salmon_quant_se(transcript, out_dir + index_name, se_fq, threads, out_dir)
    corset_salmon_eq_classes(transcript, out_dir + quant_name + "/aux_info/eq_classes.txt", out_dir)



if __name__ == "__main__":
    if len(sys.argv) == 5:
        run_se_salmon(transcript = sys.argv[1], se_fq = sys.argv[2],
                      threads = sys.argv[3], out_dir = sys.argv[4])
    elif len(sys.argv) == 6:
        run_pe_salmon(transcript = sys.argv[1], pe_fq1 = sys.argv[2],
                      pe_fq2 = sys.argv[3], threads = sys.argv[4],
                      out_dir = sys.argv[5])
    else:
        print("\nUsage:")
        print("For single-ended reads: python3 corset_wrapper.py transcript fastq_file threads output_directory")
        print("For paired-ended reads: python3 corset_wrapper.py transcript fastq_file1 fastq_file2 threads output_directory")
        sys.exit()
