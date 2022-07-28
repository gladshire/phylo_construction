import sys
import os

TRANSRATE_LOC = "~/miles/packages/transrate-1.0.3/bin/"
TRANSRATE_CMD = TRANSRATE_LOC + "transrate"

def transrate_ref(assembly, pe_fq1, pe_fq2, threads, out_dir, reference):
    if out_dir == ".": out_dir = os.getcwd()
    if out_dir != os.getcwd(): os.chdir(out_dir)
    if out_dir[-1] != "/": out_dir += "/"
    if out_dir[0] != "/": out_dir = "/" + str(out_dir)

    path_assembly, file_assembly = os.path.split(assembly)
    assembly_base_name = os.path.splitext(file_assembly)[0]
    results_name = assembly_base_name + "_transrate_results"

    cmd = [TRANSRATE_CMD, "--assembly", assembly, "--left", pe_fq1,
           "--right", pe_fq2, "--threads", str(threads), "--reference",
           reference, "--output", results_name]

    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)

def transrate_no_ref(assembly, pe_fq1, pe_fq2, threads, out_dir):
    if out_dir == ".": out_dir = os.getcwd()
    if out_dir != os.getcwd(): os.chdir(out_dir)
    if out_dir[-1] != "/": out_dir += "/"
    if out_dir[0] != "/": out_dir = "/" + str(out_dir)

    path_assembly, file_assembly = os.path.split(assembly)
    assembly_base_name = os.path.splitext(file_assembly)[0]
    results_name = assembly_base_name + "_transrate_results"

    cmd = [TRANSRATE_CMD, "--assembly", assembly, "--left", pe_fq1,
           "--right", pe_fq2, "--threads", str(threads),
           "--output", results_name]

    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)


if __name__ == "__main__":
    if len(sys.argv) == 6:
        transrate_no_ref(assembly = sys.argv[1], pe_fq1 = sys.argv[2],
                         pe_fq2 = sys.argv[3], threads = sys.argv[4],
                         out_dir = sys.argv[5])
    elif len(sys.argv) == 7:
        transrate_ref(assembly = sys.argv[1], pe_fq1 = sys.argv[2],
                      pe_fq2 = sys.argv[3], threads = sys.argv[4],
                      out_dir = sys.argv[5], reference = sys.argv[6])
    else:
        print("Usage:")
        print("python3 transrate_wrapper.py assembly_fasta fastq_file1 fastq_file2 threads output_directory [optional: reference]")
