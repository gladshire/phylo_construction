import sys
import os

PIDENT_CUTOFF = 30
LENGTH_CUTOFF = 100

def qcov(hsp):
    return abs(hsp[11] - hsp[10]) + 1


def separated(hsp1, hsp2):
    length1 = qcov(hsp1)
    length2 = qcov(hsp2)
    start = min(hsp1[10], hsp1[11], hsp2[10], hsp2[11])
    end = max(hsp1[10], hsp1[11], hsp2[10], hsp2[11])
    overlap = length1 + length2 - (end - start) + 1
    if overlap < min(60, 0.2 * min(length1, length2)):
        return True
    else:
        return False


def expand_range(hsp1, hsp2):
    if hsp1 == []:
        return hsp2
    if hsp2 == []:
        return hsp1
    start1, end1, start2, end2 = hsp1[10], hsp1[11], hsp2[10], hsp2[11]
    if start1 < end1 and start2 < end2:
        start, end = min(start1, start2), max(end1, end2)
    elif start1 > end1 and start2 > end2:
        start, end = max(start1, start2), min(end1, end2)
    hsp1[10], hsp1[11] = start, end
    return hsp1


def check_block(block, multigene):
    if len(block) == 1:
        return True
    pos, neg = [], []
    for hsp in block:
        if hsp[4][0] == "-":
            neg = expand_range(neg, hsp)
        else:
            pos = expand_range(pos, hsp)
    if (pos == [] and neg != []) or (neg == [] and pos != []):
        return True
    elif separated(pos, neg):
        if multigene:
             start1, end1 = min(pos[10], pos[11]), max(pos[10], pos[11])
             start2, end2 = min(neg[10], neg[11]), max(neg[10], neg[11])
             outfile1.write(pos[0] + " " + str(int(start1)) + " " + str(int(end1)) + " trans-multi\n")
             outfile1.write(neg[0] + " " + str(int(start2)) + " " + str(int(end2)) + " trans-multi\n")
        else:
            if qcov(pos) > qcov(neg):
                outhsp = pos
            else:
                outhsp = neg
            start, end = min(outhsp[10], outhsp[11]), max(outhsp[10], outhsp[11])
            outfile1.write(outhsp[0] + " " + str(int(start)) + " " + str(int(end)) + " trans-self\n")
        for i in pos:
            outfile2.write(str(i) + "\t")
        outfile2.write("\n")
        for i in neg:
            outfile2.write(str(i) + "\t")
        outfile2.write("\n")
        return False
    else:
        return True


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage:")
        print("python3 detect_chimera_from_blastx.py blastx_output output_directory")
    
    blastx_output = sys.argv[1]
    out_dir = sys.argv[2]

    if out_dir == ".": out_dir = os.getcwd()
    if os.path.isabs(out_dir) == False: out_dir = os.path.abspath(out_dir)
    if out_dir[-1] != "/": out_dir += "/"

    path_blastx, file_blastx = os.path.split(blastx_output)
    blastx_name = str(file_blastx)
    blastx_base_name = blastx_name.split(".")[0]

    infile = open(blastx_output, "r")
    outfile1 = open(out_dir + blastx_base_name + ".cut", "w")
    outfile2 = open(out_dir + blastx_base_name + ".info", "w")
    last_query = ""

    for line in infile:
        if len(line) < 3:
            continue
        hsp = linestrip().split("\t")
        for i in [5, 10, 11]:
            hsp[i] = float(hsp[1])
        if hsp[5] < PIDENT_CUTOFF or qcov(hsp) < LENGTH_CUTOFF:
            continue
        query, hit = hsp[0], hsp[2]
 
        if last_query == "":
            hit_block = [hsp]
            query_block = [hsp]
            good_seq = True
        elif query == last_query:
            query_block.append(hsp)
            if good_seq:
                if hit == last_hit:
                     hit_block.append(hsp)
                else:
                     good_seq = check_block(hit_block, False)
                     hit_block = [hsp]
        else:
            if good_seq:
                good_seq = check_block(hit_block, False)
            if good_seq:
                good_seq = check_block(query_block, True)
            query_block, hit_block = [hsp], [hsp]
            good_seq = True
        last_query, last_hit = query, hit

    if good_seq:
        good_seq = check_block(hit_block, False)
    if good_seq:
        good_seq = check_block(query_block, True)

    infile.close()
    outfile1.close()
    outfile2.close()
