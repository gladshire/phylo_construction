# Uses kraken2 to label reads


import sys, os, re
from os.path import expanduser

home = expanduser("~")

if __name == "__main__":
	if len(sys.argv) != 4:
		print("Usage:")
		print("For seq extraction: python3 extract_seqs.py Order_name genome[cp,mt,both] 
