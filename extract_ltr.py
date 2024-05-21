#!/usr/local/env python3

from Bio import SeqIO
import sys

lib = sys.argv[1]

ltrs = []

for te in SeqIO.parse(lib, "fasta"):
	print(te.id.split("#")[1])
	if "LTR" in te.id.split("#")[1]:
		ltrs.append(te)
		print("added !!!!")

SeqIO.write(ltrs, lib+".ltr", "fasta")
