#!/bin/env python

import numpy as np
from numpy import save
import sys
from Bio import SeqIO

"""
This function takes as input classified TEs at lineage level and transforms nucleotide sequences in
a 2D representation using one hot coding (methods 1 to 5)
"""
def conversion2d(file, maxlength):

	sequences = [x for x in SeqIO.parse(file, "fasta")]
	numSeqs = len(sequences)
	langu = ['A', 'C', 'G', 'T', 'N']
	rep2d = np.zeros((numSeqs, 5, maxlength))
	posSeq = 0
	labels = np.zeros((numSeqs, 1))

	print("Doing conversion into 2D of "+file+" using the centering method")

	for te in sequences:
		seq = str(te.seq)
		posNucl = 0

		order = -1
		if str(te.id).upper().find("_NEGATIVE") != -1:
			order = 0
		elif str(te.id).upper().find("_POSITIVE") != -1:
			order = 1

		if order != -1:	
			if len(seq) < maxlength:

				# to complete TEs with NNs centering the sequence
				times = int((maxlength-len(seq))/2)
				seq = str('N'*times+str(seq)+'N'*(times+1))[0:maxlength]
			else:
				seq = seq[0:maxlength]

			for nucl in seq:
				posLang = langu.index(nucl.upper())
				rep2d[posSeq][posLang][posNucl] = 1
				posNucl += 1

			labels[posSeq][0] = order
		else:
			print("---------- error: --------------")
			print(te.id)
			print("--------------------------------")

		if posSeq % 500 == 0:
			print("Doing "+str(posSeq) +" of "+str(numSeqs))
		posSeq += 1

	print("saving features file...")
	save(file+'_center.npy', rep2d.astype(bool))
	print("done!!")
	print("saving labels file...")
	save(file+'_center_labels.npy', labels.astype(bool))
	print("done!!")



"""
This function deletes all characters that are no DNA (A, C, G, T, N)
"""
def filter(file):
	newFile = open(file+".filtered", "w")
	for te in SeqIO.parse(file, "fasta"):
		seq = str(te.seq)
		filterDna = [x for x in seq if x.upper() in ['A', 'C', 'G', 'T', 'N']]
		newSeq = "".join(filterDna)
		newFile.write(">"+str(te.id)+"\n"+newSeq+"\n")


if __name__ == '__main__':
	seqfile = sys.argv[1]
	filter(seqfile)
	maxLen=23200
	conversion2d(seqfile+".filtered", maxLen)

