#!/usr/bin/python

import sys
from Bio import SeqIO
from itertools import izip

fh1 = SeqIO.parse(open(sys.argv[1], "r"), "fastq-illumina")
fh2 = SeqIO.parse(open(sys.argv[2], "r"), "fastq-illumina")

for rec1, rec2 in izip(fh1, fh2):
	SeqIO.write([rec1[:-1]], sys.stdout, "fastq-illumina")
	SeqIO.write([rec1[-1:] + rec2], sys.stdout, "fastq-illumina")
