#!/usr/bin/env python

import sys
import numpy
from Bio import SeqIO

print "base\tqual"
for rec in SeqIO.parse(open(sys.argv[1]), sys.argv[2]):
	for n, qual in enumerate(rec.letter_annotations["phred_quality"]):
		print "%s\t%s" % (n+1, qual)
