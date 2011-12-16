import sys
import numpy
from Bio import SeqIO

bases = []
for n in xrange(0,2000):
	bases.append(list())

for rec in SeqIO.parse(open(sys.argv[1]), sys.argv[2]):
	for n, qual in enumerate(rec.letter_annotations["phred_quality"]):
		bases[n].append(qual)

for n, scores in enumerate(bases):
	if scores:
		print "%s\t%s" % (n, numpy.mean(scores))

