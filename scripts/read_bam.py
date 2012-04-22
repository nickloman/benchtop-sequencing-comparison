#!/usr/bin/env python

import pysam
import sys
from runutils import read_run_details
from Bio import SeqIO

reference = dict([(rec.id, rec) for rec in SeqIO.parse(sys.argv[2], "fasta")])

def has_masked(s):
        return len([c for c in s if c.islower()])

MINIMUM_MAPPING_QUALITY = 1

print "sample\tref\trid\tmapped\tmapq\tinsertions\tl_insertions\tdeletions\tl_deletions\trlen"

samples = read_run_details(sys.argv[1])
for sample in samples:
	mapped = 0
	unmapped = 0

	samfile = pysam.Samfile(sample['Path'], "rb")
	id = 1
	for read in samfile:
		if read.is_unmapped or \
		read.mapq < MINIMUM_MAPPING_QUALITY or \
		has_masked(str(reference[samfile.getrname(read.tid)][read.pos : read.pos + read.alen].seq)):
			unmapped += 1

			print "%s\t%s\t%s\t0\t0\t0\t0\t0\t0\t%s" % (sample['Description'], sample['Reference'], id, read.qlen)
		else:
			mapped += 1
	
			discrete_insertions = 0
			discrete_deletions = 0

			length_insertions = 0
			length_deletions = 0

			for flag, bases in read.cigar:
				if flag == 1:
					discrete_insertions += 1
					length_insertions += bases
				if flag == 2:
					discrete_deletions += 1
					length_deletions += bases

			print "%s\t%s\t%s\t1\t%s\t%s\t%s\t%s\t%s\t%s" % (sample['Description'], sample['Reference'], id, read.mapq, discrete_insertions, length_insertions, discrete_deletions, length_deletions, read.qlen)
		id += 1
	print >>sys.stderr, sample['Path'], mapped, unmapped

