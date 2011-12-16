import pysam
import sys

mapped = 0
unmapped = 0

print "rid\tmapped\tmapq\tinsertions\tl_insertions\tdeletions\tl_deletions\trlen"
samfile = pysam.Samfile(sys.argv[1], "rb")
id = 1
for read in samfile:
	if read.is_unmapped:
		unmapped += 1

		print "%s\t0\t0\t0\t0\t0\t0\t%s" % (id, read.qlen)
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

		print "%s\t1\t%s\t%s\t%s\t%s\t%s\t%s" % (id, read.mapq, discrete_insertions, length_insertions, discrete_deletions, length_deletions, read.qlen)
	id += 1

print >>sys.stderr, mapped, unmapped
