import pysam
import sys

mapped = 0
unmapped = 0
dup = 0

ids = set()

samfile = pysam.Samfile(sys.argv[1], "rb")
for read in samfile:
	if read.is_unmapped:
		unmapped += 1
	else:
		mapped += 1

	if read.qname in ids:
		dup += 1
	else:
		ids.add(read.qname)

print "%d %d %d %d" % (mapped, unmapped, (100 / float(mapped + unmapped) * float(mapped)), dup)
