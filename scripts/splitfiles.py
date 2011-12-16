import sys
import os
import tempfile
from Bio import SeqIO

def read_sff_lines(f, size, sffaccs, counter):
	for rec in SeqIO.parse(open(f), sys.argv[3]):
		clip_left = rec.annotations["clip_qual_left"]
		clip_right = rec.annotations["clip_qual_right"]

		rec_len = clip_right - clip_left

		counter += rec_len

		print >>sffaccs, rec.id
		if counter > size:
			return 0
	return counter

def go(fn, size):
	sffaccs = tempfile.NamedTemporaryFile(delete=False)

	fh = open(fn, "w")
	files = [ln.rstrip() for ln in open(sys.argv[2])]
	counter = 0
	for f in files:
		counter = read_sff_lines(f, size, sffaccs, counter)
		if not counter:
			break

	sffaccs.close()
	print sffaccs.name

	os.system("sfffile -i %s -o %s %s" % (sffaccs.name, fn, " ".join(files)))
	os.system("rm %s" % (sffaccs.name))

MULTIPLE = float(4700000)

tag = sys.argv[1]

for n in [12.5]:
#1, 2, 5, 10, 15, 20, 25, 30, 50]:
	print n, go("%s%sx_in.iontor.%s" % (tag, n, sys.argv[3]), n * MULTIPLE)
