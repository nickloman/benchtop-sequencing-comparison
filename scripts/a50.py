from Bio import SeqIO
import sys

class Manifest:
	def __init__(self, cols):
		self.id = cols[0]
		self.path = cols[1]
		self.description = cols[2]
		self.extra = cols[3:]

def read_manifest(fn):
	samples = []
	for ln in open(fn):
		if ln.startswith('#'):
			continue
		cols = ln.rstrip().split("\t")
		samples.append(Manifest(cols))
	return samples

def stats(s, fh):
	contig_lengths = [len(rec) for rec in SeqIO.parse(fh, "fasta") if len(rec) >= MIN_LENGTH]
	contig_lengths.sort(reverse=True)
	print "\n".join([s.id + "\t" + str(l) + "\t" + s.description for l in contig_lengths])

MIN_LENGTH = 0
if len(sys.argv) == 3:
	MIN_LENGTH = int(sys.argv[2])

print "assembly\tlength\tdescription"

samples = read_manifest(sys.argv[1])
for s in samples:
	stats(s, open(s.path))

