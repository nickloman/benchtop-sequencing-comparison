
from Bio import SeqIO
import sys
import re

hps = {}
for rec in SeqIO.parse(open(sys.argv[1]), "fasta"):
	homopolymers = re.findall("(A{2,})", str(rec.seq))
	homopolymers.extend(re.findall("(C{2,})", str(rec.seq)))
	homopolymers.extend(re.findall("(T{2,})", str(rec.seq)))
	homopolymers.extend(re.findall("(G{2,})", str(rec.seq)))

	for h in homopolymers:
		try:
			hps[len(h)] += 1		
		except:
			hps[len(h)] = 1

print hps

counts = {}

for ln in open(sys.argv[2]):
	cols = ln.rstrip().split("\t")
	if cols[2] != "Homopolymer gap":
		continue

	hplen = int(cols[3])
	try:
		counts[(cols[0], cols[1])][hplen] += 1
	except:
		counts[(cols[0], cols[1])] = {}
		for k in hps.keys():
			counts[(cols[0], cols[1])][k] = 0
		counts[(cols[0], cols[1])][hplen] += 1

for tag in counts.keys():
	for k in hps.keys():
		print "%s\t%s\t%s\t%s" % (tag[0], tag[1], k, float(counts[tag][k]) / float(hps[k]))
