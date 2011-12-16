import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

d = SeqIO.to_dict(SeqIO.parse(open(sys.argv[1]), "fasta"))

for ln in sys.stdin:
	cols = ln.rstrip().split()
	if cols[4] == 'D':
		start = int(cols[1])
		end = int(cols[2])

		seq = d[cols[0]][start - 1 : end]
		rec = SeqRecord(seq=seq.seq, id=cols[5], description='')

		SeqIO.write([rec], sys.stdout, "fasta")


