import sys
from Bio import SeqIO

fh1 = open(sys.argv[2], "w")
fh2 = open(sys.argv[3], "w")

itr = SeqIO.parse(open(sys.argv[1]), "fastq")
while 1:
	rec1 = itr.next()
	rec2 = itr.next()

	SeqIO.write([rec1], fh1, "fastq")
	SeqIO.write([rec2], fh2, "fastq")
