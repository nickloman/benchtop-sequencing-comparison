from runutils import read_assemblies
from stats import get_stats
from Bio import SeqIO
import sys

# usage: <ref file> <assemblies>

ref_length = sum([len(rec) for rec in SeqIO.parse(open(sys.argv[1]), "fasta")])

assemblies = read_assemblies(sys.argv[2])
stats = [(a, get_stats(a['Desc'], "../../" + a['Path'], "fasta", ref_length)) for a in assemblies]

for a, s in stats:
	print "%s\t%s\t%s" % (a['Desc'], a['AssemblySoftware'], s.n_vals[0.5])
