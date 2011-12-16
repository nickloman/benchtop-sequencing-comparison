from runutils import read_run_details_hash, get_stats
from operator import itemgetter
import numpy
import sys

def modal(lVals):
	#http://stackoverflow.com/questions/6252280/find-the-most-frequent-number-in-a-numpy-vector
	return max(map(lambda val: (lVals.count(val), val), set(lVals)))[1]

runs = read_run_details_hash(sys.argv[1])
print runs.keys()

to_show = {
	'HPA2'  : '454 Junior (1)',
	'HPA4'  : '454 Junior (2)',
	'BHAM5_1.5' : 'Ion Torrent (1)', 
	'BHAM6_1.5' : 'Ion Torrent (2)',
	'ILMN1_filtered' : 'MiSeq (1)',
	'ILMN1_L5_280_filtered' : 'MiSeq (strain 280)',
	'ILMN1_L6_282_filtered' : 'MiSeq (strain 282)',
	'ILMN1_L1_627_filtered' : 'MiSeq (strain 627)',
	'ILMN1_L2_283_filtered' : 'MiSeq (strain 283)',
	'ILMN1_L3_541_filtered' : 'MiSeq (strain 541)',
	'ILMN1_L4_518_filtered' : 'MiSeq (strain 518)',
	'ILMN1_L7_540_filtered' : 'MiSeq (strain 540)',
	}

def get_alignment_stats(aln):
	try:
		fn = "mapping/bwa/%s.coverage.txt.R.stats" % (aln['RunID'],)
		fh = open(fn)
		cols = fh.readline().split()
		return ["%.02f" % (float(f)) for f in cols[2:5]]
	except:
		return ['-', '-', '-']

for id, description in sorted(to_show.items(), key=itemgetter(1)):
	s = get_stats(runs[id])
	cov = get_alignment_stats(runs[id])
	fields = [ runs[id]['Display'], len(s.contig_lengths), sum(s.contig_lengths), modal(s.contig_lengths), "%s (%s)" % (int(numpy.mean(s.contig_lengths)), int(numpy.std(s.contig_lengths))), cov[0], cov[1] ]
	print " & ".join([str(s) for s in fields]) + "\\\\"

