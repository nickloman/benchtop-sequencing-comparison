from Bio import SeqIO
from runutils import read_run_details
import sys

class Stats:
    def __init__(self, contig_lengths, n_vals):
        self.contig_lengths = contig_lengths
	self.n_vals = n_vals

def get_stats(seq_name, fh, format):
    contig_lengths = []
    for rec in SeqIO.parse(fh, format):
        contig_lengths.append(len(rec))

    contig_lengths.sort(reverse=True)

    n_vals = {0.5 : 0, 0.75 : 0, 0.9 : 0}

    total = sum(contig_lengths)
    for key in n_vals.keys():
        l = 0
        for n in contig_lengths:
            l += n
            if l > total * key:
                n_vals[key] = n
                break

    return Stats(contig_lengths, n_vals)

def stats(seq_name, fh, format):
    h = get_stats(seq_name, fh, format)
    print "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d" % (len(h.contig_lengths),
                                      min(h.contig_lengths),
                                      max(h.contig_lengths),
                                      sum(h.contig_lengths),
                                      sum(h.contig_lengths) / len(h.contig_lengths),
                                      h.n_vals[0.5],
                                      h.n_vals[0.75],
                                      h.n_vals[0.9])

x = read_run_details(sys.argv[1])
try:
	filter = sys.argv[2]
except:
	filter = None

print "Centre\tRunID\tStrain\tTech\tNotes\tReads\tMin\tMax\tSum\tAvg\tN50\tN75\tN90"
for r in x:
	if filter and filter != r['RunID']:
		continue

	print "%s\t%s\t%s\t%s\t%s\t" % (r['Centre'], r['RunID'], r['Strain'], r['Tech'], r['Notes']),

	if r['Tech'] == '454':
		format = 'sff'
	else:
		format = 'fastq'

	stats(r['Strain'], open("reads/%s" % (r['Filename'])), format)


