import sys
from Bio import SeqIO

class Stats:
    def __init__(self, contig_lengths, n_vals):
        self.contig_lengths = contig_lengths
        self.n_vals = n_vals

def get_stats(run):
    if run['Tech'] == '454':
        format = 'sff'
    else:
        format = 'fastq'

    fh = open("reads/%s" % (run['Filename']))
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

def read_run_details(fn):
        runs = []
        colnames = []
	headers = False
        for ln in open(fn):
                ln = ln.rstrip()
                if ln.startswith('#') and headers == False:
                        colnames = ln[1:].split("\t")
			headers = True
                        continue
		if ln.startswith('#'):
			continue

                cols = ln.split("\t")

                d = {}
                for key in colnames:
                        d[key] = ''

                for n, v in enumerate(cols):
			try:
                        	d[colnames[n]] = v
			except:
				print >>sys.stderr, "Warning misformed %s" % (ln,)

                runs.append(d)
        return runs

def hashit(list, key):
	return dict([(item[key], item) for item in list])

def read_references(fn):
	return hashit(read_run_details(fn), 'RefID')

def read_run_details_hash(fn):
	return hashit(read_run_details(fn), 'RunID')

def read_assemblies(fn):
	return read_run_details(fn)
