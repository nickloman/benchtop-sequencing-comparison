#!/usr/bin/env python

from Bio import SeqIO
import re
import sys
import glob
import os.path
from runutils import read_assemblies
from collections import defaultdict

if len(sys.argv) != 3:
	print "Usage: %s reference_genbank_file assemblies" % ( sys.argv[0] )
	raise SystemExit

assemblies = read_assemblies(sys.argv[2])

refhash = {}
for rec in SeqIO.parse(open(sys.argv[1]), "genbank"):
	# remove version
	id = rec.id.split(".")[0]
	refhash[id] = rec

def go(assembly_fn, gap_fn, out_gap_fn):
	outfh = open(out_gap_fn, "w")

	asshash = {}
	for rec in SeqIO.parse(open(assembly_fn), "fasta"):
		asshash[rec.id] = rec

	stati = defaultdict(int)
	homopolymers = defaultdict(int)

	results  = []
	for ln in open(gap_fn):
		cols = ln.rstrip().split("\t")

		if cols[0] == "reference":
			position = int(cols[2]) - 1
			contig = cols[1]
		elif cols[0] == "assembly":
			position = int(cols[7]) - 1 
			contig = cols[6]
		else:
			continue

		contig_pos = int(cols[10]) - 1
		distance_to_edge = min(contig_pos, len(asshash[cols[9]]) - contig_pos)

		homopolymer = str(refhash[contig].seq[position:position + 15])

		assembly_seq = asshash[cols[9]][contig_pos:][0:15]

		if len(homopolymer) > 1:
			regex = "(%s+)" % (homopolymer[1])
			m = re.match(regex, homopolymer[1:], re.IGNORECASE)

			homopolymer_len = len(m.group(1))
		else:
			homopolymer_len = 1

		if distance_to_edge <= 10:
			# status = "contig/scaffold gap"
			status = "Assembly gap"
		elif 'N' in homopolymer or 'n' in homopolymer:
			# skip these
			continue
			status = "reference gap"
		elif int(cols[4]) > 5:
			if 'n' in assembly_seq or 'N' in assembly_seq:
				status = "Assembly gap"
				# status = "contig/scaffold gap"
			else:
				status = "Assembly gap"
				# status = "gap"
		elif homopolymer_len == 1:
			status = "Assembly gap"
		else:
			status = "Homopolymer gap"
			homopolymers[homopolymer_len] += 1

		cols.extend([m.group(1), len(m.group(1)), distance_to_edge, status, homopolymer, homopolymer_len, str(assembly_seq.seq)])
		results.append(cols)

		print >> outfh, "\t".join([str(s) for s in cols])
#		stati[status] += 1

#	return stati, homopolymers
	return results

print "assembly\tgaptype\thplen"
for a in assemblies:
	dirname = a['Name'] + '_scores'
	gapfile = glob.glob("%s/*__gaps.txt" % (dirname,))
	if len(gapfile) != 1:
		raise SystemExit
	m = re.search('(\d+)', os.path.basename(gapfile[0]))
	alignment_number = m.group(1)

	assembly_fn = "%s/alignment%s/%s.fas" % (dirname, alignment_number, a['Name'])
	gap_fn = "%s/alignment%s__gaps.txt" % (dirname, alignment_number)
	out_gap_fn = "%s/gaps.txt" % (dirname,)

#	stati, homopolymers = go(assembly_fn, gap_fn, out_gap_fn)
	results = go(assembly_fn, gap_fn, out_gap_fn)
	for cols in results:
		print "\t".join([str(x) for x in [a['Desc'], cols[14], cols[16]]])

#	for status, value in stati.iteritems():
#		print "%s\tcategory\t%s\t%s" % (a['Name'], status, value)
#	for hplen, value in homopolymers.iteritems():
#		print "%s\thomopolymer\t%s\t%s" % (a['Name'], hplen, value)

