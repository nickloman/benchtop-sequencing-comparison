#!/usr/bin/env python

# usage: ace_to_alignment gap_file orientation_file


from Bio.Sequencing import Ace
from Bio.Align.Generic import Alignment
from Bio.Alphabet import IUPAC, Gapped
import sys

#!/usr/bin/env python

from Bio import SeqIO
import re
import sys
import gc

def cut_ends(read, start, end):
	'''Replace residues on either end of a sequence with gaps.
 
	In this case we want to cut out the sections of each read which the assembler has 
	decided are not good enough to include in the contig and replace them with gaps 
	so the other information about positions in the read is maintained
	''' 
	return (start-1) * '-' + read[start-1:end] + (len(read)-end) * '-'
 
def pad_read(read, start, conlength):
	''' Pad out either end of a read so it fits into an alignment.
 
	The start argument is the position of the first base of the reads sequence in
	the contig it is part of. If the start value is negative (or 0 since ACE
	files count from 1, not 0) we need to take some sequence off the start
	otherwise each end is padded to the length of the consensus with gaps.
	'''
	if start < 1:
		seq = read[-1*start+1:]
	else:
		seq = (start-1) * '-' + read
	seq = seq + (conlength-len(seq)) * '-'
	return seq

def parse_ace(ace_file):
	ace_gen = Ace.parse(open(ace_file, 'r'))
	contig = ace_gen.next()
	align = Alignment(Gapped(IUPAC.ambiguous_dna, "-"))
	align.add_sequence(contig.name, contig.sequence)
 
	for readn in range(len(contig.reads)):
		clipst = contig.reads[readn].qa.qual_clipping_start
		clipe = contig.reads[readn].qa.qual_clipping_end
		start = contig.af[readn].padded_start
		seq = cut_ends(contig.reads[readn].rd.sequence, clipst, clipe)

		seq = pad_read(seq, start, len(contig.sequence))
		align.add_sequence(contig.reads[readn].rd.name + "_" + contig.af[readn].coru, seq)

	return contig, align

def find_padded(unpadded_coord, read):
	for n, c in enumerate(read.seq):
		if c == '*':
			continue
		unpadded_coord -= 1
		if not unpadded_coord:
			break
	return n

def unpadded_len(sequence):
	 return len(sequence) - sequence.count('*')

def get_col(align, padded_to_find):
	col = ''
	orientation = ''
	for read in align[1:]:
		letter = str(read[padded_to_find:padded_to_find+1].seq)
		if letter != '-':
			col += letter
			orientation += str(read.id.split("_")[1])
	return col, orientation

def stars_vs_letters(col):
	count1 = col.count('*')
	count2 = len(col) - count1
	return abs(count1 - count2)

def variation(c1, c2):
	return cmp( stars_vs_letters(c1[0]), stars_vs_letters(c2[0]) )

def count_columns(align, padded_to_find):
	cols = []
	for n in xrange(0,10,1):
		cols.append( get_col(align, padded_to_find+n) )
	cols.sort(cmp = variation)
	return cols

def forward_reverse_bias(column, orientation):
	star_f = 0
	star_r = 0
	letter_f = 0
	letter_r = 0

	for n, c in enumerate(column):
		o = orientation[n]
		if c == '*' and o == 'U':
			star_f += 1
		if c == '*' and o == 'C':
			star_r += 1
		if c != '*' and o == 'U':
			letter_f += 1
		if c != '*' and o == 'C':
			letter_r += 1
	return [star_f, star_r, letter_f, letter_r]

def go(acefile, col_to_find, reverse_complement, n_to_return, line_to_print):
	global my_cached_ace
	global my_cached_ace_file

	if acefile != my_cached_ace_file:
		 my_cached_ace_file = acefile
		 my_cached_ace = parse_ace(acefile)
	contig_info, align = my_cached_ace

	if reverse_complement:
		 col_to_find = unpadded_len(contig_info.sequence) - col_to_find

	# print "find %s" % (col_to_find,)

	padded_to_find = find_padded(col_to_find, align[0])
	# print padded_to_find

	#n_to_find = 7
	#for read in align:
	#	txt = str(read[padded_to_find:padded_to_find+n_to_find].seq)
	#	if txt != "-" * n_to_find:
	#		 print "%-30s %s" % (read.id, str(read[padded_to_find:padded_to_find+7].seq))

	cols = count_columns(align, padded_to_find)
	for n, col in enumerate(cols[0:n_to_return]):
		stats = forward_reverse_bias(col[0], col[1])
		
		outlist = line_to_print
	#	outlist.extend([acefile, col_to_find, n_to_return, n, len(col[0]), "\t".join([str(s) for s in stats]), col[0].count('*'), stars_vs_letters(col[0]), col[0], col[1]])
		outlist.extend([stats[0], stats[1], float(stats[0]) / float(stats[0] + stats[1]), stats[2], stats[3], float(stats[2]) / float(stats[2] + stats[3])])

		print "\t".join(str(s) for s in outlist)

def read_orientations(fn):
	orientations = {}
	for ln in open(fn):
		 cols = ln.rstrip().split()
		 if not cols or cols[0] != "contig":
			continue

		 if cols[3] == "forward":
			orientations[cols[1]] = 0
		 else:
			orientations[cols[1]] = 1
	return orientations

my_cached_ace = None
my_cached_ace_file = 'hello'

def main():
	ori = read_orientations(sys.argv[2])

	print "star_f, star_r, letter_f, letter_r, num_stars, stars_vs_letters"

	for ln in open(sys.argv[1]):
		cols = ln.rstrip().split("\t")
		if cols[0] == "reference":
			position = int(cols[10])
			contig = cols[9]
		elif cols[0] == "assembly":
			position = int(cols[2])
			contig = cols[1]
		else:
			continue

		if cols[14] == "homopolymer gap":
			go(contig + ".ace", position, ori[contig], int(cols[4]), cols)
		else:	
			print ln,
	
		gc.collect()

main()
