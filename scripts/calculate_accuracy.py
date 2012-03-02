#!/usr/bin/env python

from Bio import SeqIO
import pysam
import re
import sys
import random
from array import array

MINIMUM_MAPPING_QUALITY = 1

def to_dict(tags):
	return dict(
		[(tag, val) for tag, val in tags]
	)


class Results:
	def __init__(self):
		self.results = {}

	def update_result(self, zero_index_base, base_quality, bin_to_update):
		base_quality = ord(base_quality) - 33
		base = zero_index_base + 1
		try:
			base = zero_index_base + 1
			self.results[base][base_quality][bin_to_update] += 1
		except:
			if base not in self.results:
				self.results[base] = {}
			if base_quality not in self.results[base]:
				self.results[base][base_quality] = {}
				self.results[base][base_quality]['matches'] = 0
				self.results[base][base_quality]['indel_mismatches'] = 0
				self.results[base][base_quality]['subst_mismatches'] = 0
				self.results[base][base_quality]['unaligned'] = 0

				self.results[base][base_quality][bin_to_update] += 1

	def add_match(self, zero_index_base, base_quality):
		self.update_result(zero_index_base, base_quality, 'matches')

	def add_indelmismatch(self, zero_index_base, base_quality):
		self.update_result(zero_index_base, base_quality, 'indel_mismatches')

	def add_substmismatch(self, zero_index_base, base_quality):
		self.update_result(zero_index_base, base_quality, 'subst_mismatches')

	def add_unaligned(self, zero_index_base, base_quality):
		self.update_result(zero_index_base, base_quality, 'unaligned')

def has_masked(s):
	return len([c for c in s if c.islower()])

results = Results()

reference = dict([(rec.id, rec) for rec in SeqIO.parse(sys.argv[2], "fasta")])

masked = 0
reads = 0

samfile = pysam.Samfile(sys.argv[1], "rb")
for read in samfile:
#	print read.pos
#	print read.aend, read.alen
#	print samfile.getrname(read.tid)
	# unmapped

	reads += 1

	if read.is_unmapped:
		# unmapped reads get added to unaligned counts
		for n, base in enumerate(read.seq):
			results.add_unaligned(n, read.qual[n])
		continue

	if read.mapq < MINIMUM_MAPPING_QUALITY:
		continue

	if has_masked(
		str(reference[samfile.getrname(read.tid)][read.pos : read.pos + read.alen].seq)
	):
		masked += 1
		continue

	tags = to_dict(read.tags)
	num_mismatches = tags['NM']

	align_str = []
	for flag, bases in read.cigar:
		if flag == 4: # soft clipped base
			align_str += 'S' * bases
		elif flag == 1: # insertion
			align_str += 'I' * bases
		elif flag == 2: # deletion
			pass
		elif flag == 0: # regular match
			align_str += 'M' * bases
		else:
			print "Unsupported CIGAR flag"
			raise SystemExit

	counter = 0	
	for flag, bases in read.cigar:
		if flag == 2: # deletion
			r = random.randint(0, 1) # randomly assign deletion to an adjacent base as per Green / Ewing
			try:
				if r == 1:
					if align_str[counter + 1] != 'M':
						print >>sys.stderr, "over write error"
					align_str[counter + 1] = 'D'
				else:
					if align_str[counter] != 'M':
						print >>sys.stderr, "over write error"
					align_str[counter] = 'D'
			except Exception, e:
				print >>sys.stderr, e
				print >>sys.stderr, read
		else:
			counter += bases

	if read.is_reverse:
		seq = read.seq[::-1]
		qual = read.qual[::-1]
		align_str.reverse()
		cigar = read.cigar[::-1]
	else:
		seq = read.seq
		qual = read.qual
		cigar = read.cigar

        if cigar[0][0] == 4:
                print read
		print read.is_reverse
                print "first base soft clipped"

	x = 0
	for n, base in enumerate(seq):
		flag = align_str[n]

		if flag == 'S':
			# soft masked bases get added to unaligned count
			results.add_unaligned(n, qual[n])
		else:
			if base == '=' and flag == 'M':
				results.add_match(n, qual[n])
			elif flag == 'M': # SNP
				results.add_substmismatch(n, qual[n])
				x += 1
			else:
				results.add_indelmismatch(n, qual[n])
				x += 1

#	print x, num_mismatches
#	if x != num_mismatches:
#		print read
#		print seq
#		print qual
#		print cigar
#		print align_str
#
#		print "WOOPS"
#		raise SystemExit

print >>sys.stderr, "reads: %d, masked: %d" % (reads, masked,)
print "pos\tqual\tvariable\tvalue"
for pos in sorted(results.results.keys()):
	for qual in sorted(results.results[pos].keys()):
		for key in sorted(results.results[pos][qual].keys()):
			print "%s\t%s\t%s\t%s" % (pos, qual, key, results.results[pos][qual][key])
