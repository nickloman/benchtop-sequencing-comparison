
import sys
import os.path

ROOT_PATH = "../"

if len(sys.argv) != 3:
	print "Usage: %s assemblies reference" % (sys.argv[0],)
	raise SystemExit

print "PATH := ${PATH}:../mauve_snapshot_2011-08-19"
print

files = []
for ln in open(sys.argv[1]):
	if ln.startswith('#'):
		continue

	cols = ln.rstrip().split("\t")
	tag, fn = cols[0:2]

	print "%s: %s %s" % (tag, ROOT_PATH+fn, sys.argv[2])
	print "\trm -f %s" % (tag,)
	print "\tln -s %s %s" % (ROOT_PATH+fn, tag)
	print

	print "%s_scores/chromosomes.txt: %s %s" % (tag, tag, sys.argv[2])
	print "\tmauvego %s %s %s" % (tag, sys.argv[2], tag)
	print

#	print "%s_scores/gaps.txt: %s_scores/chromosomes.txt" % (tag,)

	files.append("%s_scores/chromosomes.txt" % (tag,))

for f in files:
	print "all: %s" % (f,)
