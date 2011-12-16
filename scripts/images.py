from runutils import read_run_details
import os
import sys

runs = read_run_details(sys.argv[1])
for run in runs:
	print run['RunID']

	if run['Tech'] == "454":
		filetype = 'sff'
	elif run['Tech'] == "MiSeq":
		filetype = 'fastq-illumina'
	elif run['Tech'] == "Ion Torrent":
		filetype = 'fastq'
	else:
		print "unsupported format"
		raise SystemExit

	title = "%s - %s - Quality Scores" % (run['RunID'], run['Tech'])

	cmd = "python qual.py reads/%s %s \"%s\" > images/%s_qual.png" % (
		run['Filename'],
		filetype,
		title,
		run['RunID'])

	os.system(cmd)

