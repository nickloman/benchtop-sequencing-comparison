from runutils import read_run_details
import pysam
import numpy
import sys

print "Description\tReference\tExtra\tBases\t0-9\t10-19\t20-29\t30-39\t40-49\t50-59"
runs = read_run_details(sys.argv[1])
for run in runs:
	# idea from BioPython
	SCORE_OFFSET = int(run['PhredOffset'])
	q_mapping = dict() 
	for letter in range(0, 255): 
		q_mapping[chr(letter)] = letter-SCORE_OFFSET 

	samfile = pysam.Samfile(run['Path'], "rb")

	quals = []
	for read in samfile:
		quals.extend([q_mapping[q] for q in read.qual])

	output = []
	output.append(len(quals))
	output.extend(numpy.histogram(quals, bins=[0, 10, 20, 30, 40, 50, 60])[0])

	print "%s\t%s\t%s\t" % (run['Description'], run['Reference'], run['Extra']),
	print "\t".join([str(x) for x in output])


