from Bio import SeqIO
import sys
import numpy as np
import matplotlib.mlab as mlab

format = sys.argv[2]

lengths = []
for rec in SeqIO.parse(open(sys.argv[1]), format):
	lengths.append(len(rec))

import matplotlib.pyplot as plt

n, bins, patches = plt.hist(lengths, bins=range(min(lengths), max(lengths)))

# add a 'best fit' line
y = mlab.normpdf(bins, np.mean(lengths), np.std(lengths))
print >>sys.stderr, y
#, numpy.mean(lengths), numpy.sd(lengths))
l = plt.plot(bins, y, 'r--', linewidth=1)

plt.xlabel('Read length')
plt.ylabel('Frequency')
plt.title(sys.argv[3])
plt.grid(True)

plt.savefig(sys.stdout, format='png')
