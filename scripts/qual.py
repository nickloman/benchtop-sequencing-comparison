from Bio import SeqIO
import sys
import numpy as np
import matplotlib.mlab as mlab

format = sys.argv[2]

quals = []
for rec in SeqIO.parse(open(sys.argv[1]), format):
	quals.extend(rec.letter_annotations["phred_quality"])

import matplotlib.pyplot as plt

n, bins, patches = plt.hist(quals, bins=range(0, 40))

# add a 'best fit' line
y = mlab.normpdf(bins, np.mean(quals), np.std(quals))
print >>sys.stderr, y
#, numpy.mean(quals), numpy.sd(quals))
# l = plt.plot(bins, y, 'r--', linewidth=1)

plt.xlabel('Phred-scaled quality score')
plt.ylabel('Frequency (bases)')
plt.title(sys.argv[3])
plt.grid(True)

plt.savefig(sys.stdout, format='png')
