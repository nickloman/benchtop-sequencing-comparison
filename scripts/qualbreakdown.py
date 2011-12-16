from Bio import SeqIO

import sys
import numpy as np
import matplotlib.mlab as mlab

# -*- coding: utf-8 -*-
class Histogram(object):
    """
    Ascii histogram
    """
    def __init__(self, data, bins=10):
        """
        Class constructor
        
        :Parameters:
            - `data`: array like object
        """
        self.data = data
        self.bins = bins
        self.h = np.histogram(self.data, bins=self.bins)
    def horizontal(self, height=4, character ='|'):
        """Returns a multiline string containing a
        a horizontal histogram representation of self.data
        :Parameters:
            - `height`: Height of the histogram in characters
            - `character`: Character to use
        >>> d = normal(size=1000)
        >>> h = Histogram(d,bins=25)
        >>> print h.horizontal(5,'|')
        106            |||
                      |||||
                      |||||||
                    ||||||||||
                   |||||||||||||
        -3.42                         3.09
        """
        his = """"""
        bars = self.h[0]/max(self.h[0])*height
        for l in reversed(range(1,height+1)):
            line = ""
            if l == height:
                line = '%s '%max(self.h[0]) #histogram top count
            else:
                line = ' '*(len(str(max(self.h[0])))+1) #add leading spaces
            for c in bars:
                if c >= np.ceil(l):
                    line += character
                else:
                    line += ' '
            line +='\n'
            his += line
        his += '%.2f'%self.h[1][0] + ' '*(self.bins) +'%.2f'%self.h[1][-1] + '\n'
        return his
    def vertical(self,height=20, character ='|'):
        """
        Returns a Multi-line string containing a
        a vertical histogram representation of self.data
        :Parameters:
            - `height`: Height of the histogram in characters
            - `character`: Character to use
        >>> d = normal(size=1000)
        >>> Histogram(d,bins=10)
        >>> print h.vertical(15,'*')
                              236
        -3.42:
        -2.78:
        -2.14: ***
        -1.51: *********
        -0.87: *************
        -0.23: ***************
        0.41 : ***********
        1.04 : ********
        1.68 : *
        2.32 :
        """
        his = """"""
        xl = ['%.2f'%n for n in self.h[1]]
        lxl = [len(l) for l in xl]
        bars = self.h[0]/max(self.h[0])*height
        his += ' '*(max(bars)+2+max(lxl))+'%s\n'%max(self.h[0])
        for i,c in enumerate(bars):
            line = xl[i] +' '*(max(lxl)-lxl[i])+': '+ character*c+'\n'
            his += line
        return his
            

format = sys.argv[2]

quals = []
longest =0 
bases  = 0
seqs   = 0
meanquals = []

bases50 = []
for rec in SeqIO.parse(open(sys.argv[1]), format):
	if len(rec) > longest:
		longest = len(rec)
	seqs += 1
	quals.extend(rec.letter_annotations["phred_quality"])

	bases50.extend(rec.letter_annotations["phred_quality"][0:50])

	meanquals.append(np.mean(rec.letter_annotations["phred_quality"]))

h = Histogram(meanquals,bins=10)
print h.horizontal(20)
print h.vertical(20)

q30_bases = float(len([q for q in quals if q >= 30]))
q20_bases = float(len([q for q in quals if q >= 20]))
q17_bases = float(len([q for q in quals if q >= 17]))
bases = float(len(quals))

print "Bases = %d " %  (bases,)
print "First 50 bases mean (%s)" % (np.mean(bases50))
print "All bases mean (%s)" % (np.mean(quals))
print "Q17+  = %d (%s)" %  (q17_bases, q17_bases / bases)
print "Q20+  = %d (%s)" %  (q20_bases, q20_bases / bases)
print "Q30+  = %d (%s)" %  (q30_bases, q30_bases / bases)
print "Longest=%d" % (longest)
print "Seqs  = %d" % (seqs,)

