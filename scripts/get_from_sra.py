#!/usr/bin/python

import os
import sys

fullfn = sys.argv[1]
dir = fullfn[0:6]
fn = fullfn[0:9]

cmd = "/home/nick/.aspera/connect/bin/ascp -QT -i /home/nick/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/%s/%s/%s ." % (dir, fn, fullfn)
print cmd
os.system(cmd)


