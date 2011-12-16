import matplotlib.pyplot as plt
import sys

plt.plot([1,2,3,2,3,4])

plt.savefig(sys.stdout, format='png')
