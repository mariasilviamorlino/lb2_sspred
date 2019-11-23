'''Usage: python3 check_profiles.py <pssm> <logfile>'''

import sys
import numpy as np

filename = sys.argv[1]
outpath = sys.argv[2]

pssm = np.loadtxt(filename)
with open(outpath,'a') as out:
    if pssm.sum() != 0:
        out.write("{0}\n".format(filename))
