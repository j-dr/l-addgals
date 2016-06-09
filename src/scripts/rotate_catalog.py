import numpy as np
import os
import sys
import pickle
from rot_mock_tools import *

#rotation matrix from Molly
rfile = sys.argv[3]

with open(rfile, 'r') as fp:
    rmat = pickle.load(fp)

#make sure to remove old file
try:
    os.remove(sys.argv[2])
except OSError:
    pass

#do the work
rot_mock_file(sys.argv[1],rmat,sys.argv[2])
