import numpy as np
import os
import sys
from rot_mock_tools import *

#rotation matrix from Molly
rmat = np.array([[0.382192160012073,0.054545737487712,0.922471633898423],
                 [-0.924078667040543,0.025570715732622,0.381345978892521],
                 [-0.002787462265158,-0.998183801220418,0.060177479469294]])

#make sure to remove old file
try:
    os.remove(sys.argv[2])
except OSError:
    pass

#do the work
rot_mock_file(sys.argv[1],rmat,sys.argv[2])
