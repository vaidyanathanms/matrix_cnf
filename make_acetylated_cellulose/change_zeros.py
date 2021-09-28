# To convert zeros in input pdb files

import os
import sys
import numpy as np
import re
import shutil
import glob
import math
import struct
import subprocess
from aux_change_zeros import *

## Input data

if len(sys.argv) != 2:
    raise RuntimeError("Requires exactly one input file name")
else:
    inp_fname = str(sys.argv[1])

## Main analysis

if not os.path.exists(inp_fname):
    raise RuntimeError(inp_fname + " does not exist")

out_fname = "processed_" + inp_fname
read_and_process_file(inp_fname,out_fname)


