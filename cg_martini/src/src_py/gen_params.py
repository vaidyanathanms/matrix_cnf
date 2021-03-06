# To generate itp files for certain polymers

import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
from auxgen_params import *
#------------------BEGIN INPUTS-----------------------

# Input data - Polymer matrix
matrix   = 'petg' #pla/petg/p3hb
nmons    = 20 # number of matrix monomers in output

if matrix == 'petg':
    design_petg(nmons,matrix.upper())
elif matrix == 'pla':
    design_pla(nmons,matrix.upper())
elif matrix == 'p3hb':
    design_p3hb(nmons,matrix.upper())
else:
    raise RuntimeError('Unknown matrix input: ' + matrix)

