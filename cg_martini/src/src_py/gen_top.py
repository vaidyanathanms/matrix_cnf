# To generate topology for cellulose/matrices

from auxgen_top import *

#----Read input file - filename, ncnf_fibers, acetfrac--------------
ch_per_cnf   = 18 # default
if len(sys.argv) == 6:
    ch_per_cnf = int(sys.argv[5]) # number of chains per bundle
elif len(sys.argv) != 5: 
    print('Unknown number of arguments: ', len(sys.argv),\
          str(sys.argv))
    exit()

print('Input file name: ',sys.argv[1])

#-----Process input data - Cellulose/Acetylated Cellulose-----------
fname    = str(sys.argv[1]) # Coarse-grained file
ncnf     = int(sys.argv[2]) # num of cellulose bundles
acetfrac = float(sys.argv[3]) # acetylated fraction
cell_dp  = int(sys.argv[4]) # degree of polymerization of cellulose

#-----Check input file is present-----------------------------------
if not os.path.exists(fname):
    raise RuntimeError(fname + ' not found in path!')

#-----Generate log file---------------------------------------------
fout = gen_logfile(fname,ncnf,acetfrac,cell_dp,ch_per_cnf)
residarr,resnamearr,aidarr,anamearr,rxarr,ryarr,rzarr,massarr = read_gro_file(fname)
                                                                
#-----Generate bead list--------------------------------------------
glycan_list = create_martini_beads(cell_dp,ncnf,ch_per_cnf,residarr,\
                                   aidarr,anamearr)
# CREATE BOND LIST
bond_list = create_bond_list(cell_dp,ncnf,ch_per_cnf,glycan_list)

# CREATE ANGLE LIST
# CREATE DIHEDRAL LIST
# CREATE ATOMTYPE LIST
# CREATE NONBONDED PARAMETER LIST
# CREATE BONDED PARAMETER LIST
# COMBINE AND WRITE

