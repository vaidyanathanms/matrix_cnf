#!/bin/bash

#BSUB -P BIP189
#BSUB -W 00:10
#BSUB -nnodes 1
#BSUB -J pp_WT_T_600
#BSUB -o outdir/out.%J
#BSUB -e outdir/err.%J

module load caascade/1.1.beta .gcc/6.4.0 spectrum-mpi/10.3.1.2-20200121 gcc/6.4.0 spectrum-mpi/10.3.1.2-20200121
module load gromacs/2020.2-rdtscp_off

export GMX_MAXBACKUP=-1;
export OMP_NUM_THREADS=2;

echo "begin job.."
echo $PWD

mkdir -p initdir

# editconf box
jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi editconf -f WT_switchgrass_nch_20.pdb -bt cubic -d 1.0 -o initconf.gro
wait

# solvate with solvent - organic solvent
# no solvation

# solvate with cosolvent - water
# no cosolvation

# make enermin_tpr file
jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi grompp -f minim.mdp -p WT_switchgrass_nch_20.top -c initconf.gro -o enermin.tpr
wait

# generate temperature_coupling files
jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi select -s initconf.gro -sf tcgrp_inp.txt -on tcgrp_indx.ndx
wait

cp *.pdb initdir/
cp *.psf initdir/
cp *.top initdir/
cp *.txt initdir/
cp *.ndx initdir/
cp *.gro initdir/
