#!/bin/bash

#SBATCH -A bsd
#SBATCH -p burst
#SBATCH -t 0-00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH -J py_jobname
#SBATCH -o outdir/out.%J
#SBATCH -e outdir/err.%J

module load PE-gnu/3.0
module load gromacs/2020.6
module load vmd

export OMP_NUM_THREADS=32;

echo "begin job.."
echo $PWD

mkdir -p initdir

# editconf box
srun gmx_mpi editconf -f py_meltconf -bt cubic -d py_dval -o py_boxmeltconf
wait

# make enermin_tpr file
srun gmx_mpi grompp -f minim.mdp -p py_topol -c py_finconf -o enermin.tpr
wait

cp *.pdb initdir/
cp *.psf initdir/
cp *.top initdir/
cp *.txt initdir/
cp *.ndx initdir/
cp *.gro initdir/
