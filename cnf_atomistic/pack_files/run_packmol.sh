#!/bin/bash

#SBATCH -A bsd
#SBATCH -p burst
#SBATCH -t 0-03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH -J packmol_gen
#SBATCH -o outdir/out.%J
#SBATCH -e outdir/err.%J

module load vmd
module load python

echo "begin job.."
echo $PWD

~/packmol/packmol < pack_cellulose.inp
