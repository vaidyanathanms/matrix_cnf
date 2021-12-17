#!/bin/bash

#SBATCH -A bsd
#SBATCH -p burst
#SBATCH -t 0-23:00:00
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

mkdir -p outdir
mkdir -p trajfiles

#---------------------------------------------------------Generate enermin files----------------------------------
finit_inp=./enermin.tpr
if ! test -f "$finit_inp"; then
	echo "begin generating enermin.tpr.."
	# generate enermin files
	srun gmx_mpi grompp -f minim.mdp -c py_finconf -p py_topol -o enermin.tpr
	wait

fi
wait

#-----------------------------------Minimize and generate NVT_high files-------------------------------------------
fnvthigh_inp=./nvt_high.tpr 
if ! test -f "$fnvthigh_inp"; then

	echo "begin running enermin.tpr.."
	# run enermin.tpr
	srun gmx_mpi mdrun -s enermin.tpr -cpo state_min.cpt -cpi state_min.cpt -cpt 2 -g md_min.log -o traj_min.trr -e ener_min.edr -c confout_min.gro -maxh 23
	wait

	# generate nvt high files
	echo "begin generating nvt_high.tpr.."

	srun gmx_mpi grompp -f nvt_high.mdp -c confout_min.gro -p py_topol -o nvt_high.tpr 
	wait

        cp md_min.log trajfiles/md_min.log
	cp traj_min.trr trajfiles/traj_min.trr
	cp ener_min.edr trajfiles/ener_min.edr
	cp confout_min.gro trajfiles/confout_min.gro
fi
wait


#-----------------------------------Run NVT_high and generate target NVT files--------------------------------------
fnvt_inp=./nvt.tpr
if ! test -f "$fnvt_inp"; then


        echo "begin running high temperature NVT: nvt_high.tpr.."
        # run nvt_high.tpr
        srun gmx_mpi mdrun -s nvt_high.tpr -cpo state_nvt_high.cpt -cpi state_nvt_high.cpt -cpt 5 -g md_nvt_high.log -o traj_nvt_high.trr -e ener_nvt_high.edr -c confout_nvt_high.gro -pin off -maxh 23
        wait


        echo "begin generating target temperature NVT: nvt.tpr.."
        # generate target nvt files

        srun gmx_mpi grompp -f nvt.mdp -c confout_nvt_high.gro -p py_topol -o nvt.tpr
        wait

        cp md_nvt_high.log trajfiles/md_nvt_high.log
        cp traj_nvt_high.trr trajfiles/traj_nvt_high.trr
        cp ener_nvt_high.edr trajfiles/ener_nvt_high.edr
        cp confout_nvt_high.gro trajfiles/confout_nvt_high.gro

fi
wait

#------------------------------------Run target NVT and generate NPT Berendsen---------------------------------------
fnpt_berend_inp=./npt_berendsen.tpr
if ! test -f "$fnpt_berend_inp"; then

	echo "begin running nvt.tpr.."
	# run target nvt.tpr

	srun gmx_mpi mdrun -s nvt.tpr -cpo state_nvt.cpt -cpi state_nvt.cpt -cpt 5 -g md_nvt.log -o traj_nvt.trr -e ener_nvt.edr -c confout_nvt.gro -pin off -maxh 23
	wait

	echo "begin generating npt_berendsen.tpr.."
	# generate npt_berendsen files
	srun gmx_mpi grompp -f npt_berendsen.mdp -c confout_nvt.gro -p py_topol -o npt_berendsen.tpr
	wait

	cp md_nvt.log trajfiles/md_nvt.log
	cp traj_nvt.trr trajfiles/traj_nvt.trr
	cp ener_nvt.edr trajfiles/ener_nvt.edr
	cp confout_nvt.gro trajfiles/confout_nvt.gro
fi
wait

#-----------------------------------Run NPT Berendsen and generate NPT_main files------------------------------------
fnpt_inp=./npt_main.tpr
if ! test -f "$fnpt_inp"; then

	echo "begin running npt_berendsen.tpr.."
	# run npt_berendsen.tpr
	srun gmx_mpi mdrun -s npt_berendsen.tpr -cpo state_npt_berendsen.cpt -cpi state_npt_berendsen.cpt -cpt 5 -g md_npt_berendsen.log -o traj_npt_berendsen.trr -e ener_npt_berendsen.edr -c confout_npt_berendsen.gro  -pin off -maxh 23
	wait

	echo "begin generating npt_main.tpr.."
	# generate npt_main files

	srun gmx_mpi grompp -f npt_main.mdp -c confout_npt_berendsen.gro -p py_topol -o npt_main.tpr
	wait

	cp md_npt_berendsen.log trajfiles/md_npt_berendsen.log
	cp traj_npt_berendsen.trr trajfiles/traj_npt_berendsen.trr
	cp ener_npt_berendsen.edr trajfiles/ener_npt_berendsen.edr
	cp confout_npt_berendsen.gro trajfiles/confout_npt_berendsen.gro
else

        echo "begin running npt_main.tpr.."
        # run npt_main.tpr

        srun gmx_mpi mdrun -s npt_main.tpr -cpo state_npt_main.cpt -cpi state_npt_main.cpt -cpt 5 -g md_npt_main.log -o traj_npt_main.trr -e ener_npt_main.edr -c confout_npt_main.gro -pin off -maxh 23

        wait

        cp md_npt_main.log trajfiles/md_npt_main.log
        cp traj_npt_main.trr trajfiles/traj_npt_main.trr
        cp ener_npt_main.edr trajfiles/ener_npt_main.edr
        cp confout_npt_main.gro trajfiles/confout_npt_main.gro
fi
wait
#------------------------------------------------------------------------------------------------------------------------------------

echo "End of run.."
