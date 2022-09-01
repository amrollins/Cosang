#!/bin/bash

#SBATCH -N 24 #check nodes for this number

#SBATCH -n 384 #check nodes for this number n/N determines the amount of available memory

#SBATCH -o /work2/09034/tg883333/stampede2/GCosmo/gadget.log

#SBATCH -e /work2/09034/tg883333/stampede2/GCosmo/gadget.err

#SBATCH -p normal #normal for KNL

#SBATCH -A TG-PHY200076 #allocation ID

#SBATCH --mail-user=amrollins@crimson.ua.edu

#SBATCH --mail-type=ALL #END if you just want the last notification

#SBATCH -t 47:00:00 # set maximum run time

echo “Loading”

module load TACC intel impi hdf5 gsl fftw2

export I_MPI_COMPATIBILITY=4

cd /work2/09034/tg883333/stampede2/CoSANG_22

echo “Start the simulation”

ibrun /work2/09034/tg883333/stampede2/CoSANG_22/P_Gadget3 /work2/09034/tg883333/stampede2/CoSANG_22/gadget.param /work2/09034/tg883333/stampede2/CoSANG_22/sage.par
