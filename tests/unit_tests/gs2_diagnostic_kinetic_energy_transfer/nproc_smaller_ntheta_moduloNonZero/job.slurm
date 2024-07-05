#!/usr/bin/env bash
#SBATCH --job-name=tprim4p0              # Job name
#SBATCH --mail-type=NONE                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tms535@york.ac.uk    # Where to send mail
#SBATCH --ntasks=4                      # Number of MPI tasks to request
#SBATCH --cpus-per-task=1                # Number of CPU cores per MPI task
#SBATCH --mem=8G                         # Total memory to request
#SBATCH --time=0-00:30:00                # Time limit (DD-HH:MM:SS)
#SBATCH --account=PET-GSPT-2019         # Project account
#SBATCH --output=job_%j.log              # Standard output and error log

#setting up GK environment
export GK_SYSTEM=viking
export MAKEFLAGS=-IMakefiles
export GK_PROJECT_CODE=PET-GSPT-2019
export OMP_NUM_THREADS=1
export PATH=${PATH}:${HOME}/bin # location of gs2 executable

# Load modules
module purge
module load gompi/2022b OpenMPI/4.1.4-GCC-12.2.0 # Requirements
module load netCDF-Fortran/4.6.0-gompi-2022b # For netcdf
module load FFTW/3.3.10-GCC-12.2.0 # For FFTW
module load OpenBLAS/0.3.21-GCC-12.2.0 # For Lapack/BLAS
module load Python/3.10.8-GCCcore-12.2.0 # For testing
module load GCC/7.3.0-2.30
module load FFTW/3.3.8-gompi-2018b
module load netCDF-Fortran/4.6.0-gompi-2022b
module load git/2.38.1-GCCcore-12.2.0-nodocs
# Often required when testing things on the front end
ulimit -s unlimited

#running GK simulation
mpirun -np 4 /users/tms535/Applications/gs2-implement/bin/gs2 test_kinetic_energy_transfer.in
