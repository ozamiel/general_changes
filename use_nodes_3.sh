#!/bin/bash
##PBS -l walltime=100:00:00
### number of process and nodes to use. each node contains up to 24 processors
##PBS -l nodes=16:ppn=24 
## 42:ppn=24+1:ppn=16

## PBS -l nodes=compute-0-30:ppn=24+compute-0-31:ppn=24+compute-0-33:ppn=24+compute-0-34:ppn=24+compute-0-35:ppn=24+compute-0-36:ppn=24+compute-0-37:ppn=24+compute-0-38:ppn=24+compute-0-39:ppn=24+compute-0-42:ppn=24+compute-0-43:ppn=24+compute-0-44:ppn=24+compute-0-45:ppn=24+compute-0-47:ppn=24+compute-0-54:ppn=24+compute-0-55:ppn=24

#PBS -l nodes=compute-0-23:ppn=24+compute-0-55:ppn=24+compute-0-51:ppn=24+compute-0-52:ppn=24+compute-0-53:ppn=24+compute-0-54:ppn=24+compute-0-32:ppn=24+compute-0-33:ppn=24+compute-0-34:ppn=24+compute-0-35:ppn=24+compute-0-36:ppn=24+compute-0-57:ppn=24+compute-0-58:ppn=24+compute-0-59:ppn=24+compute-0-61:ppn=24+compute-0-62:ppn=24+compute-0-65:ppn=24+compute-0-66:ppn=24+compute-0-67:ppn=24+compute-0-68:ppn=24+compute-0-69:ppn=24+compute-0-70:ppn=24+compute-0-71:ppn=24+compute-0-72:ppn=24
##+compute-0-3:ppn=24+compute-0-50:ppn=24+compute-0-27:ppn=24+compute-0-28:ppn=24+compute-0-22:ppn=24+compute-0-26:ppn=24+compute-0-46:ppn=24+compute-0-48:ppn=24

##+compute-0-26:ppn=24+compute-0-32:ppn=24+compute-0-33:ppn=24+compute-0-34:ppn=24+compute-0-35:ppn=24
##+compute-0-12:ppn=24+compute-0-13:ppn=24+compute-0-14:ppn=16
##+compute-0-50:ppn=24+compute-0-52:ppn=24+compute-0-56:ppn=24+compute-0-57:ppn=24+compute-0-58:ppn=24

#PBS -l pvmem=15gb
#PBS -m abe 
#PBS -M ozamiel@mail.tau.ac.il  
#PBS -N srjetm_3

hostname

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

module purge
module load mpi/openmpi-1.10.4 hdf5

mpirun -mca btl_openib_allow_ib 1 /gamma_home/ozamiel/athena2/bin/athena -r /gamma_home/ozamiel/athena2/runs_magfld/B_method/output_files_pr_cor_o/jet.00047.rst -i /gamma_home/ozamiel/athena2/runs_magfld/B_method/athinput.ejetm  -d /gamma_home/ozamiel/athena2/runs_magfld/B_method/output_files_pr_cor_o/refined_3

##mpirun -mca btl_openib_allow_ib 1 /gamma_home/ozamiel/athena/bin/athena -r /gamma_home/ozamiel/athena/runs_srjet/runs_magfld/B_method/output_files_pr_cor_n/jet.00047.rst mesh/refinement=adaptive mesh/numlevel=3 -d /gamma_home/ozamiel/athena/runs_srjet/runs_magfld/B_method/output_files_pr_cor_n/refined_3 

##mpirun -mca btl_openib_allow_ib 1 /gamma_home/ozamiel/athena2/bin/athena -i /gamma_home/ozamiel/athena2/runs_magfld/B_method/athinput.ejetm -d /gamma_home/ozamiel/athena2/runs_magfld/B_method/output_files_pr_cor_o

##mpirun -mca btl_openib_allow_ib 1 /gamma_home/pabolmasov/athena/bin/athena -i /gamma_home/pabolmasov/athena/models/athinput.ejet -d /gamma_home/pabolmasov/athena/models/ejet

## mpirun -mca btl_openib_allow_ib 1 /gamma_home/pabolmasov/athena/bin/athena -i /gamma_home/pabolmasov/athena/models/athinput.djet  -d /gamma_home/pabolmasov/athena/models/djet

## mpirun -mca btl_openib_allow_ib 1 /gamma_home/pabolmasov/athena/bin/athena -r /gamma_home/pabolmasov/athena/models/djet/jet.00049.rst -d /gamma_home/pabolmasov/athena/models/djet

