#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-03:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=main
#SBATCH --partition=batch
##SBATCH --mail-type=ALL
##SBATCH --mail-user=artem.fediai@kit.edu
#SBATCH --error out_%j_%J.err
#SBATCH --output out_%j_%J.out
#SBATCH --array=10-17

eval "$(conda shell.bash hook)"
#cd ${SLURM_SUBMIT_DIR}
conda activate my-rdkit-env

python kmean_test.py $SLURM_ARRAY_TASK_ID 10000
