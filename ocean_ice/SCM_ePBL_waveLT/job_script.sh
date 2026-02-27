#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name="SCM_iH_ePBL"
#SBATCH --output=SCM_hurri_o.%j
#SBATCH --error=SCM_hurri_e.%j
#SBATCH --qos=normal
#SBATCH --partition=batch
#SBATCH --clusters=c5
#SBATCH --account=your division
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your email


# Avoid job errors because of filesystem synchronization delays
sync && sleep 1
export FI_VERBS_PREFER_XRC=0
srun --ntasks=1 --cpus-per-task=1 --export=ALL ../../build/ocean_only/MOM6

