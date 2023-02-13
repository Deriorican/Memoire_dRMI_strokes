#!/bin/bash
# Submission script for Lemaitre3
#SBATCH --job-name=registration
#SBATCH --array=0-46
#SBATCH --time=00:2:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5000 # megabytes
#SBATCH --partition=batch,debug
#
#SBATCH --mail-user=louis.lovat@student.uclouvain.be
#SBATCH --mail-type=ALL
#
#SBATCH --output=registration.out

PATIENTS_PATH="/home/users/l/l/llovat/stroke_pilab/study/subjects/"

STATIC_REF_FILE="/home/users/l/l/llovat/Memoire/Atlas_Maps/"

ATLAS_SUB_FOLD="Harvard"

NUM_FILES=2

FILES=("harvardoxford-subcortical_prob_Right_Thalamus" "harvardoxford-subcortical_prob_Left_Thalamus")

NUM_PATIENTS=1

PATIENTS=("20_11_02_E1" "20_11_02_E0" "20_11_02_E2" "20_01_01_E0" "20_02_01_E0" "20_08_02_E1" "20_10_01_E1" "20_10_02_E0" "20_06_02_E3" "20_05_01_E1" "20_06_02_E1" "20_05_01_E2" "20_01_01_E1" "20_05_01_E3" "20_03_01_E0" "20_02_02_E3" "20_08_02_E0" "20_02_02_E2" "20_05_02_E0" "20_08_02_E2" "20_06_01_E1" "20_06_02_E2" "20_07_01_E2" "20_06_01_E0" "20_10_02_E3" "20_10_02_E1" "20_08_01_E0" "20_10_01_E2" "20_05_02_E1" "20_02_01_E2" "20_08_02_E3" "20_07_01_E0" "20_03_01_E1" "20_11_02_E3" "20_06_01_E2" "20_08_01_E2" "20_10_02_E2" "20_01_01_E2" "20_02_02_E1" "20_02_02_E0" "20_02_01_E1" "20_07_01_E1" "20_05_01_E0" "20_06_02_E0" "20_05_02_E2" "20_08_01_E1" "20_03_01_E2")

srun python3 ./registration_fixedPatient.py $PATIENTS_PATH $STATIC_REF_FILE $ATLAS_SUB_FOLD $NUM_FILES ${FILES[@]} $NUM_PATIENTS ${PATIENTS[$SLURM_ARRAY_TASK_ID]}
echo "Task ID: $SLURM_ARRAY_TASK_ID"