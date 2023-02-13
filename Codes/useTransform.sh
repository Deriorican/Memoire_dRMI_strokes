#!/bin/bash
# Submission script for Lemaitre3
#SBATCH --job-name=transformation
#SBATCH --array=/*\ 0-NUMBER_OF_PATIENTS_YOU_HAVE-1 /*\
#SBATCH --time=00:2:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5000 # megabytes
#SBATCH --partition=batch,debug
#
#SBATCH --mail-user=/*\ YOUR_EMAIL_ADRESS /*\
#SBATCH --mail-type=ALL
#
#SBATCH --output=transformation.out

PATIENTS_PATH=#/*\ "YOUR_STUDY" the path to the folder containing your patients' directories /*\

ATLAS_PATH=#/*\ "YOUR_ATLAS" the path to your atlas directory /*\

ATLAS_SUB_FOLD=#/*\ "YOUR_ATLAS_SUB_FOLD" the path, from your atlas directory, to the reference file of the atlas from which you want to register. If the reference file is in the main directory of the atlas write "-" /*\

NUM_FILES=#/*\ N the number of files you want to transform (same as the length of the following list) /*\

FILES=#/*\ ("YOUR_FILE1" "YOUR_FILE2" "YOUR_FILE3" ...) all your files to transform names/*\

NUM_PATIENTS=1

PATIENTS=#/*\ ("YOUR_PATIENT1" "YOUR_PATIENT2" "YOUR_PATIENT3" ...) all your patients' name/*\

srun python3 ./useTransform.py $PATIENTS_PATH $ATLAS_PATH $ATLAS_SUB_FOLD $NUM_FILES ${FILES[@]} $NUM_PATIENTS ${PATIENTS[$SLURM_ARRAY_TASK_ID]}
echo "Task ID: $SLURM_ARRAY_TASK_ID"