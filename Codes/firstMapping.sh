#!/bin/bash
# Submission script for Lemaitre3
#SBATCH --job-name=registration
#SBATCH --array=/*\ 0-NUMBER_OF_PATIENTS_YOU_HAVE-1 /*\
#SBATCH --time=00:45:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5000 # megabytes
#SBATCH --partition=batch,debug
#
#SBATCH --mail-user=/*\ YOUR_EMAIL_ADRESS /*\
#SBATCH --mail-type=ALL
#
#SBATCH --output=registration.out

PATIENTS_PATH=#/*\ "YOUR_STUDY" the path to the folder containing your patients' directories /*\

ATLAS_PATH=#/*\ "YOUR_ATLAS" the path to your atlas directory /*\

ATLAS_SUB_FOLD=#/*\ "YOUR_ATLAS_SUB_FOLD" the path, from your atlas directory, to the reference file of the atlas from which you want to register. If the reference file is in the main directory of the atlas write "-" /*\

ATLAS_REF=#/*\ "YOUR_ATLAS_REF" the name of the reference file of the atlas to use for registration. WITHOUT .nii.gz ! /*\ 

SUB_FOLD=#/*\ "YOUR_SUB_FOLD" the path, from your patient directory, to the reference file of the patient from which you want to register. If the reference file is in the main directory of the patient write "-" /*\ 

STATIC_REF_FILE=#/*\ "YOUR_STATIC_REF_FILE" the name of the reference file to use. WITHOUT the patient name and WITHOUT .nii.gz !/*\

NUM_PATIENTS=1 # DO NOT CHANGE

PATIENTS=#/*\ ("YOUR_PATIENT1" "YOUR_PATIENT2" "YOUR_PATIENT3" ...) all your patients' name/*\

srun python3 ./firstMapping.py $PATIENTS_PATH $ATLAS_PATH $ATLAS_SUB_FOLD $ATLAS_REF $SUB_FOLD $STATIC_REF_FILE $NUM_PATIENTS ${PATIENTS[$SLURM_ARRAY_TASK_ID]}
echo "Task ID: $SLURM_ARRAY_TASK_ID"