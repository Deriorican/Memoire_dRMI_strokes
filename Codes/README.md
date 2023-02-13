# Registration User's manual ?

## What is needed ?

* registration_fixedPatient.py
* firstMapping.py
* useTransform.py
* firstMapping.sh
* useTransform.sh
* an atlas (or another patient if you want to register a patient onto another)

All other files are either out-dated versions of registrating (you can use them but nothing will be explained here) or codes that are used for my own analysis.

## What is mandatory for the code to work ?

* The organisation of your files :
Even though directories' name can be different they should have the same structure. 
That is : A directory containing another directory for each patient. 
* The naming of your files :
All files relative to your patients should also be named as PATIENTNAME_FILENAME.nii.gz

Don't hesitate to look in the github in the "TestStudy" folder to have an example.

## How to use the code ?

### Without slurm :

Run first the following command in your terminal. Of course you will have to adapt the arguments to your needs. 
* arg 1 is the path to the folder containing your patients' directories.
* arg 2 is the path to your atlas directory
* arg 3 is the path, from your atlas directory, to the reference file of the atlas from which you want to register. If the reference file is in the main directory of the atlas write "-"
* arg 4 the name of the reference file of the atlas to use for registration. WITHOUT .nii.gz !
* arg 5 is the path, from your patient directory, to the reference file of the patient from which you want to register. If the reference file is in the main directory of the patient write "-"
* arg 6 the name of the reference file to use. WITHOUT the patient name and WITHOUT .nii.gz !
* arg 7 the number of patients on which you want to use the code
* arg 8,9,... the patients on which you want to use the code
The following code is the one to use in my case. Check the git to better understand the use of the command.

py firstMapping.py "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\TestStudy" "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\Atlas_Maps" "-" "FSL_HCP1065_FA_1mm" "microstructure\\dti" "_FA" 4 "20_01_01_E0" "20_01_01_E1" "20_01_01_E2" "20_11_02_E1"

This command takes 20 minutes per patient to run in my case. 
At the end, in your patient's directory there should be :

* PATIENT_affine.pickle file
* PATIENT_header.pickle file
* PATIENT_mapping.pickle file
* PATIENT_inverse_mapping.pickle file
* transformed folder containing one .nii.gz file and a QC folder with the quality check of the registration


After running this command, you can proceed to the registration of your ROI with the following command (Once again to be adapted) :
* arg 1 is the path to the folder containing your patients' directories.
* arg 2 is the path to your atlas directory
* arg 3 is the path, from your atlas directory, to the reference file of the atlas from which you want to register. If the reference file is in the main directory of the atlas write "-"
* arg 4 the number of files from your atlas you want to register
* arg 5,6,... the name of all files from your atlas you want to register. WITHOUT .nii.gz !
* arg n the number of patients on which you want to use the code
* arg n+1,n+2,... the patients on which you want to use the code
The following code is the one to use in my case. Check the git to better understand the use of the command.

py useTransform.py "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\TestStudy" "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\Atlas_Maps" "Cerebellar" 2 "Cortico_Ponto_Cerebellum_Left" "Cortico_Ponto_Cerebellum_Right" 4 "20_01_01_E0" "20_01_01_E1" "20_01_01_E2" "20_11_02_E1"

This command should be very fast. It will add all registered files in the transformed folder of the patients.
Note that the atlas maps are registered onto the patient's space.


### With slurm

Please look at the .sh files and replace /*\ this is a commentary /*\ commentaries by your own values.

Run firstMapping.sh and wait it to finish before running useTransform.sh.

You might have to reformat the files before running them :
* Use "dos2unix YOUR_FILE.sh" to reformat the YOUR_FILE.sh files
* Use "sbatch YOUR_FILE.sh" to run YOUR_FILE.sh files.


