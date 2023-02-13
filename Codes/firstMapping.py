from registration_fixedPatient import getMapping_MultiplePatients, transformMultipleFiles_MultiplePatients, qualityCheck_MultiplePatients
from os.path import join as pjoin
import sys

def main(args):
    print(args)
    study_path = args[1]
    atlas_path = args[2]
    atlas_sub_folder = args[3]
    if atlas_sub_folder == "-":
        atlas_sub_folder = ""
    file_to_transform = args[4]
    sub_folder = args[5]
    if sub_folder == "-":
        sub_folder = ""
    patient_ref = args[6]
    number_of_patients = int(args[7])
    patient_list = []
    for i in range(8 , 8 + number_of_patients):
        patient_list.append(args[i])
        
    suffixe = ".nii.gz"
    atlas_ref_path = pjoin(atlas_path, atlas_sub_folder, file_to_transform + suffixe)

    path_to_file = pjoin(atlas_path, atlas_sub_folder)
    

    mapping = getMapping_MultiplePatients(atlas_ref_path, study_path, patient_list, sub_folder, patient_ref)
    transformMultipleFiles_MultiplePatients([path_to_file], [file_to_transform], study_path, patient_list)
    qualityCheck_MultiplePatients(study_path, patient_list, sub_folder, patient_ref, file_to_transform + suffixe)

    return 1


if __name__ == "__main__":
    # Example of command : 
    # py firstMapping.py "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\TestStudy" "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\Atlas_Maps" "-" "FSL_HCP1065_FA_1mm" "microstructure\\dti" "_FA" 4 "20_01_01_E0" "20_01_01_E1" "20_01_01_E2" "20_11_02_E1"
    main(sys.argv)