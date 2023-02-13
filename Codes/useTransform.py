from registration_fixedPatient import transformMultipleFiles_MultiplePatients
from os.path import join as pjoin
import sys

def main(args):
    print(args)
    study_path = args[1]
    atlas_path = args[2]
    atlas_sub_folder = args[3]
    if atlas_sub_folder == "-":
        atlas_sub_folder = ""
    number_of_file_to_transform = int(args[4])
    file_to_transform_list = []
    for i in range(5, 5 + number_of_file_to_transform):
        file_to_transform_list.append(args[i])
    number_of_patients = int(args[5 + number_of_file_to_transform])
    patient_list = []
    for i in range(6 + number_of_file_to_transform, 6 + number_of_file_to_transform + number_of_patients):
        patient_list.append(args[i])

    path_to_file = pjoin(atlas_path, atlas_sub_folder)
    path_to_file_list = [path_to_file for file_to_transform in file_to_transform_list]

    transformMultipleFiles_MultiplePatients(path_to_file_list, file_to_transform_list, study_path, patient_list)

    return 1
    

if __name__ == "__main__":
    # Example of command : 
    # py useTransform.py "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\TestStudy" "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\Atlas_Maps" "Cerebellar" 2 "Cortico_Ponto_Cerebellum_Left" "Cortico_Ponto_Cerebellum_Right" 4 "20_01_01_E0" "20_01_01_E1" "20_01_01_E2" "20_11_02_E1"
    main(sys.argv)