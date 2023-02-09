from os.path import join as pjoin
import numpy as np
from dipy.io.image import load_nifti, save_nifti

import nibabel as nib

path_to_patients = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study"
patient = "20_01_01_E0"
sub_folder = "microstructure\\diamond"
file_type =  "20_01_01_E0_diamond_fractions.nii.gz"
test_file = pjoin(path_to_patients, patient, sub_folder, file_type)

other_sub_folder = "microstructure\\dti"
other_file_type = "20_01_01_E0_FA.nii.gz"
other_test_file = pjoin(path_to_patients, patient, other_sub_folder, other_file_type)

ref_file = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\Atlas_Maps\\FSL_HCP1065_FA_1mm.nii.gz"

"""ref_data, ref_affine, ref_image = load_nifti(ref_file, return_img=True)
test_data, test_affine, test_image = load_nifti(test_file, return_img=True)
other_test_data, other_test_affine, other_test_image = load_nifti(other_test_file, return_img=True)


print(other_test_data.shape)
print(test_data.shape)
print(ref_data.shape)"""

final_file = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\transformed\\20_01_01_E0_diamond_fractions_transformed.nii.gz"
final_image = nib.load(final_file)
test_image = nib.load(test_file)
print(test_image)
print("==============================")
print(final_image)
