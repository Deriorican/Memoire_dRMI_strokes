from os.path import join as pjoin
from os.path import exists
from os import makedirs
import numpy as np
from dipy.viz import regtools
from dipy.data import fetch_stanford_hardi
from dipy.data.fetcher import fetch_syn_data
from dipy.io.image import load_nifti, save_nifti
from dipy.align.imaffine import (transform_centers_of_mass,
                                 AffineMap,
                                 MutualInformationMetric,
                                 AffineRegistration)
from dipy.align.transforms import (TranslationTransform3D,
                                   RigidTransform3D,
                                   AffineTransform3D)
from dipy.align import affine_registration
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.metrics import CCMetric
import os
import psutil
import pickle
import sys
import nibabel as nib


# inner psutil function
def process_memory():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss


def profile(func):
    def wrapper(*args, **kwargs):
 
        mem_before = process_memory()
        result = func(*args, **kwargs)
        mem_after = process_memory()
        print("{}:consumed memory: {:,}".format(
            func.__name__,
            mem_before, mem_after, mem_after - mem_before))
 
        return result
    return wrapper


@profile
def getMapping(atlas_ref_path, study_path, patient, sub_folder, patient_ref, nbins=32, sampling_prop=None, level_iters=[10000, 1000, 100], sigmas = [3.0, 1.0, 0.0], factors = [4, 2, 1], free_deform_level_iters = [10, 10, 5]):
    suffixe = ".nii.gz"
    patient_ref_path = pjoin(study_path, patient, sub_folder, patient + patient_ref + suffixe)
    save_direct_path = pjoin(study_path, patient, patient + "_mapping.pickle")
    save_inverse_path = pjoin(study_path, patient, patient + "_inverse_mapping.pickle")
    save_header_path = pjoin(study_path, patient, patient + "_header.pickle")
    save_affine_path = pjoin(study_path, patient, patient + "_affine.pickle")

    static_data, static_affine, static_image = load_nifti(patient_ref_path, return_img=True)
    static = static_data

    moving_data, moving_affine = load_nifti(atlas_ref_path)
    moving = moving_data

    pipeline = ["center_of_mass", "translation", "rigid", "affine"]
    xformed_img, reg_affine = affine_registration(
        moving,
        static,
        moving_affine=moving_affine,
        static_affine=static_affine,
        nbins=nbins,
        metric="MI",
        pipeline=pipeline,
        level_iters=level_iters,
        sigmas=sigmas,
        factors=factors)

    prealign = reg_affine
    free_deform_metric = CCMetric(3)
    free_deform_level_iters = [10, 10, 5]
    sdr = SymmetricDiffeomorphicRegistration(free_deform_metric, free_deform_level_iters)
    mapping = sdr.optimize(static, moving, static_affine, moving_affine, prealign=prealign)

    with open(save_direct_path, 'wb') as handle:
        pickle.dump(mapping, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(save_header_path, 'wb') as handle:
        pickle.dump(static_image.header, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(save_affine_path, 'wb') as handle:
        pickle.dump(static_affine, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(save_header_path, 'wb') as handle:
        pickle.dump(static_image.header, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Mapping computed for patient {patient}")

    static_data, static_affine = load_nifti(atlas_ref_path)
    static = static_data

    moving_data, moving_affine = load_nifti(patient_ref_path)
    moving = moving_data

    pipeline = ["center_of_mass", "translation", "rigid", "affine"]
    xformed_img, reg_affine = affine_registration(
        moving,
        static,
        moving_affine=moving_affine,
        static_affine=static_affine,
        nbins=nbins,
        metric="MI",
        pipeline=pipeline,
        level_iters=level_iters,
        sigmas=sigmas,
        factors=factors)

    prealign = reg_affine
    free_deform_metric = CCMetric(3)
    free_deform_level_iters = [10, 10, 5]
    sdr = SymmetricDiffeomorphicRegistration(free_deform_metric, free_deform_level_iters)
    mapping = sdr.optimize(static, moving, static_affine, moving_affine, prealign=prealign)

    with open(save_inverse_path, 'wb') as handle:
        pickle.dump(mapping, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Inverse mapping computed for patient {patient}")
    return mapping


def getMapping_MultiplePatients(atlas_ref_path, study_path, patient_list, sub_folder, patient_ref, nbins=32, sampling_prop=None, level_iters=[10000, 1000, 100], sigmas = [3.0, 1.0, 0.0], factors = [4, 2, 1], free_deform_level_iters = [10, 10, 5]):
    all_mappings = []
    for patient in patient_list:
        all_mappings.append(getMapping(atlas_ref_path, study_path, patient, sub_folder, patient_ref, nbins, sampling_prop, level_iters, sigmas, factors, free_deform_level_iters))
    return all_mappings


def transformFile(path_to_file, file_to_transform, study_path, patient, mapping=None, specified_save=None, isMask=False):
    suffixe = ".nii.gz"
    complete_path_to_file = pjoin(path_to_file, file_to_transform + suffixe)
    if not exists(pjoin(study_path, patient, "transformed")):
        makedirs(pjoin(study_path, patient, "transformed"))
    if specified_save is None:
        save_path = pjoin(study_path, patient, "transformed", patient + "_" + file_to_transform + suffixe)
    else:
        save_path = specified_save
    if mapping is None:
        mapping_path = pjoin(study_path, patient, patient + "_mapping.pickle")
        with open(mapping_path, 'rb') as handle:
            mapping = pickle.load(handle)
    
    header_path = pjoin(study_path, patient, patient + "_header.pickle")
    with open(header_path, 'rb') as handle:
        prev_header = pickle.load(handle)

    affine_path = pjoin(study_path, patient, patient + "_affine.pickle")
    with open(affine_path, 'rb') as handle:
        affine = pickle.load(handle)
    
    data_to_transform, affine_to_transform = load_nifti(complete_path_to_file)
    transformed = mapping.transform(data_to_transform)
    if isMask:
        transformed[transformed > 0] = 1

    new_img = nib.Nifti1Image(transformed, affine, header=prev_header)
    header = new_img.header
    header['pixdim'] = prev_header['pixdim'] # for example
    # Update the header information in the NIfTI file
    new_img.update_header()
    nib.save(new_img, save_path)
    print(f"{file_to_transform + suffixe} transformed for patient {patient}.")
    return transformed


def transformMultipleFiles(path_to_file_list, file_to_transform_list, study_path, patient, mapping=None):
    all_transformed = []
    for (path_to_file, file_to_transform) in zip(path_to_file_list, file_to_transform_list):
        all_transformed.append(transformFile(path_to_file, file_to_transform, study_path, patient, mapping))
    return all_transformed


def transformMultipleFiles_MultiplePatients(path_to_file_list, file_to_transform_list, study_path, patient_list, mapping=None):
    all_transformed = []
    for i, patient in enumerate(patient_list):
        if mapping is not None:
            mapping_patient = mapping[i]
        else :
            mapping_patient = None
        all_transformed.append(transformMultipleFiles(path_to_file_list, file_to_transform_list, study_path, patient, mapping_patient))
    return all_transformed

def qualityCheck(study_path, patient, sub_folder, patient_ref, atlas_ref):
    suffixe = ".nii.gz"
    patient_ref_path = pjoin(study_path, patient, sub_folder, patient + patient_ref + suffixe)
    atlas_ref_path = pjoin(study_path, patient, "transformed", patient + "_" + atlas_ref)
    save_root = pjoin(study_path, patient, "transformed", "QC", patient + "_qc")
    if not exists(pjoin(study_path, patient, "transformed", "QC")):
        makedirs(pjoin(study_path, patient, "transformed", "QC"))

    atlasReference, atlasReference_affine = load_nifti(atlas_ref_path)
    patientReference, patientReference_affine= load_nifti(patient_ref_path)

    regtools.overlay_slices(patientReference, atlasReference, None, 0,
                            "Global Reference", "Patient Reference", save_root + "_0.png")
    regtools.overlay_slices(patientReference, atlasReference, None, 1,
                            "Global Reference", "Patient Reference", save_root + "_1.png")
    regtools.overlay_slices(patientReference, atlasReference, None, 2,
                            "Global Reference", "Patient Reference", save_root + "_2.png")
    print(f"patient {patient} checked.")
    return 1


def qualityCheck_MultiplePatients(study_path, patient_list, sub_folder, patient_ref, atlas_ref):
    for patient in patient_list:
        qualityCheck(study_path, patient, sub_folder, patient_ref, atlas_ref)
    return 1

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

    sub_folder = pjoin("microstructure", "dti")
    patient_ref = "_FA"
    atlas_ref = "FSL_HCP1065_FA_1mm.nii.gz"
    atlas_ref_path = pjoin(atlas_path, atlas_ref)

    path_to_file = pjoin(atlas_path, atlas_sub_folder)
    path_to_file_list = [path_to_file for file_to_transform in file_to_transform_list]

    #mapping = getMapping_MultiplePatients(atlas_ref_path, study_path, patient_list, sub_folder, patient_ref)
    transformMultipleFiles_MultiplePatients(path_to_file_list, file_to_transform_list, study_path, patient_list)
    #qualityCheck_MultiplePatients(study_path, patient_list, sub_folder, patient_ref, atlas_ref)

    return 1


# Commands : 
# py registration_fixedPatient.py "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study" "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\Atlas_Maps" "-" 1 "FSL_HCP1065_FA_1mm" 4 "20_01_01_E0" "20_01_01_E1" "20_01_01_E2" "20_11_02_E1"

# XTRACT (12)
# "xtract_prob_Corticospinal_Tract_L"
# "xtract_prob_Corticospinal_Tract_R"
# "xtract_prob_Frontal_Aslant_Tract_L"
# "xtract_prob_Frontal_Aslant_Tract_R"
# "xtract_prob_Superior_Longitudinal_Fasciculus_1_L"
# "xtract_prob_Superior_Longitudinal_Fasciculus_2_L"
# "xtract_prob_Superior_Longitudinal_Fasciculus_1_R"
# "xtract_prob_Superior_Longitudinal_Fasciculus_2_R"
# "xtract_prob_Superior_Longitudinal_Fasciculus_3_L"
# "xtract_prob_Superior_Longitudinal_Fasciculus_3_R"
# "xtract_prob_Superior_Thalamic_Radiation_L"
# "xtract_prob_Superior_Thalamic_Radiation_R"

# Cerebellar (2)
# "Cortico_Ponto_Cerebellum_Left"
# "Cortico_Ponto_Cerebellum_Right"

# Lobes (2)
# "mni_prob_Frontal_Lobe"
# "mni_prob_Parietal_Lobe"

# Harvard (2)
# "harvardoxford-subcortical_prob_Right_Thalamus"
# "harvardoxford-subcortical_prob_Left_Thalamus"


if __name__ == "__main__":
    main(sys.argv)

    
    
