from os.path import join as pjoin
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
def getTransform(moving_file, static_file, nbins=32, sampling_prop=None, level_iters=[10000, 1000, 100], sigmas = [3.0, 1.0, 0.0], factors = [4, 2, 1], free_deform_level_iters = [10, 10, 5], savetransform_file = None):
    static_data, static_affine, static_image = load_nifti(static_file, return_img=True)
    static = static_data

    moving_data, moving_affine, moving_image = load_nifti(moving_file, return_img=True)
    moving = moving_data

    metric = MutualInformationMetric(nbins, sampling_prop)
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

    if savetransform_file is not None:
        with open(savetransform_file, 'wb') as handle:
            pickle.dump(mapping.transform, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return mapping.transform


def transformFiles(file_list, transform, transformed_file_list, ref_shape=(182,218,182)):
    for file, transformed_file in zip(file_list, transformed_file_list):
        moving_data, moving_affine, moving_image = load_nifti(file, return_img=True)
        if len(moving_data.shape) == 3:
            transformed = transform(moving_data)
        elif len(moving_data.shape) == 5:
            final_shape = (ref_shape[0], ref_shape[1], ref_shape[2], 1, moving_data.shape[-1])
            transformed = np.zeros(final_shape)
            for i in range(moving_data.shape[-1]):
                transformed[:,:,:,0,i] = transform(moving_data[:,:,:,0,i])
                save_nifti(transformed_file + "_" + str(i) +".nii.gz", transformed[:,:,:,0,i], np.eye(4))
            print(moving_data.shape)
            print(transformed.shape)


def getTransform_MultiplePatients(path_to_patients, patient_list, patient_reference_file, file_subpath, static_reference_file):
    for patient in patient_list:
        patient_path = pjoin(path_to_patients, patient)
        complete_patient_reference_file = pjoin(file_subpath, patient + patient_reference_file + ".nii.gz")
        global_path_patient_reference_file = pjoin(patient_path, complete_patient_reference_file)
        transform_saving_file = pjoin(patient_path, patient + "_transform.pickle")
        transform = getTransform(global_path_patient_reference_file, static_reference_file, savetransform_file=transform_saving_file)
        print(f"Transform of patient {patient} calculated.")


def transformFiles_MultiplePatients(path_to_patients, patient_list, file_list, file_subpath, ref_shape=(182,218,182)):
    for patient in patient_list:
        patient_path = pjoin(path_to_patients, patient)
        transform_saving_file = pjoin(patient_path, patient + "_transform.pickle")
        with open(transform_saving_file, "rb") as input_file:
            transform = pickle.load(input_file)
        complete_file_list = [pjoin(patient_path, file_subpath, patient + file + ".nii.gz") for file in file_list]
        transformed_file_folder = pjoin(patient_path, "transformed")
        transformed_file_list = [pjoin(transformed_file_folder, patient + file + "_transformed") for file in file_list]
        transformFiles(complete_file_list, transform, transformed_file_list, ref_shape)
        print(f"All patient {patient}'s files transformed and saved.")


def qualityCheck(path_to_patients, patient_list, static_reference_file):
    for patient in patient_list:
        patient_reference_file_path = pjoin(path_to_patients, patient, "transformed", patient + "_FA_transformed.nii.gz")
        globalReference, globalReference_affine, globalReference_image = load_nifti(static_reference_file, return_img=True)
        patientReference, patientReference_affine, patientReference_image = load_nifti(patient_reference_file_path, return_img=True)

        regtools.overlay_slices(globalReference, patientReference, None, 0,
                                "Global Reference", "Patient Reference", pjoin(path_to_patients, patient, "transformed", patient + "_check_0.png"))
        regtools.overlay_slices(globalReference, patientReference, None, 1,
                                "Global Reference", "Patient Reference", pjoin(path_to_patients, patient, "transformed", patient + "_check_1.png"))
        regtools.overlay_slices(globalReference, patientReference, None, 2,
                                "Global Reference", "Patient Reference", pjoin(path_to_patients, patient, "transformed", patient + "_check_2.png"))
        print(f"patient {patient} checked.")
    print(f"All patients checked.")


def main(args):
    path_to_patients = args[1]
    static_reference_file = args[2]
    patient_list = []
    for i in range(3, len(args)):
        patient_list.append(args[i])
    file_subpath = "microstructure\\diamond"
    file_list = ["_diamond_logKappa", "_diamond_fractions"]
    transformFiles_MultiplePatients(path_to_patients, patient_list, file_list, file_subpath)



if __name__ == "__main__":
    main(sys.argv)