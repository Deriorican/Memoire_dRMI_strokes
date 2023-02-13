import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from os.path import join as pjoin
from os.path import exists
from os import makedirs
import numpy as np
from dipy.viz import regtools
from dipy.data import fetch_stanford_hardi
from dipy.data.fetcher import fetch_syn_data
from dipy.io.image import load_nifti, save_nifti
import os
import pickle
from registration_fixedPatient import *


def getCommonMask(study_path, patient_list, sub_folder, patient_ref, atlas_ref_path, threshold = 1e-8):
    suffixe = ".nii.gz"
    save_file = pjoin(study_path, "common_mask" + suffixe)
    atlas_ref_data, _ = load_nifti(atlas_ref_path)
    common_mask = np.ones(atlas_ref_data.shape)

    for patient in patient_list:
        inverse_mapping_path = pjoin(study_path, patient, patient + "_inverse_mapping.pickle")
        with open(inverse_mapping_path, 'rb') as handle:
            inverse_mapping = pickle.load(handle)
        patient_ref_path = pjoin(study_path, patient, sub_folder, patient + patient_ref + suffixe)
        ref_data, _ = load_nifti(patient_ref_path)
        inverted_ref_data = inverse_mapping.transform(ref_data)
        inverted_mask = inverted_ref_data < threshold
        common_mask[inverted_mask] = 0
    
    save_nifti(save_file, common_mask, np.eye(4))
    return common_mask


def applyCommonMask(study_path, patient_list):
    suffixe = ".nii.gz"
    for patient in patient_list:
        specified_save = pjoin(study_path, patient, patient + "_common_mask" + suffixe)
        transformFile(study_path, "common_mask", study_path, patient, mapping=None, specified_save=specified_save, isMask=True)
    return 1


def getHistogram(study_path, patient, sub_folder, metric, mask, nBins=100):
    suffixe = ".nii.gz"
    if not exists(pjoin(study_path, patient, "transformed")):
        makedirs(pjoin(study_path, patient, "transformed"))
    image_path = pjoin(study_path, patient, sub_folder, patient + "_" + metric + suffixe)
    mask_path = pjoin(study_path, patient, "transformed", patient + "_" + mask + suffixe)
    common_mask_path = pjoin(study_path, patient, patient + "_common_mask" + suffixe)
    
    histo_file = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + mask + "_histo" + ".csv")
    bins_file = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + mask + "_bins" + ".csv")

    common_mask_data, _ = load_nifti(common_mask_path)
    mask_data, _ = load_nifti(mask_path)
    mask_data = (mask_data - np.min(mask_data)) * common_mask_data
    mask_data = mask_data / np.max(mask_data)
    image_data, _ = load_nifti(image_path)

    maxBins = np.max(image_data[mask_data > 0])
    bins = np.linspace(0, maxBins, nBins+1)
    histoPoints = np.zeros(nBins)
    for i in range(nBins-1):
        if bins[i] == 0:
            indices1 = bins[i] < image_data
        else:
            indices1 = bins[i] <= image_data 
        indices2 = image_data < bins[i+1]
        indices = indices1 * indices2
        histoPoints[i] = np.sum(mask_data[indices])
    indices1 = bins[nBins-1] < image_data
    histoPoints[nBins-1] = np.sum(mask_data[indices1])

    histoPoints = histoPoints / (np.sum(histoPoints) * maxBins)
    bins = np.array(bins)
    histoPoints = np.array(histoPoints)

    with open(histo_file, 'w') as my_file:
        np.savetxt(my_file, histoPoints)
    with open(bins_file, 'w') as my_file:
        np.savetxt(my_file, bins)

    print('Histogram exported.')
    return histoPoints, bins


def getHistogramMultipleFiles(study_path, patient, sub_folder, metric, mask_list, nBins=100):
    all_histo = []
    all_bins = []
    for mask in mask_list:
        new_histo, new_bins = getHistogram(study_path, patient, sub_folder, metric, mask, nBins)
        all_histo.append(new_histo)
        all_bins.append(new_bins)
    return all_histo, all_bins


def getHistogramMultipleFiles_MultiplePatients(study_path, patient_list, sub_folder, metric, mask_list, nBins=100):
    all_histo = []
    all_bins = []
    for patient in patient_list:
        new_histo, new_bins = getHistogramMultipleFiles(study_path, patient, sub_folder, metric, mask_list, nBins)
        all_histo.append(new_histo)
        all_bins.append(new_bins)
    return all_histo, all_bins


def main(args):
    print(args)
    study_path = args[1]
    atlas_path = args[2]
    number_of_masks = int(args[3])
    mask_list = []
    for i in range(4, 4 + number_of_masks):
        mask_list.append(args[i])
    number_of_patients = int(args[4 + number_of_masks])
    patient_list = []
    for i in range(5 + number_of_masks, 5 + number_of_masks + number_of_patients):
        patient_list.append(args[i])

    sub_folder = pjoin("microstructure", "dti")
    patient_ref = "_FA"
    atlas_ref = "FSL_HCP1065_FA_1mm.nii.gz"
    atlas_ref_path = pjoin(atlas_path, atlas_ref)
    #getCommonMask(study_path, patient_list, sub_folder, patient_ref, atlas_ref_path)
    #applyCommonMask(study_path, patient_list)

    metric_sub_folder = pjoin("microstructure", "dti")
    metric_list = ["FA", "MD", "RD", "AD"]
    for metric in metric_list:
        getHistogramMultipleFiles_MultiplePatients(study_path, patient_list, metric_sub_folder, metric, mask_list)

    return 1


# Command : 
# py analysis_fixedPatient.py "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\study" "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\Atlas_Maps" 1 "FSL_HCP1065_FA_1mm" 4 "20_01_01_E0" "20_01_01_E1" "20_01_01_E2" "20_11_02_E1"

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

# METRICS :
# 
if __name__ == "__main__":
    main(sys.argv)