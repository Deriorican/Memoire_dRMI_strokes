import matplotlib.pyplot as plt
from os.path import join as pjoin
from os.path import exists
from os import makedirs
import numpy as np
from dipy.viz import regtools
from dipy.data import fetch_stanford_hardi
from dipy.data.fetcher import fetch_syn_data
from dipy.io.image import load_nifti, save_nifti
import os

patientPath = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E1\\microstructure\\dti\\transformed\\20_01_01_E1_transformed_FA.nii.gz"
maskPath = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\Atlas_Maps\\XTRACT\\xtract_prob_Corticospinal_Tract_L.nii.gz"

def getHisto(study_path, patient, image_type, mask_path, mask_type, nBins):
    complete_mask_path = pjoin(mask_path, mask_type + ".nii.gz")
    mask, mask_affine = load_nifti(complete_mask_path, return_img=False)
    image, image_affine = load_nifti(pjoin(study_path, patient, "transformed", patient + "_" + image_type + "_transformed" + ".nii.gz"), return_img=False)

    common_mask, common_affine = load_nifti(pjoin(study_path, "common_mask.nii.gz"))
    mask = mask * common_mask
    mask -= np.min(mask)
    mask = mask / np.max(mask)

    image -= np.min(image)
    image = image / np.max(image)
    maxBins = np.max(image[mask > 0])

    bins = np.linspace(0, maxBins, nBins+1)
    histoPoints = np.zeros(nBins)
    for i in range(nBins-1):
        indices1 = bins[i]< image 
        indices2 = image <= bins[i+1]
        indices = indices1 * indices2
        histoPoints[i] = np.sum(mask[indices])
    indices1 = bins[nBins-1] < image
    histoPoints[nBins-1] = np.sum(mask[indices1])

    histoPoints = histoPoints / np.sum(histoPoints)
    bins = np.array(bins)
    histoPoints = np.array(histoPoints)
    save_file = pjoin(study_path, patient, "histograms", patient + "_" + image_type + "_" + mask_type + "_histo" + ".csv")
    bins_file = pjoin(study_path, patient, "histograms", patient + "_" + image_type + "_" + mask_type + "_bins" + ".csv")
    if not exists(pjoin(study_path, patient, "histograms")):
        makedirs(pjoin(study_path, patient, "histograms"))
    with open(save_file, 'w') as my_file:
        np.savetxt(my_file, histoPoints)
    with open(bins_file, 'w') as my_file:
        np.savetxt(my_file, bins)
    print('Histogram exported.')
    return bins, histoPoints


def getHistoMultiplePatientsImagesMasks(study_path, patient_list, image_list, mask_path, mask_type_list, nBins):
    Histos = []
    for patient in patient_list:
        patient_histos = []
        for image in image_list:
            image_histos = []
            for mask_type in mask_type_list:
                bins, newHisto = getHisto(study_path, patient, image, mask_path, mask_type, nBins)
                image_histos.append(newHisto)
            patient_histos.append(image_histos)
        Histos.append(patient_histos)
        print(f"Histograms for patient {patient} finished.")
    return bins, Histos



def compareHisto(bins, histos, labels, ncols, titles):
    n_histos = len(histos)
    fig, axs = plt.subplots(1, ncols)
    for i in range(n_histos):
        if ncols != 1:
            a = i // (n_histos // ncols)
            axs[a].bar(bins[i][:-1], histos[i], width=bins[i][1]-bins[i][0], align="edge", alpha = 1 / n_histos * ncols, label=labels[i])
            axs[a].legend()
        else:
            axs.bar(bins[i][:-1], histos[i], width=bins[i][1]-bins[i][0], align="edge", alpha = 1 / n_histos * ncols, label=labels[i])
            axs.legend()
    for i in range(len(titles)):
        if ncols != 1:
            axs[i].set_title(titles[i])
        else:
            axs.set_title(titles[0])
    plt.show()


def getCommonMask(study_path, patient_list):
    for i in range(len(patient_list)):
        total_path = pjoin(study_path, patient_list[i], "transformed", patient_list[i] + "_FA_transformed.nii.gz")
        if i == 0:
            common_mask, common_affine = load_nifti(total_path, return_img=False)
            common_mask = common_mask > 0
        else:
            other_mask, other_affine = load_nifti(total_path, return_img=False)
            other_mask = other_mask > 0
            common_mask = common_mask * other_mask
    to_save = np.zeros(common_mask.shape)
    to_save[common_mask] = 1
    save_nifti(pjoin(study_path, "common_mask.nii.gz"), to_save, np.eye(4))
    return to_save