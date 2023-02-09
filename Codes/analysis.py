import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from os.path import join as pjoin
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
            #axs[a].plot(bins[i][:-1], savgol_filter(histos[i], len(bins[i]) - 1, 20)) 
            axs[a].legend()
        else:
            axs.bar(bins[i][:-1], histos[i], width=bins[i][1]-bins[i][0], align="edge", alpha = 1 / n_histos * ncols, label=labels[i])
            #axs.plot(bins[i][:-1], savgol_filter(histos[i], len(bins[i]) - 1, 20)) 
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



if False:
    study_path = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study"
    patient_list = ["20_01_01_E0", "20_01_01_E1", "20_01_01_E2", "20_11_02_E1"]
    getCommonMask(study_path, patient_list)

if False:
    study_path = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study"
    patient_list = ["20_01_01_E0", "20_01_01_E1", "20_01_01_E2", "20_11_02_E1"]
    image_path_list = ["AD", "FA", "MD", "RD"]
    mask_path = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\Atlas_Maps\\XTRACT"
    mask_type_list = ["xtract_prob_Corticospinal_Tract_L", "xtract_prob_Corticospinal_Tract_R"]

    bins, Histos = getHistoMultiplePatientsImagesMasks(study_path, patient_list, image_path_list, mask_path, mask_type_list, 100)

if True:
    histo1 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\histograms\\20_01_01_E0_FA_xtract_prob_Corticospinal_Tract_L_histo.csv")
    histo2 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E1\\histograms\\20_01_01_E1_FA_xtract_prob_Corticospinal_Tract_L_histo.csv")
    histo3 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E2\\histograms\\20_01_01_E2_FA_xtract_prob_Corticospinal_Tract_L_histo.csv")
    
    histo4 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\histograms\\20_01_01_E0_AD_xtract_prob_Corticospinal_Tract_L_histo.csv")
    histo5 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E1\\histograms\\20_01_01_E1_AD_xtract_prob_Corticospinal_Tract_L_histo.csv")
    histo6 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E2\\histograms\\20_01_01_E2_AD_xtract_prob_Corticospinal_Tract_L_histo.csv")

    histo7 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\histograms\\20_01_01_E0_RD_xtract_prob_Corticospinal_Tract_L_histo.csv")
    histo8 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E1\\histograms\\20_01_01_E1_RD_xtract_prob_Corticospinal_Tract_L_histo.csv")
    histo9 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E2\\histograms\\20_01_01_E2_RD_xtract_prob_Corticospinal_Tract_L_histo.csv")
    
    histo10 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\histograms\\20_01_01_E0_MD_xtract_prob_Corticospinal_Tract_L_histo.csv")
    histo11 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E1\\histograms\\20_01_01_E1_MD_xtract_prob_Corticospinal_Tract_L_histo.csv")
    histo12 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E2\\histograms\\20_01_01_E2_MD_xtract_prob_Corticospinal_Tract_L_histo.csv")


    histos = [histo1, histo2, histo3, histo4, histo5, histo6, histo7, histo8, histo9, histo10, histo11, histo12]

    bins1 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\histograms\\20_01_01_E0_FA_xtract_prob_Corticospinal_Tract_L_bins.csv")
    bins2 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E1\\histograms\\20_01_01_E1_FA_xtract_prob_Corticospinal_Tract_L_bins.csv")
    bins3 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E2\\histograms\\20_01_01_E2_FA_xtract_prob_Corticospinal_Tract_L_bins.csv")
    
    bins4 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\histograms\\20_01_01_E0_AD_xtract_prob_Corticospinal_Tract_L_bins.csv")
    bins5 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E1\\histograms\\20_01_01_E1_AD_xtract_prob_Corticospinal_Tract_L_bins.csv")
    bins6 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E2\\histograms\\20_01_01_E2_AD_xtract_prob_Corticospinal_Tract_L_bins.csv")

    bins7 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\histograms\\20_01_01_E0_RD_xtract_prob_Corticospinal_Tract_L_bins.csv")
    bins8 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E1\\histograms\\20_01_01_E1_RD_xtract_prob_Corticospinal_Tract_L_bins.csv")
    bins9 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E2\\histograms\\20_01_01_E2_RD_xtract_prob_Corticospinal_Tract_L_bins.csv")

    bins10 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\histograms\\20_01_01_E0_MD_xtract_prob_Corticospinal_Tract_L_bins.csv")
    bins11 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E1\\histograms\\20_01_01_E1_MD_xtract_prob_Corticospinal_Tract_L_bins.csv")
    bins12 = np.loadtxt("C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E2\\histograms\\20_01_01_E2_MD_xtract_prob_Corticospinal_Tract_L_bins.csv")

    bins = [bins1, bins2, bins3, bins4, bins5, bins6, bins7, bins8, bins9, bins10, bins11, bins12]

    labels = ["20_01_01_E0", "20_01_01_E1", "20_01_01_E2", "20_01_01_E0", "20_01_01_E1", "20_01_01_E2", "20_01_01_E0", "20_01_01_E1", "20_01_01_E2", "20_01_01_E0", "20_01_01_E1", "20_01_01_E2"]
    titles = ["FA", "AD", "RD", "MD"]
    compareHisto(bins, histos, labels, 4, titles)
