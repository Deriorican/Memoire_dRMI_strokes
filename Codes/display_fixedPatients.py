import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from os.path import join as pjoin
import numpy as np
from dipy.viz import regtools
from dipy.data import fetch_stanford_hardi
from dipy.data.fetcher import fetch_syn_data
from dipy.io.image import load_nifti, save_nifti
import os
import pickle
from registration_fixedPatient import *
from gaussian_fit import *

def displayHistogram(study_path, patient_list, metric, mask, plot_folder_path=None, interp=True, supplementary_name = ""):
    N_patients = len(patient_list)
    alpha = 1/N_patients
    if plot_folder_path is None:
        plot_folder_path = pjoin(study_path, "Plots", mask)
    save_path = pjoin(plot_folder_path, metric + "_" + mask + supplementary_name + ".svg")
    
    fig, ax = plt.subplots()

    for patient in patient_list:
        histogram_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + mask + "_histo.csv")
        bins_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + mask + "_bins.csv")

        histo = np.loadtxt(histogram_path)
        bins = np.loadtxt(bins_path)
        
        area = np.sum(histo) * bins[-1]
        if interp:
            gauss = fitNGauss(bins, histo, 25)
            xs = np.linspace(bins[0], bins[-1], len(bins)*3)
            gaussianApprox = evaluateNGauss(xs, gauss)
            ax.plot(xs, gaussianApprox, label=patient + f"[{area}]")
            ax.fill_between(xs, np.zeros(len(bins)*3), gaussianApprox, alpha=alpha)
            #ax.plot(bins[:-1], savgol_filter(histo, len(bins) - 1, 10), label=patient) 
            #ax.fill_between(bins[:-1], np.zeros(len(bins) - 1), savgol_filter(histo, len(bins) - 1, 10), alpha=alpha) 
        else:
            ax.bar(bins[:-1], histo, width=bins[1]-bins[0], align="edge", alpha=alpha, label=patient + f"[{area}]")
    
    ax.legend()
    ax.set_title(metric + " - " + mask)
    ax.grid(linestyle = '--', linewidth = 0.5)

    ax.set_xlabel("Pixel values")
    ax.set_ylabel("Frequency")

    fig.savefig(save_path)

    print(f"{metric} - {mask} plot done.")
    #plt.show()


def main(args):
    """print(args)
    study_path = args[1]
    metric = args[2]
    mask = args[3]
    number_of_patients = int(args[4])
    patient_list = []
    for i in range(5, 5 + number_of_patients):
        patient_list.append(args[i])    """

    masks = ["xtract_prob_Corticospinal_Tract_L",
             "xtract_prob_Corticospinal_Tract_R",
             "xtract_prob_Frontal_Aslant_Tract_L",
             "xtract_prob_Frontal_Aslant_Tract_R",
             "xtract_prob_Superior_Longitudinal_Fasciculus_1_L",
             "xtract_prob_Superior_Longitudinal_Fasciculus_2_L",
             "xtract_prob_Superior_Longitudinal_Fasciculus_1_R",
             "xtract_prob_Superior_Longitudinal_Fasciculus_2_R",
             "xtract_prob_Superior_Longitudinal_Fasciculus_3_L",
             "xtract_prob_Superior_Longitudinal_Fasciculus_3_R",
             "xtract_prob_Superior_Thalamic_Radiation_L",
             "xtract_prob_Superior_Thalamic_Radiation_R",
             "Cortico_Ponto_Cerebellum_Left",
             "Cortico_Ponto_Cerebellum_Right",
             "mni_prob_Frontal_Lobe",
             "mni_prob_Parietal_Lobe",  
             "harvardoxford-subcortical_prob_Right_Thalamus",
             "harvardoxford-subcortical_prob_Left_Thalamus"]

    metrics = ['FA', "MD", "AD", "RD"]

    patient_list = ["20_01_01_E0", "20_01_01_E1", "20_01_01_E2", "20_11_02_E1"]

    study_path = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study"

    for metric in metrics:
        for mask in masks:
            displayHistogram(study_path, patient_list, metric, mask)


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