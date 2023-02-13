import matplotlib.pyplot as plt
import numpy as np
from os.path import join as pjoin
from os.path import exists
from os import makedirs
from gaussian_fit import evaluateNGauss, fitNGauss
from display_fixedPatients import interpolateHisto
from scipy.optimize import minimize
import sys

def evaluateSigmoid(x, sigmoid):
    A = sigmoid[0]
    mu = sigmoid[1]
    offset = sigmoid[2]
    l = sigmoid[3]
    value =   A * (2 / (1 + np.exp(-l * (x - mu))) - 1 - offset)
    return value


def objective_function(params, l, xs, ref):
    sigmoid = [params[0], params[1], params[2], l]
    sigmoidEval = evaluateSigmoid(xs, sigmoid)
    RMSE = (np.sum((ref - sigmoidEval)**2) / len(xs))**(1/2)
    return RMSE


def fitSigmoid(xs, ref, NGauss=25, r=0.99):
    X = (xs[-1] - xs[0]) / (NGauss - 1)
    l = - 2 * np.log(1 / (r / 2 + 1 / 2) - 1) / X
    A = 0
    mus = np.linspace(xs[0], xs[-1], 10)
    offset = 0
    best_res = []
    best_score = np.inf
    for mu in mus:
        params = [A, mu, offset]
        res = minimize(objective_function, params, args=(l, xs, ref))
        if res.success:
            if res.fun < best_score:
                best_score = res.fun
                best_res = res.x
    sigmoid = [best_res[0], best_res[1], best_res[2], l]
    return sigmoid


def displayContralateral(study_path, patient, metric, left_mask, right_mask, plot_folder_path=None):
    alpha = 1/3
    if plot_folder_path is None:
        plot_folder_path = pjoin(study_path, "Contralateral_Comparison_Plots", patient)
    if not exists(plot_folder_path):
        makedirs(plot_folder_path)

    save_path = pjoin(plot_folder_path, metric + "_" + left_mask + ".svg")
    
    fig, axs = plt.subplots(3, 1)

    left_bins_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + left_mask + "_bins.csv")
    left_histogram_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + left_mask + "_histo.csv")
    left_bins = np.loadtxt(left_bins_path)
    left_histo = np.loadtxt(left_histogram_path)
    right_bins_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + right_mask + "_bins.csv")
    right_histogram_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + right_mask + "_histo.csv")
    right_bins = np.loadtxt(right_bins_path)
    right_histo = np.loadtxt(right_histogram_path)

    if left_bins[-1] > right_bins[-1]:
        right_histo = interpolateHisto(right_bins, right_histo, left_bins)
        bins = left_bins
    else:
        left_histo = interpolateHisto(left_bins, left_histo, right_bins)
        bins = right_bins

    xs = np.linspace(bins[0], bins[-1], len(bins)*3)
    left_gauss = fitNGauss(bins, left_histo, 25)
    left_gaussianApprox = evaluateNGauss(xs, left_gauss)
    right_gauss = fitNGauss(bins, right_histo, 25)
    right_gaussianApprox = evaluateNGauss(xs, right_gauss)
    left_x = np.sum(left_gaussianApprox * xs) / np.sum(left_gaussianApprox)
    left_y = np.mean(left_gaussianApprox)
    right_x = np.sum(right_gaussianApprox * xs) / np.sum(right_gaussianApprox)
    right_y = np.mean(right_gaussianApprox)

    axs[0].plot(xs, left_gaussianApprox, label="Left", color="blue")
    axs[0].fill_between(xs, np.zeros(len(bins)*3), left_gaussianApprox, alpha=alpha, color="blue")
    axs[0].scatter(left_x, left_y, color="blue", marker="1", s=(30*72./fig.dpi)**2)
    axs[0].plot(xs, right_gaussianApprox, label="Right", color="green")
    axs[0].fill_between(xs, np.zeros(len(bins)*3), right_gaussianApprox, alpha=alpha, color="green")
    axs[0].scatter(right_x, right_y, color="green", marker="2", s=(30*72./fig.dpi)**2)
    axs[0].legend()
    axs[0].set_title(metric + " - " + left_mask)
    axs[0].grid(linestyle = '--', linewidth = 0.5)
    axs[0].set_xlabel("Pixel values")
    axs[0].set_ylabel("Frequency")

    diff = (left_gaussianApprox - right_gaussianApprox) / np.max([np.max(left_gaussianApprox ), np.max(right_gaussianApprox)])
    #sigmoid = fitSigmoid(xs, diff)
    #sigmoidEval = evaluateSigmoid(xs, sigmoid)
    #area = round(np.mean(sigmoidEval), 4)
    axs[1].plot(xs, diff, label="Left - Right", color="magenta", alpha = 2 * alpha)
    #axs[1].plot(xs, sigmoidEval, label="Eval" + str(area), color="red")
    axs[1].fill_between(xs, np.zeros(len(bins)*3), diff, alpha=alpha, color="magenta")
    axs[1].legend()
    axs[1].set_title(metric + " - " + left_mask)
    axs[1].grid(linestyle = '--', linewidth = 0.5)
    axs[1].set_xlabel("Pixel values")
    axs[1].set_ylabel("Difference between contralerals")


    axs[2].plot(xs, np.cumsum(diff), color="red", alpha=2*alpha)
    axs[2].fill_between(xs, np.zeros(len(bins)*3), np.cumsum(diff), alpha=alpha, color="red")
    axs[2].set_title(metric + " - " + left_mask)
    axs[2].grid(linestyle = '--', linewidth = 0.5)
    axs[2].set_xlabel("Pixel values")
    axs[2].set_ylabel("Cumulative area")

    fig.savefig(save_path)
    plt.close()

    print(f"{metric} - {left_mask} plot done.")

def meanContralateralComparison(study_path, patient_list, metric, injured_masks, healthy_masks, mask_name, plot_folder_path=None):
    alpha = 1/3
    if plot_folder_path is None:
        plot_folder_path = pjoin(study_path, "Contralateral_Comparison_Plots", "Mean")
    if not exists(plot_folder_path):
        makedirs(plot_folder_path)
    save_path = pjoin(plot_folder_path, metric + "_" + mask_name + ".svg")
    
    fig, axs = plt.subplots(3, 1, sharex=True)

    ref_bins = 0
    max_bin = 0
    for patient, injured_mask, healthy_mask in zip(patient_list, injured_masks, healthy_masks):
        injured_bins_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + injured_mask + "_bins.csv")
        injured_bins = np.loadtxt(injured_bins_path)
        healthy_bins_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + healthy_mask + "_bins.csv")
        healthy_bins = np.loadtxt(healthy_bins_path)
        if injured_bins[-1] > max_bin:
            max_bin = injured_bins[-1]
            ref_bins = np.copy(injured_bins)
        if healthy_bins[-1] > max_bin:
            max_bin = healthy_bins[-1]
            ref_bins = np.copy(healthy_bins)

    mean_injured = np.zeros(len(ref_bins) - 1)
    mean_healthy = np.zeros(len(ref_bins)- 1)
    for patient, injured_mask, healthy_mask in zip(patient_list, injured_masks, healthy_masks):
        injured_bins_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + injured_mask + "_bins.csv")
        injured_histogram_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + injured_mask + "_histo.csv")
        injured_bins = np.loadtxt(injured_bins_path)
        injured_histo = np.loadtxt(injured_histogram_path)
        mean_injured += interpolateHisto(injured_bins, injured_histo, ref_bins)
        healthy_bins_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + healthy_mask + "_bins.csv")
        healthy_histogram_path = pjoin(study_path, patient, "histograms", patient + "_" + metric + "_" + healthy_mask + "_histo.csv")
        healthy_bins = np.loadtxt(healthy_bins_path)
        healthy_histo = np.loadtxt(healthy_histogram_path)
        mean_healthy += interpolateHisto(healthy_bins, healthy_histo, ref_bins)

    xs = np.linspace(ref_bins[0], ref_bins[-1], len(ref_bins)*3)
    injured_gauss = fitNGauss(ref_bins, mean_injured, 25)
    injured_gaussianApprox = evaluateNGauss(xs, injured_gauss)
    healthy_gauss = fitNGauss(ref_bins, mean_healthy, 25)
    healthy_gaussianApprox = evaluateNGauss(xs, healthy_gauss)
    injured_x = np.sum(injured_gaussianApprox * xs) / np.sum(injured_gaussianApprox)
    injured_y = np.mean(injured_gaussianApprox)
    healthy_x = np.sum(healthy_gaussianApprox * xs) / np.sum(healthy_gaussianApprox)
    healthy_y = np.mean(healthy_gaussianApprox)

    axs[0].plot(xs, injured_gaussianApprox, label="Injured", color="blue")
    axs[0].fill_between(xs, np.zeros(len(ref_bins)*3), injured_gaussianApprox, alpha=alpha, color="blue")
    axs[0].scatter(injured_x, injured_y, color="blue", marker="1", s=(30*72./fig.dpi)**2)
    axs[0].plot(xs, healthy_gaussianApprox, label="Right", color="green")
    axs[0].fill_between(xs, np.zeros(len(ref_bins)*3), healthy_gaussianApprox, alpha=alpha, color="green")
    axs[0].scatter(healthy_x, healthy_y, color="green", marker="2", s=(30*72./fig.dpi)**2)
    axs[0].legend()
    axs[0].set_title(metric + " - " + mask_name)
    axs[0].grid(linestyle = '--', linewidth = 0.5)
    axs[0].set_ylabel("Frequency")

    diff = (injured_gaussianApprox - healthy_gaussianApprox) / np.max([np.max(injured_gaussianApprox ), np.max(healthy_gaussianApprox)])
    #sigmoid = fitSigmoid(xs, diff)
    #sigmoidEval = evaluateSigmoid(xs, sigmoid)
    #area = round(np.mean(sigmoidEval), 4)
    axs[1].plot(xs, diff, label="Injured - Healthy", color="magenta", alpha = 2 * alpha)
    #axs[1].plot(xs, sigmoidEval, label="Eval" + str(area), color="red")
    axs[1].fill_between(xs, np.zeros(len(ref_bins)*3), diff, alpha=alpha, color="magenta")
    axs[1].legend()
    axs[1].grid(linestyle = '--', linewidth = 0.5)
    axs[1].set_ylabel("Difference between contralerals")


    axs[2].plot(xs, np.cumsum(diff), color="red", alpha=2*alpha)
    axs[2].fill_between(xs, np.zeros(len(ref_bins)*3), np.cumsum(diff), alpha=alpha, color="red")
    axs[2].grid(linestyle = '--', linewidth = 0.5)
    axs[2].set_xlabel("Pixel values")
    axs[2].set_ylabel("Cumulative area")

    fig.savefig(save_path)
    plt.close()

    print(f"{metric} - {mask_name} plot done.")


def main(args):
    left_masks = ["xtract_prob_Corticospinal_Tract_L",
             "xtract_prob_Frontal_Aslant_Tract_L",
             "xtract_prob_Superior_Longitudinal_Fasciculus_1_L",
             "xtract_prob_Superior_Longitudinal_Fasciculus_2_L",
             "xtract_prob_Superior_Longitudinal_Fasciculus_3_L",
             "xtract_prob_Superior_Thalamic_Radiation_L",
             "Cortico_Ponto_Cerebellum_Left",
             "harvardoxford-subcortical_prob_Left_Thalamus"]
    right_masks = ["xtract_prob_Corticospinal_Tract_R",
             "xtract_prob_Frontal_Aslant_Tract_R",
             "xtract_prob_Superior_Longitudinal_Fasciculus_1_R",
             "xtract_prob_Superior_Longitudinal_Fasciculus_2_R",
             "xtract_prob_Superior_Longitudinal_Fasciculus_3_R",
             "xtract_prob_Superior_Thalamic_Radiation_R",
             "Cortico_Ponto_Cerebellum_Right",
             "harvardoxford-subcortical_prob_Right_Thalamus"]

    mask_names = ["xtract_prob_Corticospinal_Tract",
             "xtract_prob_Frontal_Aslant_Tract",
             "xtract_prob_Superior_Longitudinal_Fasciculus_1",
             "xtract_prob_Superior_Longitudinal_Fasciculus_2",
             "xtract_prob_Superior_Longitudinal_Fasciculus_3",
             "xtract_prob_Superior_Thalamic_Radiation",
             "Cortico_Ponto_Cerebellum",
             "harvardoxford-subcortical_prob_Thalamus"]

    metrics = ['FA', "MD", "AD", "RD"]

    patient_list = ["20_01_01_E0", "20_01_01_E1", "20_01_01_E2", "20_11_02_E1"]

    study_path = "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\TestStudy"

    """for metric in metrics:
        for left_mask, right_mask in zip(left_masks, right_masks):
            for patient in patient_list:
                print(patient, metric, left_mask, right_mask)
                displayContralateral(study_path, patient, metric, left_mask, right_mask)"""

    for metric in metrics:
        for left_mask, right_mask, mask_name in zip(left_masks, right_masks, mask_names):
            injured_masks = [left_mask] * len(patient_list)
            healthy_masks = [right_mask] * len(patient_list)
            meanContralateralComparison(study_path, patient_list, metric, injured_masks, healthy_masks, mask_name)


if __name__ == "__main__":
    main(sys.argv)