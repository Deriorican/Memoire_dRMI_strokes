import numpy as np
import matplotlib.pyplot as plt
from os.path import join as pjoin
import sys

def evaluateNGauss(x, gauss):
    if type(x) == np.ndarray:
        value = np.zeros(x.shape)
    else:
        value = 0
    for (A, mu, sigma) in gauss:
        value += A * np.exp(- (1 / 2) * ((x - mu) / sigma)**2)
    return value


def fitNGauss(bins, histo, NGauss):
    ratio = 1.65
    xs = bins[:-1] + (bins[1]-bins[0])/2
    mus = np.linspace(xs[0], xs[-1], num=NGauss)
    if NGauss > 1:
        s = (mus[1] - mus[0]) / (2 * np.sqrt(2 * np.log(ratio)))
    else:
        sigmas = (xs[-1] - xs[0])/2
    B = np.zeros(NGauss)
    M = np.zeros((NGauss, NGauss))
    for a in range(NGauss):
        B[a] = np.sum(histo * np.exp(-1/2 * ((xs - mus[a]) / s) ** 2))
        for i in range(NGauss):
            M[a, i] = np.sum(np.exp(-1/2 * ((xs - mus[i]) / s) ** 2) * np.exp(-1/2 * ((xs - mus[a]) / s) ** 2))

    As = np.linalg.solve(M, B)
    gauss = []
    for i in range(NGauss):
        gauss.append((As[i], mus[i], s))
    return gauss


def displayGaussianApprox(gauss, bins, histo, label="", title=""):

    fig, ax = plt.subplots()

    ax.bar(bins[:-1], histo, width=bins[1]-bins[0], align="edge", alpha=1/3, label=label + " raw")
    xs = np.linspace(bins[0], bins[-1], num = 3*len(bins))
    ax.plot(xs, evaluateNGauss(xs, gauss), label=label + " smoothed")
    ax.grid(linestyle = '--', linewidth = 0.5)
    ax.set_xlabel("Pixel values")
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    ax.legend()

    plt.show()


def fitNGaussFromPath(bins_path, histo_path, NGauss):

    histo = np.loadtxt(histo_path)
    bins = np.loadtxt(bins_path)

    return fitNGauss(bins, histo, NGauss), bins, histo


def main():
    bins_path = "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\TestStudy\\20_01_01_E0\\histograms\\20_01_01_E0_FA_xtract_prob_Corticospinal_Tract_L_bins.csv"
    histo_path = "C:\\Users\\Louis Lovat\\Desktop\\Memoire_dRMI_strokes\\TestStudy\\20_01_01_E0\\histograms\\20_01_01_E0_FA_xtract_prob_Corticospinal_Tract_L_histo.csv"

    gauss, bins, histo = fitNGaussFromPath(bins_path, histo_path, 20)
    print(gauss)
    displayGaussianApprox(gauss, bins, histo, label="20_01_01_E0", title="RD - xtract_prob_Corticospinal_Tract_L")

if __name__ == "__main__":
    main()
