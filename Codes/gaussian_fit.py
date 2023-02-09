import numpy as np
import matplotlib.pyplot as plt
from os.path import join as pjoin
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint

def evaluateNGauss(x, gauss):
    if type(x) == np.ndarray:
        value = np.zeros(x.shape)
    else:
        value = 0
    for (A, mu, sigma) in gauss:
        value += A * np.exp(- (1 / 2) * ((x - mu) / sigma)**2)
    return value


def objective_function(As, mus, sigmas, xs, histo):
    gauss = []
    for A, mu, sigma in zip(As, mus, sigmas):
        gauss.append((A, mu, sigma))
    gaussEval = evaluateNGauss(xs, gauss)
    RMSE = (np.sum((histo - gaussEval)**2) / len(xs))**(1/2)
    return RMSE



def fitNGauss(bins, histo, NGauss, init=None):
    xs = bins[:-1] + (bins[1]-bins[0])/2
    if init is None:
        mus = np.linspace(xs[0], xs[-1], num=NGauss)
        As = np.ones(NGauss) * np.mean(histo)
        if NGauss > 1:
            sigmas = np.ones(NGauss) * (mus[1] - mus[0]) / (2 * np.sqrt(2 * np.log(2)))
        else:
            sigmas = [(xs[-1] - xs[0])/2]
    else:
        As = init[0]
        mus = init[1]
        sigmas = init[2]
    mat = np.eye(len(As))
    ub = np.inf
    lb = 0
    cons = LinearConstraint(mat, lb, ub)
    res = minimize(objective_function, As, args=(mus, sigmas, xs, histo), constraints=cons)
    print(res.success)
    gauss = []
    for A, mu, sigma in zip(res.x, mus, sigmas):
        gauss.append((A, mu, sigma))
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


def fitNGaussFromPath(bins_path, histo_path, NGauss, init=None):

    histo = np.loadtxt(histo_path)
    bins = np.loadtxt(bins_path)

    return fitNGauss(bins, histo, NGauss, init=init), bins, histo


def main():
    bins_path = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\histograms\\20_01_01_E0_RD_xtract_prob_Corticospinal_Tract_L_bins.csv"
    histo_path = "C:\\Users\\Louis Lovat\\Desktop\\Memoire\\study\\20_01_01_E0\\histograms\\20_01_01_E0_RD_xtract_prob_Corticospinal_Tract_L_histo.csv"

    gauss, bins, histo = fitNGaussFromPath(bins_path, histo_path, 20)
    print(gauss)
    displayGaussianApprox(gauss, bins, histo, label="20_01_01_E0", title="RD - xtract_prob_Corticospinal_Tract_L")

if __name__ == "__main__":
    main()
