import numpy as np
from skimage.io import imread
from skimage.color import rgb2gray
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

content = rgb2gray(imread("content.png")[:,:,:3])
mask = imread("mask.png")[:,:,3]

dft = np.fft.fft2(mask)
dft = np.fft.fftshift(dft)
dft[:dft.shape[0] // 2 - 10, :] = np.zeros((dft.shape[0] // 2 - 10, dft.shape[1]))
dft[dft.shape[0] // 2 + 10:, :] = np.zeros((dft.shape[0] // 2 - 10, dft.shape[1]))
dft[: ,:dft.shape[1] // 2 - 10] = np.zeros((dft.shape[0], dft.shape[1] // 2 - 10))
dft[: ,dft.shape[1] // 2 + 10:] = np.zeros((dft.shape[0], dft.shape[1] // 2 - 10))
dft = np.fft.ifftshift(dft)
mask = abs(np.fft.ifft2(dft))
mask = mask / np.max(mask)

dft = np.fft.fft2(content)
dft = np.fft.fftshift(dft)
dft[:dft.shape[0] // 2 - 10, :] = np.zeros((dft.shape[0] // 2 - 10, dft.shape[1]))
dft[dft.shape[0] // 2 + 10:, :] = np.zeros((dft.shape[0] // 2 - 10, dft.shape[1]))
dft[: ,:dft.shape[1] // 2 - 10] = np.zeros((dft.shape[0], dft.shape[1] // 2 - 10))
dft[: ,dft.shape[1] // 2 + 10:] = np.zeros((dft.shape[0], dft.shape[1] // 2 - 10))
dft = np.fft.ifftshift(dft)
content = abs(np.fft.ifft2(dft))
content = content / np.max(content)

nBins = 200
bins = np.linspace(0, 1, nBins+1)
histoPoints = np.zeros(nBins)
for i in range(nBins-1):
    indices1 = bins[i]<= content 
    indices2 = content < bins[i+1]
    indices = indices1 * indices2
    histoPoints[i] = np.sum(mask[indices])
indices1 = bins[nBins-1]<= content 
indices = indices1 * indices2
histoPoints[nBins-1] = np.sum(mask[indices])

histoPoints = histoPoints / np.sum(histoPoints)
plt.bar(bins[:-1], histoPoints, width=bins[1]-bins[0], align="edge")
plt.plot(bins[:-1], savgol_filter(histoPoints, nBins, 30), color="red")
plt.show()

"""fig, axs = plt.subplots(1, 3)
axs[0].imshow(content, cmap="gray")
axs[1].imshow(mask, cmap="gray")
axs[2].imshow(content * mask, cmap="gray")
plt.show()"""
