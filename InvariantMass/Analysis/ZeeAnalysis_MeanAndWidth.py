import matplotlib.pyplot as plt
import numpy as np

# open arrays 
mean = []
mean_error = []
sigma = []
sigma_error = []

layers = [4,5,6,7,8,9,10,11,12,13,14,16,17]

# loop over layers to extract fit results 
for layer in layers:
    # format in text file is [mean, error_mean, sigma, error_sigma]
    results = np.loadtxt(f'InvariantMassFit_{layer}layers_3ClustersBound.txt')
    mean.append(results[0])
    mean_error.append(results[1])
    sigma.append(results[2])
    sigma_error.append(results[3])

fig, (ax0, ax1) = plt.subplots(nrows=2,sharex=True,figsize=(10,15))

ax0.set_title("Plot of Mean as a Function of Longitudinal Layers")
ax0.errorbar(layers, mean, yerr=mean_error,fmt='o')
ax0.set_ylabel("Mean [GeV]")

ax1.set_title("Plot of Sigma as a Function of Longitudinal Layers")
ax1.errorbar(layers, sigma, yerr=sigma_error,fmt='o')
ax1.set_xlabel("Number of Longitudinal Layers")
ax1.set_ylabel("Sigma [GeV]")

plt.show()

fig.savefig("ZeeFit_MeanAndSigma_3Clusters_LowerBound.pdf")
