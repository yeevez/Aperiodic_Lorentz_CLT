"""
template for processing the data from various simulations. the functions defined here will calculate the displacement
from the raw sample positions, as well as sample statistics on the whole sample based on various normalization functions. 
The names of the files from the simulation run need to be specified, as well as the size of the samples and the
number of collisions they ran for. calculating normalization depends on both the sampling rate and the number of collisions.
For example, if we have collected 1 million collision samples and sampled the position at every 1000 collisions (our current
plan at the time that I write this), then we need to normalized the ith position in a sample by the square root of 1000*i
(for diffusive normalization as an example). This is how the template code will normalize samples, but if we start saving
samples with a diferent approach down the line, remember to alter the code accordingly. 

if we finalize a test for Brownian motion upon which to apply to our samples, then we can include the code for that test in
this template
"""

import numpy as np
from scipy.io import loadmat


N_SAMPLES = 10000
radius = 0.3
max_bounces = 1000000
R = 3
bordersize = 2
n_positions_recorded = 1001
sampling_rate = 1000


def sample_displacement(n_samples, norm="root t"):
    #sometimes a sample number will not be there due to a mixup in the renaming scheme
    #so here we find the actual number of samples and include the try except statements
    samples = np.zeros((n_samples,n_positions_recorded,2))
    n_samples_actual = 0
    for i in range(n_samples):
        filename = 'C:/path/to/directory/sample_name' + str(i) + ".mat"
        try:
            sample = loadmat(filename)
            sample = sample['paths']
            samples[i] = sample
            n_samples_actual += 1
        except FileNotFoundError:
            continue

    samples = samples[:n_samples_actual]
    sample_displacement = samples.copy()
    for i in range(n_samples_actual):
        for j in range(n_positions_recorded):
            sample_displacement[i][j] -= sample_displacement[i][0]
    if norm == 'none':
        pass
    elif norm == 'root t':
        for j in range(n_positions_recorded):
            sample_displacement[:,j,:] /= np.sqrt(sampling_rate*j)
    elif norm == 'root tlogt':
        for j in range(n_positions_recorded):
            sample_displacement[:,j,:] /= np.sqrt(sampling_rate*j*np.log(sampling_rate*j))
    return sample_displacement

def sample_displacement_mean_var(sample_displacement,index):
    #index specifies the point along the particle trajectory we wish to calculate the mean and covariance. for instance
    #for a 5 million collision sample, with sampling rate 1000 we might call this function at index = 1000 and index = 5000
    #to find the statistics of the sample at 1 million and 5 million collisions respectively, and compare the two 
    return sample_displacement[:,index,:].mean(axis=0,),np.cov(sample_displacement[:,index,:],rowvar = False)

def max_flight_time(n_samples):
    max_time = 0
    for i in range(1,n_samples):
        filename = 'C:/path/to/directory/2d_circular_obstacle_5milcollisions_matrixA_sample_' + str(i) + "_flight_time.npy"
        try:
            flight_time = np.load(filename)
        except FileNotFoundError:
            continue
        else:
            if flight_time[0] > max_time:
                max_time = flight_time[0]

    return max_time

def growth_rate(mean_displacements):
    #returns the ratio of mean at the ith sample with mean at the (i-1)th sample
    growth_rate = np.zeros(mean_displacements.shape)
    for i in range(1,growth_rate.shape[0]):
            growth_rate[i] = mean_displacements[i]/mean_displacements[i-1]
    return growth_rate



index = 1000
sample_unnormed = sample_displacement(N_SAMPLES,norm='none')
sample_normed = sample_displacement(N_SAMPLES)
sample_supernormed = sample_displacement(N_SAMPLES,norm='root tlogt')
print("Unnormed Sample Statistics")
mean, var = sample_displacement_mean_var(sample_unnormed,index)
print("Mean",mean)
print("Covariance", var)
print("Diffusively normed Sample Statistics")
mean, var = sample_displacement_mean_var(sample_normed,index)
print("Mean",mean)
print("Covariance", var)
print("Superdiffusively normed Sample Statistics")
mean, var = sample_displacement_mean_var(sample_supernormed,index)
print("Mean",mean)
print("Covariance", var)
