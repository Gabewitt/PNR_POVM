import shutil
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import time
from datetime import datetime
import os
import sys


mu = 0.1
try_number = 10

data_path_area = f'/Volumes/Wittenberg/PNR_POVM_Data/Waveforms_RFSoC/mu={mu}/area_mu_{mu}.npy'
data_path_max_amp = f'/Volumes/Wittenberg/PNR_POVM_Data/Waveforms_RFSoC/mu={mu}/max_amp_mu_{mu}.npy'
area = np.load(data_path_area)
max_amp = np.load(data_path_max_amp)

saving_path = f'/Volumes/Wittenberg/PNR_POVM_Results/RFSoC/RAM_analysis/mu={mu}/v{try_number}/'
os.makedirs(saving_path, exist_ok=True)


#area_threshold = [6.014 51e7, 1.63105e8, 3.27852e8, 5.62401e8, 8.63422e8, 1.2551e9, 1.72303e9, 2.27388e9, 2.85853e9] #v1
# area_threshold = [1.32906e7, 6.01451e7, 1.63105e8, 3.27852e8, 5.62401e8, 8.63422e8, 1.2551e9, 1.72303e9, 2.27388e9, 2.85853e9] #v2 splitting the weird peaks
# area_threshold = [5.91451e7, 1.93105e8, 3.57852e8, 5.62401e8, 8.63422e8, 1.2551e9, 1.72303e9, 2.27388e9, 2.85853e9] #v3
#area_threshold = [5.71451e7, 2.13105e8, 3.87852e8, 5.72401e8, 8.63422e8, 1.2551e9, 1.72303e9, 2.27388e9, 2.85853e9] #v4
max_amp_threshold = [4664, 7863, 11056, 14501, 18009, 21717, 25440, 29297, 32619]
# max_amp_threshold = [2560, 4664, 7863, 11056, 14501, 18009, 21717, 25440, 29297, 32619] #v2 

# area_threshold = [10**7, 5.71451e7, 2.13105e8, 3.87852e8, 5.72401e8, 8.63422e8, 1.2551e9, 1.72303e9, 2.27388e9, 2.85853e9] #v5 getting rid of first weird peak through blocking in single photon peak from mu 0.1
# area_threshold = [1.32906e7, 5.71451e7, 2.13105e8, 3.87852e8, 5.72401e8, 8.63422e8, 1.2551e9, 1.72303e9, 2.27388e9, 2.85853e9] #v6 Getting rid of first weird peak through the valley
# max_amp_threshold = [2560, 4664, 7863, 11056, 14501, 18009, 21717, 25440, 29297, 32619] #v7

# max_amp_threshold = [2035, 4664, 7863, 11056, 14501, 18009, 21717, 25440, 29297, 32619] #v8 by looking at where the real single photon peak is from mu = 0.1
# area_threshold = [9.57511e6, 5.71451e7, 1.63105e8, 3.27852e8, 5.62401e8, 8.63422e8, 1.2551e9, 1.72303e9, 2.27388e9, 2.85853e9] #v8 more precise single photon peak from 0.1
# area_threshold = [1.32906e7, 5.71451e7, 1.63105e8, 3.27852e8, 5.62401e8, 8.63422e8, 1.2551e9, 1.72303e9, 2.27388e9, 2.85853e9] #v9 same as v6 but better higher thresholds
# max_amp_threshold = [2250, 4664, 7863, 11056, 14501, 18009, 21717, 25440, 29297, 32619] #v9 proper seperation of single photon peak from weird peak by looking at mu = 0.1 histogram
area_threshold = [5.91451e7, 1.63105e8, 3.27852e8, 5.62401e8, 8.63422e8, 1.2551e9, 1.72303e9, 2.27388e9, 2.85853e9]

area_bins = np.logspace(np.log10(area.min()), np.log10(area.max()), 5000)

plt.figure(dpi=300)
plt.hist(area, bins=area_bins)
plt.xlabel('Area')
plt.ylabel('Counts')
plt.title(f'Histogram of Areas for mu={mu}')
plt.grid()
plt.xscale('log')
if len(area_threshold) > 0:
    for i in range(len(area_threshold)):
        plt.axvline(area_threshold[i], color='r', linestyle='--', linewidth=0.2, label=f'Threshold {area_threshold[i]}')
    plt.legend(loc='upper right', fontsize='xx-small')
plt.savefig(saving_path + f'area_histogram_mu_{mu}.png', dpi=300)
plt.show()


max_amp_bins = np.logspace(np.log10(max_amp.min()), np.log10(max_amp.max()), 400)

plt.figure(dpi=300)
plt.hist(max_amp, bins=max_amp_bins)
plt.xlabel('Max Amplitude')
plt.ylabel('Counts')
plt.title(f'Histogram of Max Amplitudes for mu={mu}')
plt.grid()
plt.xscale('log')
if len(max_amp_threshold) > 0:
    for i in range(len(max_amp_threshold)):
        plt.axvline(max_amp_threshold[i], color='r', linestyle='--', linewidth=0.2, label=f'Threshold {max_amp_threshold[i]}')
    plt.legend(loc='upper right', fontsize='xx-small')
plt.savefig(saving_path + f'max_amp_histogram_mu_{mu}.png', dpi=300)
plt.show()

if len(area_threshold) > 0:
    photon_number_area = np.digitize(area, area_threshold) + 1
    np.save(saving_path + f'photon_number_area_mu_{mu}.npy', photon_number_area)
    np.save(saving_path + f'photon_number_area_counts_mu_{mu}.npy', np.bincount(photon_number_area))
    
if len(max_amp_threshold) > 0:
    photon_number_max_amp = np.digitize(max_amp, max_amp_threshold) + 1
    np.save(saving_path + f'photon_number_max_amp_mu_{mu}.npy', photon_number_max_amp)
    np.save(saving_path + f'photon_number_max_amp_counts_mu_{mu}.npy', np.bincount(photon_number_max_amp))
    

    

