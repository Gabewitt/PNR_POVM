import shutil
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import time
from datetime import datetime
import os
import sys



# data_path = '/Volumes/RFsoc/ZCU208B/Data/20260512/PNR_rfsoc/Matterhorn/gN2a/8-10/1/1856/received_data_RAM_1.npy'

data_path = '/Volumes/RFsoc/ZCU208B/Data/20260513/PNR_rfsoc/Matterhorn/gN2a/8-10/3/1419/received_data_RAM_3.npy'

# data_path = '/Volumes/RFsoc/ZCU208B/Data/20260513/PNR_rfsoc/Matterhorn/gN2a/8-10/2/1520/received_data_RAM_2.npy'

# data_path = '/Volumes/RFsoc/ZCU208B/Data/20260513/PNR_rfsoc/Matterhorn/gN2a/8-10/0.5/1604/received_data_RAM_0.5.npy'

# data_path = '/Volumes/RFsoc/ZCU208B/Data/20260513/PNR_rfsoc/Matterhorn/gN2a/8-10/0.1/1648/received_data_RAM_0.1.npy'


# data_path = '/Volumes/RFsoc/ZCU208B/Data/20260508/PNR_rfsoc/Matterhorn/gN2a/8-10/1/received_data_RAM_1.npy'

mu = 3




saving_path = f'/Volumes/Wittenberg/PNR_POVM_Data/Waveforms_RFSoC/mu={mu}/'
os.makedirs(saving_path, exist_ok=True)


data = np.load(data_path).reshape(-1, 16)

# windowed_data = []

# for waveform in data:
#     max_idx = np.argmax(waveform)
#     start = max(0, max_idx - 4)
#     end = min(16, max_idx + 5)  
#     window = waveform[start:end]
    # windowed_data.append(window)
    


max_amp_linear = data.max(axis=1)
np.save(saving_path + f'max_amp_mu_{mu}.npy', max_amp_linear)


data_sq = data.astype(np.float64) ** 2


area = np.trapz(data_sq, axis=1)
# max_amp_sq = data_sq.max(axis=1)

# area_windowed = np.array([np.trapz(w.astype(np.float64)**2) for w in windowed_data])


np.save(saving_path + f'area_mu_{mu}.npy', area)
# np.save(saving_path + f'max_amp_sq_mu_{mu}.npy', max_amp_sq)
# np.save(saving_path + f'area_windowed_mu_{mu}.npy', area_windowed)


