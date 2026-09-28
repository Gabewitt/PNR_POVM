import shutil

import numpy as np
import sys
import matplotlib.pyplot as plt
from pathlib import Path
import time
from datetime import datetime
import os


sys.path.append("/home/rfsoc/github/RFSoC/PNR_measurement_scripts/drivers")

import receive_data_rfsoc as receiver
import find_thresholds as threshold_finder
import PNR_analysis as pnr_analyzer


now = datetime.now()
day = now.strftime("%Y%m%d")
s_live = now.strftime("%H%M")



device_name = "Matterhorn"
detector_name = "gN2a8-10"
laser_name = "Pulsed_Laser_Source"
mean_photon_number = 0.5
measurement_type = "Signal_analysis"



save_dir_data = "/home/rfsoc/Documents/SNSPD_Lab_Measurements/PNR/RFSoC/Data/" + f"{day}/" 
os.makedirs(save_dir_data, exist_ok=True)




HOST = "0.0.0.0"
PORT = 65432



points_per_waveform = 16
amount_of_waveforms = 10000000   # total waveforms you expect to receive

data = receiver.receive_data(
    host=HOST,
    port=PORT,
    points_per_waveform=points_per_waveform,
    amount_of_waveforms=amount_of_waveforms
)


save_path = save_dir_data +f"received_data_RAM_{mean_photon_number}.npy"
np.save(save_path, data)
print("Saved to:", save_path)