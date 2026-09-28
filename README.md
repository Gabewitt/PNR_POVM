# README - PNR POVM Master's Project

Author: Gabriel Wittenberg (gabwitt123@gmail.com)


## Project Description

This project looked at the implementation of an RFSoC FPGA as a new acquisition device to be used in deployable quantum experiments.

This repository contains all the code used to extract, and analyse the data from the FPGA to perform a state reconstruction experiment (POVM). Included in this repository:

- Receiving data from FPGA code through Socket Library and saving as .npy
- The code necessary to find and save the PNR thresholds for threshold, maximum amplitude and area under the curve discrimination.
- The script to sort and classify the detections depending on discrimination type
- The script to perform the state reconstruction of the incoming light



## Receiving

[Receiving Script](./Receive_Data.py)

Requirement: 

- Socket sending function on the FPGA
- Socket Receiving script on the receiving system

Make sure the **HOST**, the **PORT**, the **points_per_waveform** and the **amount_of_waveforms** are identified and identical on the sending and receiving end.

Run the Receiving script on the receiving system, it will idle and wait for the sending script on the FPGA to be run. 

After completion, it will save the data as a raw .npy file.

## Extracting Waveform Information

[Extracting Information](./Extract_Photon_information.py)

Inputs: 

- **received_data_RAM_{mean_photon_number}.npy** is the raw data saved by the above Receiving Script

It processes the data in the following way:

1. Loads the .npy
2. Reshapes the data to be a list of lists. Where each sublist contains 16 points corresponding to one waveform
3. Calculates the maximum amplitude of each waveform and saves it
4. Squares each waveform
5. Calculates the area underneath the curve of each squared waveform through the trapezium method and saves thos values of areas.

Outputs:

- **max_amp_mu{mu}.npy** containing a list of maximum amplitudes corresponding to the max amp of each recorded waveform
- **area_mu_{mu}.npy** containing a list of areas corresponding to the area of the squared waveforms


## Dsicrminating Between Photon Levels

[Discrimination Script](./Plotting_Histograms_and_Thresholds.py)

Inputs:

- **area_mu_{mu}.npy** for the area discrimination
- **max_amp_mu{mu}.npy** for the max amplitude discrimination
- **max_amp_threshold = []** for the thresholds discrimnating from maximum amplitude
- **area_threshold = []** for the thresholds discriminating from area under the curve


Script:

The script will first plot and save the histograms of the photon distributions with their appropriate thresholds.
Then it will count the different photon counts depending on the given thresholds.

Outputs:

- **photon_number_max_amp_mu_{mu}.npy** uses np.digitze to assign a photon number to each max amplitude value depending on thresholds
- **photon_number_max_amp_counts_mu_{mu}.npy** uses np.bincount to count the amount of photon numbers

- **photon_number_area_mu_{mu}.npy** uses np.digitze to assign a photon number to each max amplitude value depending on thresholds
- **photon_number_area_counts_mu_{mu}.npy** uses np.bincount to count the amount of photon numbers


## POVM Measurment

[POVM Script](./POVM_output_input_analysis_RFSoC_histograms.py)

Inputs:

- **photon_number_max_amp_counts_mu_{mu}.npy**
or
- **photon_number_area_counts_mu_{mu}.npy**

Script will perform the POVM measurement of one data file. So you have to re run the script for every couunts file

It also calculates the zero photon counts.

Outputs:

- POVM State Reconstruction plots, both linear and log based.
- Reconstructed mean photon number estimate

[POVM Script for Threshold Discrimantion](./POVM_output_input_analysis_RFSoC.py)

This povm script is used for the threshold discrimantion type.

Inputs:

- **dps_corrected_thresholds_300s_1410_mu=3.csv** contains the discriminated counts from the threshold discrimination method from the RFSoC

Outputs:

- POVM State reconsctruction data

