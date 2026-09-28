import pandas as pd
import numpy as np
from POVM_extract_data import *
from POVM_fast_recursive_releation import *
from scipy.stats import poisson
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from plot_photon_statistics import *
import time



'''poissonian state'''
def poss_stat(mu, event):
    mu = mu[0]
    input_state = []
    for i in range(event):
        input_state.append(poisson.pmf(i,mu))
    input_state.append(poisson.sf(event-1,mu))
    return np.array(input_state)

def probability(eta, MM, m, n):
    """MM is the number of pixel"""
    """n photon sent"""
    """m photon detected"""
    """eta : efficiency of the array"""
    p = 0
    for j in range(m+1):
        p += (-1)**j * ( scipy.special.factorial(m) / ( scipy.special.factorial(j) * scipy.special.factorial(m-j) ) ) * ( 1- eta + (m-j)*eta/MM)**n
    pp = scipy.special.factorial(MM) / ( scipy.special.factorial(m) * scipy.special.factorial (MM-m)) * p
    return pp

def output_probability(mu, eta, M):
    """MM is the number of pixel"""
    """n photon sent"""
    """m photon detected"""
    """eta : efficiency of the array"""
    """mu is the mean photon number per pulse"""
    Q_theo = []
    
    for m in range(0,M+1):
        p = 0
        for n in range(m,M+1):
            # p += probability(eta, M, m, n) * (eta*mu)**n * math.e**(-eta*mu) / scipy.special.factorial(n)
            p += probability(eta, M, m, n) * mu**n * math.e**(-mu) / scipy.special.factorial(n)
        Q_theo.append(float(p))
    return np.array(Q_theo)

def povm_uniform_distibution(eta, M):
    P = np.zeros((M+1,M+1))
    for m in range(0,M+1): #photon I detect
        p = []
        for n in range(m,M+1): #photon I sent
            pp = probability(eta,M,m,n)
            P[m][n]= pp
    return P


def binomial_std(Q):
    std_Q = []
    Nt = sum(Q)
    Q = np.array(Q)/Nt
    for i in range(len(Q)):
        std_Q.append( (Q[i]*(1-Q[i])* Nt)**0.5 / Nt )
    return std_Q


def func_opt(mu, eta, M, Q):
    Q_theo  = output_probability(mu, eta, M)
    return np.array(Q) - Q_theo

def func_opt2(mu, event, reconstructed_input):
    input_state = poss_stat(mu, event)
    return input_state - reconstructed_input


'''Repetition rate of the laser in MHz'''
R = 1 * 10**6 

'''Acquisition time in minutes'''
t = 60

time_window = 2 #in ns

'''maximum photon-event''' 
M = 28
eta = 0.855
P = povm_uniform_distibution(eta, M)


"""SELECT CORRECT FOLDER""" # terminate the folder with / !!
path_base = "Z:/Projects/Detectors/SNSPD/Measurements/PNR_Certification_Experiment/ID900/Imax=8_dataset_3/gN2a/gN2a8-12/1180uA/"
folder = "mu=1.1574/"
path = path_base + folder
os.makedirs(os.path.dirname(path), exist_ok=True)

filename1 = "gN2a8-12_-1180uA_ZFL500LN+_1-2splitter_ZFL2500VH+_1-4splitter_2.7Vamp_60s_level1-4_2025-10-04T02_01_04"
filename2 = "gN2a8-12_-1180uA_ZFL500LN+_1-2splitter_ZFL2500VH+_1-4splitter_2.7Vamp_60s_level1-4_2025-10-04T02_03_05"
extension = ".csv"
path_complete1 = path + filename1 + extension
path_complete2 = path + filename2 + extension

if __name__ == "__main__":
    QQ_counts = extract_histograms_data_8_levels(path_complete1, path_complete2, R, t, time_window)
    print(QQ_counts)
    
    zero_left = [0]*(M+1-len(QQ_counts))
    QQ = list(QQ_counts) + zero_left
    std_QQ = binomial_std(QQ)
    Nt = sum(QQ)
    QQ = np.array(QQ)/Nt
    mu_0 = np.array([0.1])
    print("experimental Q")
    print(QQ)

    res = least_squares(func_opt, mu_0, args = (eta, M, QQ), verbose = 1, method='trf', ftol=3e-16, xtol=3e-16, gtol=3e-16, max_nfev = 10000)

    mu_fit = "%.2f" % round(float(res.x), 2)
    theoretical_Q = output_probability(res.x, eta, M)
    print("theoretical Q")
    print(theoretical_Q)

    print("mu fit")
    print(mu_fit)

    saving_folder = path + "Output Analysis" + "/" + "SDE = " + str(eta) + "/" #" - no 4th-ph/"
    os.makedirs(saving_folder, exist_ok=True)
    timestr = time.strftime("%Y%m%d-%H%M%S")
    comment = "log"
    file = open(f"{saving_folder}OUTPUT_analysis{comment}_{timestr}.txt", "w")
    file.write(f"POVM\n{P}\n\n")
    # file.write(f"STD_POVM\n{STD_POVM}\n\n")
    file.write(f"Experimental output counts = {QQ_counts}\n")
    # file.write(f"Experimental output statistics = {QQ}\nExperimental error = {std_QQ}\nComputed output = {theoretical_Q}\n\n")
    file.write(f"DIFFERENCE_output = {QQ-theoretical_Q}\nDIFFERENCE%_output = {100*(QQ-theoretical_Q)/QQ}\n\n")
    file.write(f"Fitted mu from output is {float(res.x)}\n\n")
    file.write(f"Time window used for histogram is {float(time_window)} ns\n\n")
    file.close()
    plot_statistics_output_log_FINAL(mu_fit, theoretical_Q, QQ, std_QQ, comment, saving_folder)
    plt.figure()
    comment = "linear"
    plot_statistics_output_linear_FINAL(mu_fit, theoretical_Q, QQ, std_QQ, comment, saving_folder)


    plt.figure()
    reconstructed_input = np.dot(np.linalg.inv(P), np.array(QQ))

    res = least_squares(func_opt2, mu_0, args = (M, reconstructed_input), verbose = 1, method='trf', ftol=3e-16, xtol=3e-16, gtol=3e-16, max_nfev = 10000)
    mu_fit = "%.2f" % round(float(res.x), 2)
    theoretical_input = poss_stat(res.x, M)
    saving_folder = path + "Input Analysis" + "/" + "SDE = " + str(eta) + "/" 
    os.makedirs(saving_folder, exist_ok=True)
    timestr = time.strftime("%Y%m%d-%H%M%S")
    comment = "linear"
    file = open(f"{saving_folder}Input_analysis{comment}_{timestr}.txt", "w")
    file.write(f"POVM\n{P}\n\n")
    # file.write(f"STD_POVM\n{STD_POVM}\n\n")
    file.write(f"Experimental output counts = {QQ_counts}\n")
    file.write(f"Experimental output statistics = {QQ}\nExperimental error = {std_QQ}\nComputed output = {theoretical_Q}\n\n")
    file.write(f"Reconstructed Input statistics = {reconstructed_input}\n\n")
    file.write(f"Theoretical statistics = {theoretical_input}\n\n")
    file.write(f"DIFFERENCE_input = {reconstructed_input-theoretical_input}\nDIFFERENCE%_output = {100*(reconstructed_input-theoretical_input)/reconstructed_input}\n\n")
    file.write(f"Fitted mu from input is {float(res.x)}\n\n")
    file.write(f"Time window used for histogram is {float(time_window)} ns\n\n")
    file.close()
    plot_statistics_input_log_FINAL(mu_fit, theoretical_input, reconstructed_input, std_QQ, comment, saving_folder)
    plt.figure()
    comment = "linear"
    plot_statistics_input_linear_FINAL(mu_fit, theoretical_input, reconstructed_input, std_QQ, comment, saving_folder)
    plt.figure()
    plot_statistics_input_linear_blue_FINAL(mu_fit, theoretical_input, reconstructed_input, std_QQ, comment, saving_folder)