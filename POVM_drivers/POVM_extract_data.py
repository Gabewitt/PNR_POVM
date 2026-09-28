import pandas as pd
import os

def extract_raw_data(path, R, t, time_window):
    Nt = R * t
    timebin = 100 #in ps
    timebin_from_time_window = time_window*1000/timebin
    cols_list = ["BinCenter (ps)", "Histo01", "Histo02", "Histo03", "Histo04"]
    data =  pd.read_csv(path, sep = ";", usecols=cols_list)
    
    '''Import data from excel file'''
    h1 = data["Histo01"].tolist()
    h2 = data["Histo02"].tolist()
    h3 = data["Histo03"].tolist()
    h4 = data["Histo04"].tolist()
    max_value_1 = max(h1)
    max_index_1 = h1.index(max_value_1)
    n1=0
    for i in range(max_index_1-int(timebin_from_time_window/2), max_index_1+int(timebin_from_time_window/2)):
        n1 = n1 + h1[i]

    max_value_2 = max(h2)
    max_index_2 = h2.index(max_value_2)
    n2=0
    for i in range(max_index_2-int(timebin_from_time_window/2), max_index_2+int(timebin_from_time_window/2)):
        n2 = n2 + h2[i]

    max_value_3 = max(h3)
    max_index_3 = h3.index(max_value_3)
    n3=0
    for i in range(max_index_3-int(timebin_from_time_window/2), max_index_3+int(timebin_from_time_window/2)):
        n3 = n3 + h3[i]

    max_value_4 = max(h4)
    max_index_4 = h4.index(max_value_4)
    n4=0
    for i in range(max_index_4-int(timebin_from_time_window/2), max_index_4+int(timebin_from_time_window/2)):
        n4 = n4 + h4[i]

    '''Computing th final counts'''
    c4 = n4
    c3 = n3 - c4
    c2 = n2 - c3 - c4
    c1 = n1 - c2 - c3 - c4
    c0 = Nt-n1


    '''Creating final histogram'''
    hist = []
    hist.append(c0)
    hist.append(c1)
    hist.append(c2)
    hist.append(c3)
    hist.append(c4)

    return hist


def extract_raw_thermal_data(path, R, t):
    Nt = R * t

    cols_list = ["BinCenter (ps)", "Histo01", "Histo02", "Histo03", "Histo04"]
    data =  pd.read_csv(path, sep = ";", usecols=cols_list)
    
    '''Import data from excel file'''
    h1 = data["Histo01"].tolist()
    h2 = data["Histo02"].tolist()
    h3 = data["Histo03"].tolist()
    h4 = data["Histo04"].tolist()
    max_value_1 = max(h1)
    max_index_1 = h1.index(max_value_1)
    n1=0
    for i in range(max_index_1-3, max_index_1+4):
        n1 = n1 + h1[i]


    max_value_2 = max(h2)
    max_index_2 = h2.index(max_value_2)
    n2=0
    for i in range(max_index_2-10, max_index_2+10):
        n2 = n2 + h2[i]

    max_value_3 = max(h3)
    max_index_3 = h3.index(max_value_3)
    n3=0
    for i in range(max_index_3-10, max_index_3+10):
        n3 = n3 + h3[i]


    max_value_4 = max(h4)
    max_index_4 = h4.index(max_value_4)
    n4=0
    for i in range(max_index_4-10, max_index_4+10):
        n4 = n4 + h4[i]


    '''Computing th final counts'''
    c4 = n4
    c3 = n3 - c4
    c2 = n2 - c3 - c4
    c1 = n1 - c2 - c3 - c4
    c0 = Nt-n1


    '''Creating final histogram'''
    hist = []
    hist.append(c0)
    hist.append(c1)
    hist.append(c2)
    hist.append(c3)
    hist.append(c4)


    return hist


def counts_from_histo(histo, timebin_from_time_window):
    max_value = max(histo)
    max_index_ = histo.index(max_value)
    n = 0
    for i in range(max_index_-int(timebin_from_time_window/2), max_index_+int(timebin_from_time_window/2)):
        n = n + histo[i]
    return n


def extract_histograms_data_8_levels(path1, path2, R, t, time_window):
    Nt = R * t
    timebin = 100 #in ps
    timebin_from_time_window = time_window*1000/timebin
    cols_list = ["BinCenter (ps)", "Histo01", "Histo02", "Histo03", "Histo04"]
    data1 =  pd.read_csv(path1, sep = ";", usecols=cols_list)
    data2 =  pd.read_csv(path2, sep = ";", usecols=cols_list)
    
    '''Import data from excel file'''
    h1 = data1["Histo01"].tolist()
    h2 = data1["Histo02"].tolist()
    h3 = data1["Histo03"].tolist()
    h4 = data1["Histo04"].tolist()

    h5 = data2["Histo01"].tolist()
    h6 = data2["Histo02"].tolist()
    h7 = data2["Histo03"].tolist()
    h8 = data2["Histo04"].tolist()
    histo_list = [h1, h2, h3, h4, h5, h6, h7, h8]
    # max_value_1 = max(h1)
    # max_index_1 = h1.index(max_value_1)
    # n1=0
    # for i in range(max_index_1-int(timebin_from_time_window/2), max_index_1+int(timebin_from_time_window/2)):
    #     n1 = n1 + h1[i]

    # max_value_2 = max(h2)
    # max_index_2 = h2.index(max_value_2)
    # n2=0
    # for i in range(max_index_2-int(timebin_from_time_window/2), max_index_2+int(timebin_from_time_window/2)):
    #     n2 = n2 + h2[i]

    # max_value_3 = max(h3)
    # max_index_3 = h3.index(max_value_3)
    # n3=0
    # for i in range(max_index_3-int(timebin_from_time_window/2), max_index_3+int(timebin_from_time_window/2)):
    #     n3 = n3 + h3[i]

    # max_value_4 = max(h4)
    # max_index_4 = h4.index(max_value_4)
    # n4=0
    # for i in range(max_index_4-int(timebin_from_time_window/2), max_index_4+int(timebin_from_time_window/2)):
    #     n4 = n4 + h4[i]

    n_list = []
    for i in range(len(histo_list)):
        n = counts_from_histo(histo_list[i], timebin_from_time_window)
        n_list.append(n)



    '''Computing th final counts'''
    c8 = n_list[-1]
    c7 = n_list[-2] - n_list[-1]
    c6 = n_list[-3] - n_list[-2]
    c5 = n_list[-4] - n_list[-3]
    c4 = n_list[-5] - n_list[-4]
    c3 = n_list[-6] - n_list[-5]
    c2 = n_list[-7] - n_list[-6]
    c1 = n_list[-8] - n_list[-7]
    c0 = Nt-n_list[-8]


    '''Creating final histogram'''
    hist = []
    hist.append(c0)
    hist.append(c1)
    hist.append(c2)
    hist.append(c3)
    hist.append(c4)
    hist.append(c5)
    hist.append(c6)
    hist.append(c7)
    hist.append(c8)

    return hist


if __name__ == "__main__":
    path = "05_PNR/8-10/-1200uA/8 level measurement with 1 ID900/40MHz/20240524-PNR-8-10_T1k=0.94_12-31_R=40MHz_mu=0.5/"
    filename1 = "gN2a8-10_-1200uA_40MHz_mu=0.5_1-4levels_2.7Vamp_30s_2024-05-24T12_35_09"
    filename2 = "gN2a8-10_-1200uA_40MHz_mu=0.5_5-8levels_2.7Vamp_30s_2024-05-24T12_36_53"
    extension = ".csv"
    path_complete1 = path + filename1 + extension
    path_complete2 = path + filename2 + extension
    R = 40*10**6
    t = 30
    time_window = 2 #ns
    timebin = 100 #in ps

    histo = extract_histograms_data_8_levels(path_complete1, path_complete2, R, t, time_window)
    print(histo)
    # complete_name = os.path.join(path, short_filename + ".txt")
    # file=open(complete_name,'w')

    # file.write("Experimental output = " + str(histo) + '\n' + '\n' + "Sum(Histogram) = " + str(sum(histo)) + '\n' + '\n')