import matplotlib.pylab as plt
import pandas as pd
import numpy as np

data_dict = {}
pulse_lenght_1 = {}
pulse_lenght_3 = {}
pulse_lenght_10 = {}


########################################################################
"Short pulse, 1ns window"
###mu = 0.5 ###
reconstructed_05 = [0.592985, 0.309971, 0.080915, 0.014081, 0.001839, 0.000192, 0.000015, 0.000002, 0., 0., 0., 0., 0., 0., 0.]
error_reconstructed_05 = [0.011289, 0.013777, 0.006815, 0.001772, 0.00031,  0.000041, 0.000004, 0.000001, 0., 0., 0., 0., 0., 0., 0.]
theoretical_05 = [0.592973, 0.309891, 0.080975, 0.014106, 0.001843, 0.000193, 0.000017, 0.000001, 0., 0., 0., 0., 0., 0., 0.]
mu_05_recon = 0.52
pulse_lenght_1["mu = 0.5"] = [theoretical_05, reconstructed_05, error_reconstructed_05, mu_05_recon]


###mu = 1 ###
reconstructed_1 = [0.352458, 0.36768,  0.191594, 0.066605, 0.01732,  0.00361,  0.000629, 0.000091, 0.000012, 0.000001, 0., 0., 0., 0., 0.]
error_reconstructed_1 = [0.013563, 0.019512, 0.01706,  0.008752, 0.003076, 0.000821, 0.000177, 0.000032, 0.000005, 0.000001, 0., 0., 0., 0., 0.]
theoretical_1 = [0.352488, 0.367553, 0.191631, 0.066607, 0.017363, 0.003621, 0.000629, 0.000094, 0.000012, 0.000001, 0., 0., 0., 0., 0.]
mu_1_recon = 1.04
pulse_lenght_1["mu = 1"] = [theoretical_1, reconstructed_1, error_reconstructed_1, mu_1_recon]


### mu = 2 ###
reconstructed_2 = [0.12023,  0.254743, 0.269888, 0.190434, 0.100831, 0.042753, 0.015002, 0.004567, 0.001221, 0.000241, 0.000071, 0.000018, 0., 0., 0.]
error_reconstructed_2 = [0.010007, 0.023046, 0.036906, 0.043515, 0.037479, 0.024129, 0.01191,  0.004619, 0.001434, 0.000364, 0.000076, 0.000012, 0., 0., 0.]
theoretical_2 = [0.120273, 0.254738, 0.269766, 0.190453, 0.100845, 0.042717, 0.015079, 0.004563, 0.001208, 0.000284, 0.00006,  0.000012, 0.000002, 0., 0.]
mu_2_recon = 2.12
pulse_lenght_1["mu = 2"] = [theoretical_2, reconstructed_2, error_reconstructed_2, mu_2_recon]


pulse_lenght = 1
data_dict[f"{pulse_lenght} ns"] = pulse_lenght_1
########################################################################




def plot_statistics(mu_theo, pulse_lenght):
    dd = data_dict[f"{pulse_lenght} ns"]
    theo = dd[f"mu = {mu_theo}"][0]
    recon = dd[f"mu = {mu_theo}"][1]
    error = dd[f"mu = {mu_theo}"][2]
    mu = dd[f"mu = {mu_theo}"][3]
    plt.clf()
    plt.style.use('py_files/perso2.txt')
    x = [i for i in range(len(theo))]
    plt.bar(x, height=theo, width=0.2, bottom=None, align='center', color = "gold", edgecolor = "black", label = "Poissonian distribution")

    plt.errorbar(x, recon, yerr = error, fmt='ro', label = "Reconstructed statistics")
    plt.title(r"$\mu$" + f" = {mu} , Pulse lenght = {pulse_lenght} ns")
    plt.xlabel("Photon number")
    plt.ylabel("Probability")
    plt.legend()
    tp = max(recon) + max(error)
    plt.xlim(left = -0.3, right=7.3)
    plt.ylim(bottom=0, top=1.01*tp)
    plt.grid(False)
    comment = "_cutted_data_and_summed" #start with _ (underscore)
    plt.savefig(f"mu={mu}_pulse={pulse_lenght}ns{comment}.png", bbox_inches='tight', dpi = 1200, format='png',transparent=True)
    plt.savefig(f"mu={mu}_pulse={pulse_lenght}ns{comment}.pdf", bbox_inches='tight', dpi = 1200, format='pdf',transparent=True)

    plt.show()

def plot_statistics2(mu, pulse_lenght, theo, recon, error, comment, path):

    plt.clf()
    plt.style.use('perso2.txt')
    x = [i for i in range(len(theo))]
    plt.bar(x, height=theo, width=0.2, bottom=None, align='center', color = "gold", edgecolor = "black", label = "Poissonian distribution")

    plt.errorbar(x, recon, yerr = error, fmt='ro', label = "Reconstructed statistics")
    plt.title(r"$\mu$" + f" = {mu} , Pulse lenght = {pulse_lenght} ns")
    plt.xlabel("Photon number")
    plt.ylabel("Probability")
    plt.legend()
    tp = max(recon) + max(error)
    plt.xlim(left = -0.3, right=10.3)
    # plt.ylim(bottom=0, top=1.01*tp)
    plt.ylim(bottom=1e-6, top=1)
    plt.yscale('log')
    plt.grid(False)
    # comment = "_cutted_data_and_summed" #start with _ (underscore)
    plt.savefig(f"{path}mu={mu}_pulse={pulse_lenght}ns{comment}.png", bbox_inches='tight', dpi = 1200, format='png',transparent=True)
    plt.savefig(f"{path}mu={mu}_pulse={pulse_lenght}ns{comment}.pdf", bbox_inches='tight', dpi = 1200, format='pdf',transparent=True)

    plt.show()

def plot_statistics_output(mu, pulse_lenght, theo, recon, error, comment, path):

    plt.clf()
    plt.style.use('perso2.txt')
    x = [i for i in range(len(theo))]
    plt.bar(x, height=theo, width=0.2, bottom=None, align='center', color = "gold", edgecolor = "black", label = "Poissonian fit")

    plt.errorbar(x, recon, yerr = error, fmt='ro', label = "Detected statistics")
    plt.title(r"$\mu$" + f" = {mu} , Pulse lenght = {pulse_lenght} ns")
    plt.xlabel("Photon number")
    plt.ylabel("Probability")
    plt.legend()
    tp = max(recon) + max(error)
    plt.xlim(left = -0.3, right=9.3)
    # plt.ylim(bottom=0, top=1.05*tp)
    ####################
    plt.ylim(bottom=1e-6, top=1)
    plt.yscale('log')
    #####################
    plt.grid(False)
    # comment = "_cutted_data_and_summed" #start with _ (underscore)
    plt.savefig(f"{path}mu={mu}_pulse={pulse_lenght}ns{comment}.png", bbox_inches='tight', dpi = 1200, format='png',transparent=True)
    plt.savefig(f"{path}mu={mu}_pulse={pulse_lenght}ns{comment}.pdf", bbox_inches='tight', dpi = 1200, format='pdf',transparent=True)

    plt.show()


def plot_statistics_output_log_FINAL(mu, theo, recon, error, comment, path):

    plt.clf()
    plt.style.use('perso2.txt')
    x = [i for i in range(len(theo))]
    plt.bar(x, height=theo, width=0.2, bottom=None, align='center', color = "gold", edgecolor = "black", label = "Poissonian fit")

    plt.errorbar(x, recon, yerr = error, fmt='ro', label = "Detected statistics")
    # plt.title(r"$\mu$" + f" = {mu} , Pulse lenght = {pulse_lenght} ns")
    plt.xlabel("Photon number")
    plt.ylabel("Probability")
    plt.legend()
    tp = max(recon) + max(error)
    plt.xlim(left = -0.3, right=8.3)
    # plt.ylim(bottom=0, top=1.05*tp)
    ####################
    plt.ylim(bottom=1e-6, top=1)
    plt.yscale('log')
    #####################
    plt.grid(False)
    # comment = "_cutted_data_and_summed" #start with _ (underscore)
    plt.savefig(f"{path}mu={mu}_{comment}.png", bbox_inches='tight', dpi = 1200, format='png',transparent=True)
    plt.savefig(f"{path}mu={mu}_{comment}.pdf", bbox_inches='tight', dpi = 1200, format='pdf',transparent=True)

    plt.show()

def plot_statistics_output_linear_FINAL(mu, theo, recon, error, comment, path):

    plt.clf()
    plt.style.use('perso2.txt')
    x = [i for i in range(len(theo))]
    plt.bar(x, height=theo, width=0.2, bottom=None, align='center', color = "gold", edgecolor = "black", label = "Poissonian fit")

    plt.errorbar(x, recon, yerr = error, fmt='ro', label = "Detected statistics")
    # plt.title(r"$\mu$" + f" = {mu} , Pulse lenght = {pulse_lenght} ns")
    plt.xlabel("Photon number")
    plt.ylabel("Probability")
    plt.legend()
    tp = max(recon) + max(error)
    plt.xlim(left = -0.3, right=8.3)
    plt.ylim(bottom=0, top=1.05*tp)
    ####################
    # plt.ylim(bottom=0, top=1)
    plt.yscale('linear')
    #####################
    plt.grid(False)
    # comment = "_cutted_data_and_summed" #start with _ (underscore)
    plt.savefig(f"{path}mu={mu}_{comment}.png", bbox_inches='tight', dpi = 1200, format='png',transparent=True)
    plt.savefig(f"{path}mu={mu}_{comment}.pdf", bbox_inches='tight', dpi = 1200, format='pdf',transparent=True)

    plt.show()



def plot_statistics_input_log_FINAL(mu, theo, recon, error, comment, path):

    plt.clf()
    plt.style.use('perso2.txt')
    x = [i for i in range(len(theo))]
    plt.bar(x, height=theo, width=0.2, bottom=None, align='center', color = "gold", edgecolor = "black", label = "Poissonian statistics")

    plt.errorbar(x, recon, yerr = error, fmt='ro', label = "Reconstructed statistics")
    # plt.title(r"$\mu$" + f" = {mu} , Pulse lenght = {pulse_lenght} ns")
    plt.xlabel("Photon number")
    plt.ylabel("Probability")
    plt.legend()
    tp = max(recon) + max(error)
    plt.xlim(left = -0.3, right=8.3)
    # plt.ylim(bottom=0, top=1.05*tp)
    ####################
    plt.ylim(bottom=1e-6, top=1)
    plt.yscale('log')
    #####################
    plt.grid(False)
    # comment = "_cutted_data_and_summed" #start with _ (underscore)
    plt.savefig(f"{path}mu={mu}_{comment}.png", bbox_inches='tight', dpi = 1200, format='png',transparent=True)
    plt.savefig(f"{path}mu={mu}_{comment}.pdf", bbox_inches='tight', dpi = 1200, format='pdf',transparent=True)

    plt.show()


def plot_statistics_input_linear_FINAL(mu, theo, recon, error, comment, path):

    plt.clf()
    plt.style.use('perso2.txt')
    x = [i for i in range(len(theo))]
    plt.bar(x, height=theo, width=0.2, bottom=None, align='center', color = "gold", edgecolor = "black", label = "Poissonian statistics")

    plt.errorbar(x, recon, yerr = error, fmt='ro', label = "Reconstructed statistics")
    # plt.title(r"$\mu$" + f" = {mu} , Pulse lenght = {pulse_lenght} ns")
    plt.xlabel("Photon number")
    plt.ylabel("Probability")
    plt.legend()
    tp = max(recon) + max(error)
    plt.xlim(left = -0.3, right=8.3)
    plt.ylim(bottom=0, top=1.05*tp)
    ####################
    # plt.ylim(bottom=0, top=1)
    plt.yscale('linear')
    #####################
    plt.grid(False)
    # comment = "_cutted_data_and_summed" #start with _ (underscore)
    plt.savefig(f"{path}mu={mu}_{comment}.png", bbox_inches='tight', dpi = 1200, format='png',transparent=True)
    plt.savefig(f"{path}mu={mu}_{comment}.pdf", bbox_inches='tight', dpi = 1200, format='pdf',transparent=True)

    plt.show()


def plot_statistics_output_FINAL_14pixel(mu, pulse_lenght, theo, recon, error, comment, path):

    plt.clf()
    plt.style.use('perso2.txt')
    fig, ax1 = plt.subplots(figsize=(5.667,7.570))

    x = [i for i in range(len(theo))]
    plt.bar(x, height=theo, width=1, bottom=None, align='center', color = "#80CEFF", edgecolor = "grey", label = "Poissonian fit")

    plt.errorbar(x, recon, yerr = error, fmt='ro', ms=11, label = "Detected statistics")
    plt.tick_params(axis='x', which='both', top=False, direction='out', color = 'black', labelsize = 28)
    plt.tick_params(axis='y', which='both', right=False, direction='out', color='black', labelsize = 28)
    plt.xticks(np.arange(0,9,1))
    plt.yticks(np.arange(0,1,0.1))

    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(0.5) #0.5
        ax1.spines[axis].set_color('black')
        ax1.tick_params(axis='x', colors='black', labelsize = 28)
        ax1.tick_params(axis='y', colors='black', labelsize = 28)



    #ax = plt.axes()
    #ax.xaxis.set_major_locator(plt.MultipleLocator(2))
    #ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    font1 = {'family': 'Arial',
        'color':  'black',
        'weight': 'normal',
        'size': 28,
        }
    
    font2 = {'family': 'Arial',
        'color':  'black',
        'weight': 'normal',
        'size': 28,
        }
    
    font3 = {'family': 'Arial',
        'color':  'black',
        'weight': 'bold',
        'size': 28,
        }




    
    
    plt.rcParams["axes.axisbelow"] = True
    #plt.title(r"$\mu$" + f" = {mu}")#, Pulse duration = {pulse_lenght} ns")
    
    
    #plt.xlabel("Photon number",labelpad=10, fontdict=font1)
    #plt.ylabel("Probability",labelpad=10,fontdict=font1)
    
    plt.figtext(0.6,0.7,r"$\mu$" + f" = {mu}",fontdict=font2)
    plt.figtext(-0.05,0.91,"d",fontdict=font3)
    #plt.legend(fontsize="20")
    
    tp = max(recon) + max(error)
    plt.xlim(left = -0.5, right=8.5)
    plt.ylim(bottom=0, top=1.05*tp)
    #plt.yticks(np.arange(0,1.05*tp,0.2))
    ####################
    #plt.ylim(bottom=1e-5, top=1)
    #plt.yscale('log')
    #####################
    plt.grid(False)
    # comment = "_cutted_data_and_summed" #start with _ (underscore)
    #plt.savefig(f"{path}mu={mu}_pulse={pulse_lenght}ns{comment}_FINAL.png", bbox_inches='tight', dpi = 1200, format='png',transparent=True)
    #plt.savefig(f"{path}mu={mu}_pulse={pulse_lenght}ns{comment}_FINAL.pdf", bbox_inches='tight', dpi = 1200, format='pdf',transparent=True)
    plt.savefig(f"mu={mu}_pulse={pulse_lenght}ns{comment}_FINAL_sansSerif342.png", bbox_inches='tight', dpi = 1200, format='png',transparent=False)
    #plt.savefig(f"mu={mu}_pulse={pulse_lenght}ns{comment}_FINAL.pdf", bbox_inches='tight', dpi = 1200, format='pdf',transparent=True)
    #plt.show()




def plot_statistics_input_linear_blue_FINAL(mu, theo, recon, error, comment, path):

    plt.clf()
    plt.style.use('perso2.txt')
    fig, ax1 = plt.subplots(figsize=(5.667,7.570))

    x = [i for i in range(len(theo))]
    plt.bar(x, height=theo, width=1, bottom=None, align='center', color = "#80CEFF", edgecolor = "grey", label = "Poissonian statistics")

    plt.errorbar(x, recon, yerr = error, fmt='ro', ms=11, label = "Reconstructed statistics")
    plt.tick_params(axis='x', which='both', top=False, direction='out', color = 'black', labelsize = 28)
    plt.tick_params(axis='y', which='both', right=False, direction='out', color='black', labelsize = 28)
    plt.xticks(np.arange(0,9,1))
    plt.yticks(np.arange(0,1,0.1))
    plt.legend()

    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(0.8) #0.5
        ax1.spines[axis].set_color('black')

    ax1.tick_params(axis='x', colors='black', labelsize=20)
    for tick in ax1.get_xticklabels():
        tick.set_fontname("Arial")

    ax1.tick_params(axis='y', colors='black', labelsize=20)
    for tick in ax1.get_yticklabels():
        tick.set_fontname("Arial")

    font1 = {'family': 'Arial',
        'color':  'black',
        'weight': 'normal',
        'size': 24,
        }

    
    plt.rcParams["axes.axisbelow"] = True


    plt.xlabel("Photon number",labelpad=10, fontdict=font1)
    plt.ylabel("Probability",labelpad=10,fontdict=font1)
    


    # plt.figtext(0.6,0.7,r"$\mu$" + f" = {mu}",fontdict=font2)
    # plt.figtext(-0.05,0.91,"d",fontdict=font3)
    # #plt.legend(fontsize="20")
    plt.legend(loc=1, prop={'family': 'Arial', 'size': 16}, frameon=False)


    tp = max(recon) + max(error)
    plt.xlim(left = -0.5, right=8.5)
    plt.ylim(bottom=0, top=1.25*tp)
    #plt.yticks(np.arange(0,1.05*tp,0.2))
    ####################
    #plt.ylim(bottom=1e-5, top=1)
    #plt.yscale('log')
    #####################
    plt.grid(False)
    # comment = "_cutted_data_and_summed" #start with _ (underscore)
    plt.savefig(f"{path}mu={mu}_{comment}_blue_arial.png", bbox_inches='tight', dpi = 1200, format='png',transparent=True)
    plt.savefig(f"{path}mu={mu}_{comment}_blue_arial_lowquality.png", bbox_inches='tight', dpi = 600, format='png',transparent=True)
    plt.savefig(f"{path}mu={mu}_{comment}_blue_arial.pdf", bbox_inches='tight', dpi = 1200, format='pdf',transparent=True)


if __name__ == "__main__":
    mu = 2.05
    Reconstructed_Input_statistics = [0.12754491, 0.26424608, 0.27323563, 0.18588623, 0.09451657, 0.03729952, 0.01326559, 0.00236815, 0.00163732]
    error = [0,0,0,0,0,0,0,0,0]
    comment = "linear"
    path = ""
    Theoretical_statistics = [0.12859372, 0.26375823, 0.2704969,  0.18493849, 0.09483171, 0.03890181, 0.01329857, 0.00389666, 0.00099905]

    plot_statistics_input_linear_blue_FINAL(mu, Theoretical_statistics, Reconstructed_Input_statistics, error, comment, path)
