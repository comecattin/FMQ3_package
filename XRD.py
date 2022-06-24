# -*- coding: utf-8 -*-
"""
Created on Mon May  9 16:47:20 2022

@author: C.Cattin
"""

import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
font = {'family' : 'serif',
        'serif' : 'cm10',
        'weight' : 'bold',
        'size'   : 16}
matplotlib.rc('font', **font)
matplotlib.rc('figure',figsize=(15,7))
import matplotlib.pyplot as plt



### FILE OPENNING ###
def open_file(path):
    '''
    Extract data from .txt file from XRD experiments
    Parameters
    ----------
    path : str
        Path to the file.

    Returns
    -------
    data
    '''
    
    with open(path) as file :
        data = file.read()
        file.close()

    #Extract data
    data = data.replace(' ','').replace('\n','').split(',')
    data = data[2:-1]
    data = np.array(data,dtype=np.float32)

    angle = data[0::2]
    PSD = data[1::2]


    return data, angle, PSD


def open_ref(path):
    '''
    Extract data from .txt reference XRD patern

    Parameters
    ----------
    path : str
        Path to the file

    Returns
    -------
    data : array
        All data extracted
    angle : np.array
        Diffraction angle
    intensity : np.array
        Diffraction intensity

    '''

    data = np.loadtxt(path,skiprows=1)
    angle = data[:,0]
    intensity = data[:,1]

    return data , angle , intensity


### DATA OPERATIONS ###

def normalize(intensity):
    '''
    Normalize the intensity of a XRD patern

    Parameters
    ----------
    intensity : np.array
        intensity

    Returns
    -------
    intensity : np.array
        Normalized intensity

    '''
    return (intensity-min(intensity))/(max(intensity)-min(intensity))


### PLOT ###

def plot_XRD(angle,PSD,file_name,
             ref=False, path_ref='',
             n_powder=1,
             color=['red'],label=['Powder 1'],label_ref=r'Reference',
             over = True):
    """
    Plot the XRD powder pattern

    Parameters
    ----------
    angle : np.array
        Angle of diffraction
    PSD : np.array
        Intensity
    file_name : str
        Name to save the file
    ref : Bool , optionnal
        put a reference pattern or not , Default is False
    path_ref : str
        path to the reference patern 
    n_powder : int
        Number of powder to plot
    color : list of str , optional
        Color of each line
    label : list of str , optional
        Label of each line
    label_ref : list of str , optional
        Label of the reference
    over : bool , optional
        Plot the XRD pattern on the same plot or not
        

    Returns
    -------
    None.

    """
    #Just one powder
    if n_powder == 1 :
        #Plot also the reference
        if ref :
            #Plots share the same x-axis
            fig,ax = plt.subplots(len(path_ref) +1 ,1,sharex=True)
            fig.subplots_adjust(hspace=0)

            #Experimental plot
            ax[0].plot(angle,PSD,
                       color=color,linewidth=1,label=label)
            ax[0].grid()
            ax[0].legend()

            for i,ref_i in enumerate(path_ref) :
                #Extract ref data
                data_ref_i , angle_ref_i, intensity_ref_i = open_ref(ref_i)
                #Normalize intensity
                intensity_ref_i = normalize(intensity_ref_i)
            
            
                #Reference plot
                ax[i+1].plot(angle_ref_i,intensity_ref_i,
                           color='black',linewidth=1,label=label_ref[i])

                ax[i+1].set_xlim(angle[0],angle[-1])
                ax[i+1].set_xlabel(r'$2\theta (deg.)$')
                ax[i+1].grid()
                ax[i+1].legend()

        #Don't plot the reference
        else :
            fig,ax = plt.subplots()

            ax.plot(angle,PSD,color='red',linewidth=1)
            ax.set_xlim(angle[0],angle[-1])
            ax.set_xlabel(r'$2\theta (deg.)$')
            ax.set_ylabel(r'PSD')
            ax.grid()


    # More powder
    else:

        if ref :

            if over :
                #Plot the patterns over each other
                #Plots share the same x-axis
                fig,ax = plt.subplots(len(path_ref)+1,1,sharex=True)
                fig.subplots_adjust(hspace=0)

                #Experimental plot
                for i in range(n_powder):
                    ax[0].plot(angle[i],
                               PSD[i],
                               color=color[i],
                               label=label[i],
                               linewidth=1)
                ax[0].grid()
                ax[0].set_yticks([0.25,0.5,0.75])
                ax[0].legend()

                for i,ref_i in enumerate(path_ref) :
                    #Extract ref data
                    data_ref_i , angle_ref_i, intensity_ref_i = open_ref(ref_i)
                    #Normalize intensity
                    intensity_ref_i = normalize(intensity_ref_i)
                    #Reference plot
                    ax[i+1].plot(angle_ref_i,intensity_ref_i,
                               color='black',linewidth=1,label=label_ref[i])

                    ax[i+1].set_xlim(np.min(angle),np.max(angle))
                    ax[i+1].set_xlabel(r'$2\theta (deg.)$')
                    ax[i+1].grid()
                    ax[i+1].set_yticks([0.25,0.5,0.75])
                    ax[i+1].legend()


            else :
                #Plot the patterns in different subplots
                #Plots share the same x-axis
                fig,ax = plt.subplots(n_powder+len(path_ref),1,sharex=True)
                fig.subplots_adjust(hspace=0)

                #Experimental plot
                for i in range(n_powder):
                    ax[i].plot(angle[i],
                               PSD[i],
                               color=color[i],
                               label=label[i],
                               linewidth=1)
                    ax[i].grid()
                    ax[i].legend()
                    ax[i].set_yticks([0.25,0.5,0.75])

                #Reference plot
                for i,ref_i in enumerate(path_ref) :
                    #Extract ref data
                    data_ref_i , angle_ref_i, intensity_ref_i = open_ref(ref_i)
                    #Normalize intensity
                    intensity_ref_i = normalize(intensity_ref_i)
                    #Reference plot
                    ax[i+n_powder].plot(angle_ref_i,intensity_ref_i,
                               color='black',linewidth=1,label=label_ref[i])

                    ax[i+n_powder].set_xlim(np.min(angle),np.max(angle))
                    ax[i+n_powder].set_xlabel(r'$2\theta (deg.)$')
                    ax[i+n_powder].grid()
                    ax[i+n_powder].set_yticks([0.25,0.5,0.75])
                    ax[i+n_powder].legend()









        #Don't plot the reference
        else:
            if over :
                #Plot in the same plot
                fig,ax = plt.subplots()

                for i in range(n_powder):
                    ax.plot(angle[i],PSD[i],color=color[i],label=label[i],linewidth=1)
    
                ax.set_xlabel(r'$2\theta (deg.)$')
                ax.set_ylabel(r'PSD')
                ax.set_xlim(angle[0],angle[-1])
                ax.grid()

            else :
                #Plot the patterns on different subplots
                #Plots share the same x-axis
                fig,ax = plt.subplots(n_powder,1,sharex=True)
                fig.subplots_adjust(hspace=0)

                for i in range(n_powder):
                    ax[i].plot(angle[i],PSD[i],color=color[i],label=label[i],linewidth=1)
                    ax[i].grid()
                    ax[i].legend()
                    ax[i].set_yticks([0.25,0.5,0.75])
                    ax[i].set_xlim(angle[i][0],angle[i][-1])
                ax[-1].set_xlabel(r'$2\theta (deg.)$')
                


    plt.legend()

    plt.savefig(file_name,bbox_inches='tight')

    plt.show()


### MAIN FUNCTION ###

def XRD(path,file_name,
        ref=False, path_ref='',
        n_powder=1,color='red',label='Sample 1',label_ref=r'Reference',
        over=True):
    '''
    

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name to save the file
    ref : Bool , optionnal
        put a reference pattern or not , Default is False
    path_ref : str
        path to the reference patern 
    n_powder : int
        Number of powder to plot
    color : list of str , optional
        Color of each line
    label : list of str , optional
        Label of each line
    label_ref : list of str , optional
        Label of the reference
    over : bool , optional
         Plot the XRD pattern on the same plot or not

    Returns
    -------
    angle : np.array
        Diffraction angles
    PSD : np.array
        Diffraction intensities

    '''

    if n_powder == 1 :
        #If only one powder
        data, angle, PSD = open_file(path)
        #Normalize intensity
        PSD = normalize(PSD)
        plot_XRD(angle, PSD,
                 file_name,
                 ref,path_ref,
                 n_powder,
                 color,label,label_ref)

    else :
        #If multiple powders
        #Initializing
        data , angle , PSD = open_file(path[0])
        #Normalize intensity
        PSD = normalize(PSD)

        #Loop over the powders
        for i in range(1,n_powder):
            #Extract data
            data_i , angle_i , PSD_i = open_file(path[i])
            #Normalize intensity
            PSD_i = normalize(PSD_i)
            #Append
            angle = np.vstack((angle,angle_i))
            PSD = np.vstack((PSD,PSD_i))


        #Plot
        plot_XRD(angle = angle,
                 PSD = PSD,
                 file_name = file_name,
                 ref = ref,
                 path_ref = path_ref,
                 n_powder = n_powder,
                 color = color,
                 label = label,
                 label_ref = label_ref,
                 over = over)

    return angle,PSD



if __name__=='__main__':

    #Path
    path = {
            #Li3Ni2SbO6
        'Li3Ni2SbO6_1' : '../XRD/Li3Ni2SbO6/Sample_1/Li3Ni2SbO6_1.txt',
        'Li3Ni2SbO6_2' : '../XRD/Li3Ni2SbO6/Sample_1/Li3Ni2SbO6_2.txt',
        'Li3Ni2SbO6_sample_2_1' : '../XRD/Li3Ni2SbO6/Sample_2/Li3Ni2SbO6_sample_2_1.txt',
        'Li3Ni2SbO6_sample_2_2' : '../XRD/Li3Ni2SbO6/Sample_2/Li3Ni2SbO6_sample_2_2.txt',
        'Li3Ni2SbO6_sample_3_1' : '../XRD/Li3Ni2SbO6/Sample_3/Li3Ni2SbO6_sample_3_1.txt',
        'Li3Ni2SbO6_sample_3_2' : '../XRD/Li3Ni2SbO6/Sample_3/Li3Ni2SbO6_sample_3_2.txt',
                #Protonated
        'Li3Ni2SbO6_proton_HT_1' : '../XRD/Li3Ni2SbO6/Proton_HT/Li3Ni2SbO6_proton_1.txt',
        'Li3Ni2SbO6_proton_LT_1' : '../XRD/Li3Ni2SbO6/Proton_LT/Li3Ni2SbO6_proton_1.txt',
        'Li3Ni2SbO6_proton_50_4M_48h' : '../XRD/Li3Ni2SbO6/Proton_HT/Li3Ni2SbO6_proton_50_48h_4M.txt',
        'Li3Ni2SbO6_proton_85_4M_48h' : '../XRD/Li3Ni2SbO6/Proton_HT/Li3Ni2SbO6_proton_85_48h_4M.txt',


            #Li2MnO3
        'Li2MnO3_1' : '../XRD/Li2MnO3/Sample_1/Li2MnO3_1.txt',
        'Li2MnO3_2' : '../XRD/Li2MnO3/Sample_1/Li2MnO3_2.txt',
        'Li2MnO3_sample_2_1' : '../XRD/Li2MnO3/Sample_2/Li2MnO3_sample_2_1.txt',
        'Li2MnO3_sample_2_2' : '../XRD/Li2MnO3/Sample_2/Li2MnO3_sample_2_2.txt',
        'Li2MnO3_sample_3_1' : '../XRD/Li2MnO3/Sample_3/Li2MnO3_sample_3_1.txt',
        'Li2MnO3_sample_3_2' : '../XRD/Li2MnO3/Sample_3/Li2MnO3_sample_3_2.txt',
        'Li2MnO3_sample_4_1' : '../XRD/Li2MnO3/Sample_4/Li2MnO3_sample_4_1.txt',
        'Li2MnO3_sample_4_2' : '../XRD/Li2MnO3/Sample_4/Li2MnO3_sample_4_2.txt',

                #Protonated
        'Li2MnO3_proton_HT_1' : '../XRD/Li2MnO3/Proton_HT/Li2MnO3_proton_1.txt',
        'Li2MnO3_proton_LT_1' : '../XRD/Li2MnO3/Proton_LT/Li2MnO3_proton_1.txt',
        'Li2MnO3_proton_50_1M_48h' : '../XRD/Li2MnO3/Proton_HT/Li2MnO3_proton_50_1M_48h.txt',
        'Li2MnO3_proton_50_2M_48h' : '../XRD/Li2MnO3/Proton_HT/Li2MnO3_proton_50_2M_48h.txt',
        'Li2MnO3_proton_50_4M_48h' : '../XRD/Li2MnO3/Proton_HT/Li2MnO3_proton_50_4M_48h.txt',
        'Li2MnO3_proton_85_4M_48h' : '../XRD/Li2MnO3/Proton_HT/Li2MnO3_proton_85_4M_48h.txt',


            #Li3Co2SbO6
        'Li3Co2SbO6_1' : '../XRD/Li3Co2SbO6/Sample_1/Li3Co2SbO6_1.txt',
        'Li3Co2SbO6_2' : '../XRD/Li3Co2SbO6/Sample_1/Li3Co2SbO6_2.txt',
        'Li3Co2SbO6_sample_2_1' : '../XRD/Li3Co2SbO6/Sample_2/Li3Co2SbO6_sample_2_1.txt',
        'Li3Co2SbO6_sample_2_2' : '../XRD/Li3Co2SbO6/Sample_2/Li3Co2SbO6_sample_2_2.txt',

            #Na3Co2SbO6
        'Na3Co2SbO6_1' : '../XRD/Na3Co2SbO6/Sample_1/Na3Co2SbO6_1.txt',
        'Na3Co2SbO6_2' : '../XRD/Na3Co2SbO6/Sample_1/Na3Co2SbO6_2.txt',
                #Mistake
        'Na3Co2SbO6_mistake_sample_2' : '../XRD/Na3Co2SbO6/Sample_2/Na3Co2SbO6_sample_2_1.txt',
                #Lower temperature
        'Na3Co2SbO6_sample_3_1' : '../XRD/Na3Co2SbO6/Sample_3/Na3Co2SbO6_sample_3_1.txt',
        'Na3Co2SbO6_sample_3_2' : '../XRD/Na3Co2SbO6/Sample_3/Na3Co2SbO6_sample_3_2.txt',

            #Na3Ni2SbO6
        'Na3Ni2SbO6_1' : '../XRD/Na3Ni2SbO6/Sample_1/Na3Ni2SbO6_1.txt',
        'Na3Ni2SbO6_2' : '../XRD/Na3Ni2SbO6/Sample_1/Na3Ni2SbO6_2.txt',

            #Co2P2O7
        'Co2P2O7_1' : '../XRD/Co2P2O7/Sample_1/Co2P2O7_1.txt',
        'Co2P2O7_2' : '../XRD/Co2P2O7/Sample_1/Co2P2O7_2.txt',
        'Co2P2O7_sample_2_450' : '../XRD/Co2P2O7/Sample_2_450C/Co2P2O7_sample_2_450.txt',
        'Co2P2O7_sample_2_600' : '../XRD/Co2P2O7/Sample_2_600C/Co2P2O7_sample_2_600.txt',
        'Co2P2O7_sample_3_900' : '../XRD/Co2P2O7/Sample_3_900C/Co2P2O7_sample_3_900.txt',
        'Co2P2O7_sample_3_1000' : '../XRD/Co2P2O7/Sample_3_1000C/Co2P2O7_sample_3_1000.txt',
        'Co2P2O7_beta_sample_1' : '../XRD/Co2P2O7/Beta/Sample_1/Co2P2O7_beta_sample_1.txt',


            #Na3Cu2SbO6
        'Na3Cu2SbO6_1' : '../XRD/Na3Cu2SbO6/Sample_1/Na3Cu2SbO6_1.txt',
        'Na3Cu2SbO6_sample_2_1' : '../XRD/Na3Cu2SbO6/Sample_2/Na3Cu2SbO6_sample_2_1.txt',
        'Na3Cu2SbO6_Ag_exchange' : '../XRD/Na3Cu2SbO6/Ag_exchange/Na3Cu2SbO6_Ag_exchange.txt',


            #Li3Cu2SbO6
        'Li3Cu2SbO6_1' : '../XRD/Li3Cu2SbO6/Sample_1/Li3Cu2SbO6_1.txt',
        'Li3Cu2SbO6_Ag_exchange' : '../XRD/Li3Cu2SbO6/Ag_exchange/Li3Cu2SbO6_Ag_exchange.txt',
            }



    #Reference
    reference = {
        'Li3Ni2SbO6' : '../References/Li3Ni2SbO6_theo.xy',
        'Li2MnO3' : '../References/Li2MnO3_theo.xy',
        'Li3Co2SbO6' : '../References/Li3Co2SbO6_theo.xy',
        'Na3Co2SbO6' : '../References/Na3Co2SbO6_theo.xy',
        'Na3Co2SbO6_real' : '../References/Na3Co2SbO6_xy.xy',
        'Na3Ni2SbO6' : '../References/Na3Ni2SbO6_theo.xy',
        'Co2P2O7' : '../References/Co2P2O7_theo.xy',
        'Co2P2O7_beta' : '../References/Co2P2O7_beta_theo.xy',
        'Co2P2O7_gamma' : '../References/Co2P2O7_gamma_theo.xy',
        'Na3Cu2SbO6' : '../References/Na3Cu2SbO6_theo.xy',
        'Li3Cu2SbO6' : '../References/Li3Cu2SbO6_theo.xy',
        'Ag3Cu2SbO6_from_Cu5SbO6' : '../References/Ag3Cu2SbO6_from_Cu5SbO6_simulation.txt',
        'Ag3Cu2SbO6_from_Ag3LiRu2O6' : '../References/Ag3Cu2SbO6_from_Ag3LiRu2O6_simulation.txt',
        }


    #File name of the plot
    file_name = '(Li-Na)3Cu2SbO6_Ag_exchange.pdf'

    #colors and labels of the plot
    color = ['blue','green','orange','red','violet']

    label = [r'$\mathrm{Li_3Cu_2SbO_6}$',
             r'$\mathrm{(Ag_3)Li_3Cu_2SbO_6}$',
             r'$\mathrm{Na_3Cu_2SbO_6}$',
             r'$\mathrm{(Ag_3)Na_3Cu_2SbO_6}$',
             r'$\mathrm{Co_2P_2O_7}$' + ' ' + r'$600^{\circ}C$',
             r'$\mathrm{Co_2P_2O_7}$' + ' ' + r'$900^{\circ}C$',
             r'$\mathrm{Co_2P_2O_7}$' + ' ' + r'$1000^{\circ}C$',
             r'$\mathrm{Na_3Co_2SbO_6}$'+' ' + r'$1000^{\circ}C$',
             r'$\mathrm{Co2P2O7}$'+' ' + r'$600^{\circ}C$'
             ]
    label_ref = [r'$\mathrm{Ag_3Cu_2SbO_6}$' + r' from ' + r'$\mathrm{Cu_5SbO_6}$',
                 r'$\mathrm{Ag_3Cu_2SbO_6}$' + r' from ' + r'$\mathrm{Ag_3LiRu_2O_6}$',]

    angle,PSD = XRD([path['Li3Cu2SbO6_1'],
                     path['Li3Cu2SbO6_Ag_exchange'],
                     path['Na3Cu2SbO6_1'],
                     path['Na3Cu2SbO6_Ag_exchange'],
                     ],
                    file_name,
                    n_powder=4,
                    ref=True,
                    path_ref=[reference['Ag3Cu2SbO6_from_Cu5SbO6'],
                              reference['Ag3Cu2SbO6_from_Ag3LiRu2O6']],
                    color=color,
                    label=label,
                    label_ref=label_ref,
                    over = False)





