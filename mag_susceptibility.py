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
import matplotlib.pyplot as plt


### OPEN FILE ###

def open_file(path):
    '''
    Extract data from .dat file from MPMS experiments
    Parameters
    ----------
    path : str
        Path to the file.

    Returns
    -------
    M magnetization, array
    H magnetic field, array
    '''
    
    with open(path) as file :
        data = file.read()
        file.close()
        
    #Splitting each line
    data = data.replace('\n',',').split(',')
    #Remove the header
    data = data[138:]
    #Get the title of each collumn
    title = data[1:69]
    #Remove the title of the collumnns
    data = data[69:]
    
    #Extract data
    comment = data[::len(title)]
    Time_stamp = np.array(data[1::len(title)],dtype=np.float32)
    Temperature = np.array(data[2::len(title)],dtype=np.float32)
    Magnetic_Field = np.array(data[3::len(title)],dtype=np.float32)
    center_position = np.array(data[10::len(title)],dtype=np.float32)
    range_data = np.array(data[13::len(title)],dtype=np.float32)
    min_temp = np.array(data[15::len(title)],dtype=np.float32)
    max_temp = np.array(data[16::len(title)],dtype=np.float32)
    min_field = np.array(data[17::len(title)],dtype=np.float32)
    max_field = np.array(data[18::len(title)],dtype=np.float32)
    DC_Moment_Free_Ctr = np.array(data[39::len(title)],dtype=np.float32)
    DC_Moment_Err_Free_Ctr = np.array(data[40::len(title)],dtype=np.float32)
    DC_Moment_Fixed_Ctr = np.array(data[37::len(title)],dtype=np.float32)
    DC_Moment_Err_Fixed_Ctr= np.array(data[38::len(title)],dtype=np.float32)

    return (data,
            title,
            comment,
            Time_stamp,
            Temperature,
            Magnetic_Field,
            center_position,
            range_data,
            min_temp,max_temp,
            min_field,max_field,
            DC_Moment_Free_Ctr,DC_Moment_Err_Free_Ctr,
            DC_Moment_Fixed_Ctr,DC_Moment_Err_Fixed_Ctr
            )


### MATH ###


def compute_chi(DC_Moment_Free_Ctr,Magnetic_Field,mass,Molar_Weight):
    """
    Compute magnetic susceptibility from magnetic moment, magnetic field,
    mass of the sample and the molar weight of the species.

    Parameters
    ----------
    DC_Moment_Free_Ctr : np.array
        magnetic moment
    Magnetic_Field : np.array
        applied magnetic field
    mass : float
        Mass of the sample
    Molar_Weight : float
        Molar weight of the sample

    Returns
    -------
    chi : np.array
        magnetic suceptibility
    inv_chi : np.array
        Inverse of the magnetic suceptibility

    """
    #Compute chi
    chi = DC_Moment_Free_Ctr/mass * Molar_Weight / Magnetic_Field
    #Compute the inverse of chi
    inv_chi = 1/chi
    return chi,inv_chi



def get_C_factor(inv_chi,Temperature):
    """
    Fit 1/chi versus the temperature with the Curie-Weiss law :

    \chi^{-1} = T/C - \theta_{CW}/C


    Parameters
    ----------
    inv_chi : np.array
        Inverse of the magnetic susceptibility
    Temperature : np.array
        Temperature of the system

    Returns
    -------
    C : float
        C factor in the Curie-Weiss law
    Theta_CW : float
        Theta in the Curie-Weiss law

    """

    #Only one magnetic field
    if Temperature.size == 1:
        #Only applied above 250K
        index = np.where(Temperature<=250)
        #Remove other data
        Temperature = np.delete(Temperature,index)
        inv_chi = np.delete(inv_chi,index)

        #Fit data
        inv_C,minus_Theta_over_C = np.polyfit(Temperature,inv_chi,1)

    #Various magnetic fields
    else :
        #Only applied above 250K
        #Initializing
        index = np.where(Temperature[0]<=250)
        #Remove other data
        Temperature_new = np.delete(Temperature[0],index)
        inv_chi_new = np.delete(inv_chi[0],index)

        #Fit data
        inv_C,minus_Theta_over_C = np.polyfit(Temperature_new,inv_chi_new,1)

        #Loop over all the different magnetic fields
        for i in range(1,Temperature.size):
            #ONly apply above 250K
            index = np.where(Temperature[i]<=250)
            #Remove other data
            Temperature_new_i = np.delete(Temperature[i],index)
            inv_chi_new_i = np.delete(inv_chi[i],index)

            #Fit data
            inv_C_i,minus_Theta_over_C_i = np.polyfit(Temperature_new_i,inv_chi_new_i,1)

            #Concatenate
            inv_C = np.vstack((inv_C,inv_C_i))
            minus_Theta_over_C = np.vstack((minus_Theta_over_C,minus_Theta_over_C_i))

    #Return the correct values
    C = 1/inv_C
    Theta_CW = - minus_Theta_over_C * C

    return C , Theta_CW



def compare_mu_eff(spin,C,g=2):
    """
    Compare the factor before mu_B for the mu_eff

    Parameters
    ----------
    spin : float
        spin of the system.
    C : float
        C factor in the Curie-Weiss law
    g : float, optional
        g factor. The default is 2.

    Returns
    -------
    theory : float
        Factor in theory.
    experiment : float
        Factor in the experiment.

    """

    theory = g*np.sqrt(spin*(spin+1))
    experiment = np.sqrt(8*C)

    print('Theory : {}'.format(theory))
    print('Experiment : {}'.format(experiment))

    return theory,experiment


### PLOT ###

def plot_chi_temperature(chi,inv_chi,
                         Temperature,
                         mag_field,
                         file_save=['chi_vs_T.pdf',
                                    'inv_chi_vs_T.pdf'],
                         color=['red','blue'],
                         label = r'B=0.1T'):
    """
    Generate a plot of chi and 1/chi versus the temperature

    Parameters
    ----------
    chi : np.array
        magnetic suceptibility
    inv_chi : np.array
        Inverse of the magnetic suceptibility
    Temperature : np.array
        Temperature of the system.
    mag_field : array
        magnetic fields applied to the system
    file_save : array of str, optionnal
        file name to save the plot, the default is ['chi_vs_T.pdf','inv_chi_vs_T.pdf']
    color : array of str , optionnal
        color of the different plot, the default is ['red','blue']
    label : array of str , optionnal
        label of each data set (Chemical compound and magnetic field for example)

    Returns
    -------
    fig and ax of the different plots

    """

    #Fit C and Theta_CW via the Curie-Weiss law
    C, Theta_CW = get_C_factor(inv_chi, Temperature)

    #Plot chi versus temperature
    fig_1,ax_1 = plt.subplots(1)

    #If only one magnetic field
    if chi.size == 1 :
        ax_1.plot(Temperature,chi,'o',color='red',markersize=3,mfc='none',
                  label=label)

    #If various magnetic fields
    else :
        #Loop over all the different magnetic fields
        for i in range(chi.size):
            ax_1.plot(Temperature[i],chi[i],
                    'o',color=color[i],markersize=3,mfc='none',
                    label=label[i])

    #General properties
    ax_1.grid()
    ax_1.set_xlabel(r'$T(K)$')
    ax_1.set_ylabel(r'$\chi(emu.mol^{-1})$')
    ax_1.legend()

    #Save the figure
    plt.savefig(file_save[0],bbox_inches='tight')



    #Plot 1/chi versus the temperature and the fitting curve
    #1/chi vs T
    fig_2,ax_2 = plt.subplots(1)
    #If only one magnetic field
    if chi.size == 1:
        ax_2.plot(Temperature,inv_chi,'o',color='red',markersize=3,mfc='none',
                  label=label)

        #Fitted data
        #Only applied above 250K
        index = np.where(Temperature<=250)
        #Remove other data
        Temperature = np.delete(Temperature,index)
        inv_chi = np.delete(inv_chi,index)
        #Plot
        ax_2.plot(Temperature,Temperature/C - Theta_CW/C,color='black',label='Fitted curve')

    #If various magnetic fields
    else :
        #Loop over all the different magnetic fields
        for i in range(chi.size):
            ax_2.plot(Temperature[i],inv_chi[i],
                    'o',color=color[i],markersize=3,mfc='none',
                    label=label[i])

            #Fitted data
            #Only applied above 250K
            index = np.where(Temperature[i]<=250)
            #Remove other data
            Temperature_new = np.delete(Temperature[i],index)
            #Plot
            ax_2.plot(Temperature_new,Temperature_new/C[i] - Theta_CW[i]/C[i],
                    color='black',label='Fitted curve')

    #General properties
    ax_2.grid()
    ax_2.set_xlabel(r'$T(K)$')
    ax_2.set_ylabel(r'$\chi^{-1}(emu^{-1}.mol)$')
    ax_2.legend()
    plt.savefig(file_save[1],bbox_inches='tight')



    plt.show()
    return fig_1,ax_1,fig_2,ax_2

### MAIN FUNCTION ###

def main_susceptibility(path,
                        mass,Molar_Weight,spin,mag_field,
                        file_save,
                        abnormal_value=False,
                        color = ['red','blue'],
                        label = r'B=0.1T'):
    '''
    Main function of the programm

    Parameters
    ----------
    path : str
        Path to the files.
    mass : float
        Mass of the system in g
    Molar_Weight : float
        Molar weight of the system in g/mol
    spin : float
        Spin of the system
    mag_field : array
        Applied magnetic fields in T
    file_save : list of str
        path to save the plot : [chi_VS_T.pdf,inv_chi_VS_T]
    abnormal_value : list, optional
        List of the index of abnormal values. The default is False.
    color : array of str , optionnal
        color of the different plot, the default is ['red','blue']
    label : array of str , optionnal
        label of each data set (Chemical compound and magnetic field for example)

    Returns
    -------
    C : float
        C factor in the Curie-Weiss law
    Theta_CW : float
        Theta in the Curie-Weiss law
    theory : float
        Factor in theory.
    experiment : float
        Factor in the experiment.

    '''

    #If only one magnetic field
    if type(path) == str :
        #Extract data
        (data,
         title,
         comment,
         Time_stamp,
         Temperature,
         Magnetic_Field,
         center_position,
         range_data,
         min_temp,max_temp,
         min_field,max_field,
         DC_Moment_Free_Ctr,DC_Moment_Err_Free_Ctr,
         DC_Moment_Fixed_Ctr,DC_Moment_Err_Fixed_Ctr
         ) = open_file(path)

        #Compute chi
        chi,inv_chi = compute_chi(DC_Moment_Free_Ctr, Magnetic_Field, mass, Molar_Weight)
    
        #Delete abnormal values
        if abnormal_value != False :
            chi = np.delete(chi,abnormal_value)
            inv_chi = np.delete(inv_chi,abnormal_value)
            Temperature = np.delete(Temperature,abnormal_value)

    #Various magnetic fields
    else:
        #Extract data
        #Initialization
        (data,
         title,
         comment,
         Time_stamp,
         Temperature,
         Magnetic_Field,
         center_position,
         range_data,
         min_temp,max_temp,
         min_field,max_field,
         DC_Moment_Free_Ctr,DC_Moment_Err_Free_Ctr,
         DC_Moment_Fixed_Ctr,DC_Moment_Err_Fixed_Ctr
         ) = open_file(path[0])

        #Compute chi
        chi,inv_chi = compute_chi(DC_Moment_Free_Ctr, Magnetic_Field, mass, Molar_Weight)

        #Delete abnormal values
        if abnormal_value != False :
            chi = np.delete(chi,abnormal_value)
            inv_chi = np.delete(inv_chi,abnormal_value)
            Temperature = np.delete(Temperature,abnormal_value)

        #Loop over all the different magnetic fields
        for i in range(1,len(path)):
            #Extract data
            (data,
             title,
             comment,
             Time_stamp,
             Temperature_i,
             Magnetic_Field,
             center_position,
             range_data,
             min_temp,max_temp,
             min_field,max_field,
             DC_Moment_Free_Ctr,DC_Moment_Err_Free_Ctr,
             DC_Moment_Fixed_Ctr,DC_Moment_Err_Fixed_Ctr
             ) = open_file(path[i])

            #Compute chi
            chi_i,inv_chi_i = compute_chi(DC_Moment_Free_Ctr, Magnetic_Field, mass, Molar_Weight)

            #Delete abnormal values
            if abnormal_value != False :
                chi_i = np.delete(chi_i,abnormal_value)
                inv_chi_i = np.delete(inv_chi_i,abnormal_value)
                Temperature_i = np.delete(Temperature_i,abnormal_value)

            #Concatenate
            #chi = np.vstack((chi,chi_i))
            #inv_chi = np.vstack((inv_chi,inv_chi_i))
            #Temperature = np.vstack((Temperature,Temperature_i))
            chi = [chi]
            chi.append(chi_i)
            chi = np.array(chi,dtype=object)

            inv_chi = [inv_chi]
            inv_chi.append(inv_chi_i)
            inv_chi = np.array(inv_chi,dtype=object)

            Temperature = [Temperature]
            Temperature.append(Temperature_i)
            Temperature = np.array(Temperature,dtype=object)





    #Do the plot
    fig_1 , ax_1 , fig_2 , ax_2 = plot_chi_temperature(chi,inv_chi,Temperature,
                         mag_field=mag_field,
                         file_save=file_save,color=color,label=label)

    #Get the C factor
    C , Theta_CW = get_C_factor(inv_chi, Temperature)

    #Compare results
    theory,experiment = compare_mu_eff(spin, C)

    return C , Theta_CW , theory , experiment,chi,inv_chi,Temperature, fig_1,ax_1,fig_2,ax_2


if __name__=='__main__':

    path = {
            #06/05 Experiment on Li2MnO3
            '06/05_01T' : '../Mag/Li2MnO3/06_05/06_05_Li2MnO3_01T.dat',
            '06/05_1T' : '../Mag/Li2MnO3/06_05/06_05_Li2MnO3_1T.dat',

            #11/05 Experiment on Li3Ni2SbO6
            '11/05_01T' : '../Mag/Li3Ni2SbO6/11_05/11_05_Li3Ni2SbO6_01T.dat',
            '11/05_1T' : '../Mag/Li3Ni2SbO6/11_05/11_05_Li3Ni2SbO6_1T.dat',

            #20/05 Experiment on Li3Ni2SbO6 Protonated HT
            '20/05_01T_size' : '../Mag/Li3Ni2SbO6/20_05_Proton_HT/20_05_Li3Ni2SbO6_proton_HT_01T_size.dat',
            '20/05_01T' : '../Mag/Li3Ni2SbO6/20_05_Proton_HT/20_05_Li3Ni2SbO6_proton_HT_01T.dat',
            '20/05_1T' : '../Mag/Li3Ni2SbO6/20_05_Proton_HT/20_05_Li3Ni2SbO6_proton_HT_1T.dat',

            #25/05 Experiment on Li2MnO3 protonated 50C 4M 48h
            '25/05_01T' : '../Mag/Li2MnO3/25_05/25_05_Li2MnO3_proton_50_4M_48h_01T.dat',
            '25/05_1T' : '../Mag/Li2MnO3/25_05/25_05_Li2MnO3_proton_50_4M_48h_1T.dat',
            '25/05_1T_size' : '../Mag/Li2MnO3/25_05/25_05_Li2MnO3_proton_50_4M_48h_1T_size.dat',

            #07/06 Experiment on Na3Ni2SbO6
            '07/06_01T' : '../Mag/Na3Ni2SbO6/07_06/07_06_Na3Ni2SbO6_01T.dat',
            '07/06_1T' : '../Mag/Na3Ni2SbO6/07_06/07_06_Na3Ni2SbO6_1T.dat',
            '07/06_01T_size' : '../Mag/Na3Ni2SbO6/07_06/07_06_Na3Ni2SbO6_01T_size.dat',

            #09/06 Experiment on Li3Co2SbO6
            '09/06_01T' : '../Mag/Li3Co2SbO6/09_06/09_06_Li3Co2SbO6_01T.dat',
            '09/06_1T' : '../Mag/Li3Co2SbO6/09_06/09_06_Li3Co2SbO6_1T.dat',

            #10/06 Experiment on Co2P2O7
            '10/06_01T' : '../Mag/Co2P2O7/10_06/10_06_Co2P2O7_01T.dat',
            '10/06_01T_corrected' : '../Mag/Co2P2O7/10_06/10_06_Co2P2O7_01T_corrected.dat',
            '10/06_1T' : '../Mag/Co2P2O7/10_06/10_06_Co2P2O7_1T.dat',
            '10/06_1T_corrected' : '../Mag/Co2P2O7/10_06/10_06_Co2P2O7_1T_corrected.dat',

            #14/06 Experiment on Na3Cu2SbO6
            '14/06_01T' : '../Mag/Na3Cu2SbO6/14_06/14_06_Na3Cu2SbO6_01T.dat',
            '14/06_1T' : '../Mag/Na3Cu2SbO6/14_06/14_06_Na3Cu2SbO6_1T.dat',

            #15/06 Experiment on Li3Cu2SbO6
            '15/06_01T' : '../Mag/Li3Cu2SbO6/15_06/15_06_Li3Cu2SbO6_01T.dat',
            '15/06_1T' : '../Mag/Li3Cu2SbO6/15_06/15_06_Li3Cu2SbO6_1T.dat',
            #17/06 Experiment on Li3Cu2SbO6 to correct 1T values
            '17/06_1T' : '../Mag/Li3Cu2SbO6/17_06/17_06_Li3Cu2SbO6_1T.dat',

            #23/06 Experiment on (Ag3)Na3Cu2SbO6
            '23/06_01T' : '../Mag/Na3Cu2SbO6/23_06/23_06_(Ag3)Na3Cu2SbO6_01T.dat',
            '23/06_1T' : '../Mag/Na3Cu2SbO6/23_06/23_06_(Ag3)Na3Cu2SbO6_1T.dat',




        }

    mass = {
            #06/05 Experiment on Li2MnO3
            '06/05' : 28.60e-3,
            #11/05 Experiment on Li3Ni2SbO6
            '11/05' : 20.85e-3,
            #20/05 Experiment on Li3Ni2SbO6 Protonated HT
            '20/05' : 20.57e-3,
            #25/05 Experiment on Li2MnO3 protonated 50C 4M 48h
            '25/05' : 21.17e-3,
            #07/06 Experiment on Na3Ni2SbO6
            '07/06' : 20.88e-3,
            #09/06 Experiment on Li3Co2SbO6
            '09/06' : 23.21e-3,
            #10/06 Experiment on Co2P2O7
            '10/06' : 20.44e-3,
            #14/06 Experiment on Na3Cu2SbO6
            '14/06' : 23.50e-3,
            #15/06 Experiment on Li3Cu2SbO6
            '15/06' : 20.55e-3,
            #23/06 Experiment on (Ag3)Na3Cu2SbO6
            '23/06' : 23.91e-3,
        }

    Molar_Weight = {
        'Li2MnO3' : 116.82,
        'Li3Ni2SbO6' : 355.9662,
        'Na3Ni2SbO6' : 404.11,
        'Li3Co2SbO6' : 356.45,
        'Co2P2O7' : 291.81,
        'Na3Cu2SbO6' : 413.82,
        'Li3Cu2SbO6' : 365.67,
        'Ag3Cu2SbO6' : 668.45,
        }

    abnormal_value = {
        '06/05' : [73,74],
        '20/05' : [142,143],
        '20/05_size' : [130,131],
        '25/05' : [2,363],
        '07/06' : [327,190,220,221,224,328],
        '10/06' : [129,324],
        '14/06_1T' : [55,139]

        }

    spin = {
        '06/05' : 3/2,
        '11/05' : 1,
        '09/06' : 1/2,
        '10/06' : 3/2,
        '14/06' : 1/2,
        '15/06' : 1/2
        }

    file_save = ['chi_vs_T_01T_1T.pdf','inv_chi_vs_T_01T_1T.pdf']


    C , Theta_CW , theory , experiment, chi , inv_chi , Temperature,fig_1,ax_1,fig_2,ax_2 = main_susceptibility(
        path = [path['23/06_01T'],path['23/06_1T']],
        mass = mass['23/06'],
        Molar_Weight = Molar_Weight['Ag3Cu2SbO6'],
        spin = 1/2,
        mag_field = [0.1,1],
        file_save = file_save,
        abnormal_value = False)

#Devide by 2 when there is 2 metallic atoms
    experiment = experiment / 2

