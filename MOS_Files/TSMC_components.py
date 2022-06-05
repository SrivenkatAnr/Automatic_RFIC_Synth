#===========================================================================================================================
"""
Name                : Srivenkat A
Roll Number         : EE18B038
File Name           : TSMC_components.py
File Description    : This file will extract parameters from TSMC model file and and define other passive components

Functions structure in this file:
    --> db_to_normal
    
"""
#===========================================================================================================================
#=========================================== IMPORT FILES ==================================================================
import numpy as np
import pandas as pd
import common_functions as cff # type:ignore

#===========================================================================================================================   




#===========================================================================================================================
#------------------------------------ MODEL FILE EXTRACTION --------------------------------------------------------------------

#-----------------------------------------------------------------
# Function that extracts the MOSFET File Parameeters
# Inputs  : Optimization Input Parameters
# Outputs : MOS_Parameters
def calculate_mos_parameters(circuit_initialization_parameters):
    #circuit_initialization_parameters=self.circuit_initialization_parameters
        
    # Setting Lmin and Vdd
    Lmin=circuit_initialization_parameters['MOS']['Lmin']
    vdd=circuit_initialization_parameters['MOS']['Vdd']
    cox=circuit_initialization_parameters['MOS']['cox']
    un=circuit_initialization_parameters['MOS']['un']
    vt=circuit_initialization_parameters['MOS']['vt']

    # Extracting From File
    mos_parameters = {'un':un,'cox':cox,'vt':vt,'Lmin':Lmin,'Vdd':vdd}
        
    # Printing the MOSFET Parameters
    cff.print_MOS_parameters(mos_parameters)

    return mos_parameters

#-----------------------------------------------------------------      
# Function that converts resistance to length and width
def get_TSMC_resistor(resistance):
    sheet_resistance=124.45
    L_min=0.8e-6
    W_min=0.4e-6
    dW=0.0691e-6

    if resistance<(sheet_resistance*L_min/W_min):
	    length=L_min
	    width=L_min*sheet_resistance/resistance
    else:
	    width=W_min
	    length=W_min*resistance/sheet_resistance
	
	#width=W_min-dW
	#length=width*resistance/sheet_resistance
	
    return length,width


"""
===========================================================================================================================
------------------------------------- TSMC CAPACITOR FUNCTIONS -------------------------------------------------------------
"""

#-----------------------------------------------------------------      
# Function that converts capacitance to length and width for MOS capacitor
def calculate_MOS_capacitor(cap):
    cox=17.25*1e-3
    w_check=np.sqrt(cap/cox)
    if w_check>2e-5:
        length=2e-5
        width=cap/(cox*length)
        return width,length
    if w_check<1.2e-7:
        width=1.2e-7
        length=cap/(cox*width)
        return width,length
    return w_check,w_check

#-----------------------------------------------------------------      
# Function that converts capacitance to length, width and mf for MiM capacitor
def calculate_MiM_capacitor(cap):
    if (cap>10e-12):
        mf = 1 + int(cap/10e-12)
        w_max=100e-6
        return w_max,w_max,mf
    else:
        mf = 1 + int(cap/5.63e-15)
        w_min=2e-6
        return w_min,w_min,mf

"""
===========================================================================================================================
------------------------------------- TSMC INDUCTOR FUNCTIONS -------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Finding the location of L
def find_target_L(L_vals,L):
    return np.abs(L_vals-L).argmin()

#---------------------------------------------------------------------------------------------------------------------------
# Finding the best point
def load_TSMC_Inductor_data():

	# Reading the data
	file_directory='/home/ee18b038/Auto_Ckt_Synth_Codes/TSMC_inductor_sweep/'
	#file_directory='C:/Users/roope/Studies/IIT/Prof Projects/Circuit_Synthesis/Extra_Codes/'
	data=pd.read_csv(file_directory+'inductor_sweep_1.csv')
	return data

#---------------------------------------------------------------------------------------------------------------------------
# Finding the best point
def find_TSMC_Inductor(Q,L,data):

	#Print Finding Q and L
	#print('Finding the following inductor: Q={0}, L={1}'.format(Q,L))
	L_vals=data['L']

	# Finding the location of 0.98L and 1.02L
	n1=find_target_L(L_vals,L*0.98)
	n2=find_target_L(L_vals,L*1.02)

	# Finding the loss of all the points from n1 and n2
	loss_iter=[]
	for i in range(n1,1+n2):
		Qi=data.at[i,'Q']
		Li=data.at[i,'L']
		loss=((Qi-Q)/Q)**2
		loss+=((Li-L)/L)**2
		loss_iter.append(loss)

	# Finding the location of the smallest value in loss iter
	min_value=min(loss_iter)
	min_index=loss_iter.index(min_value)
	min_index+=n1

	return min_index,data.at[min_index,'Width'],data.at[min_index,'Radius'],data.at[min_index,'N_Turns'],data.at[min_index,'gdis'],data.at[min_index,'spc']

