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

    if resistance<sheet_resistance:
	    length=L_min
	    width=L_min*sheet_resistance/resistance
    else:
	    width=W_min
	    length=W_min*resistance/sheet_resistance
	
	#width=W_min-dW
	#length=width*resistance/sheet_resistance
	
    return length,width

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
