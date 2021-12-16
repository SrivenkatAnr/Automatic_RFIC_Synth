#===========================================================================================================================
"""
Name                : Srivenkat A
Roll Number         : EE18B038
File Name           : hand_calc.py
File Description    : This file will fix intial condition of circuit parameters using hand calculations 
                      for gradient descent algorithm

Functions structure in this file:
    --> db_to_normal
    
"""
#===========================================================================================================================
#=========================================== IMPORT FILES ==================================================================
import numpy as np
import common_functions as cff # type:ignore

#===========================================================================================================================


"""
===========================================================================================================================
------------------------------------Defining the functions for simple calculations-----------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# Converting from dB to normal
# Outputs : normal values
def db_to_normal(val_db):
    return 10**(val_db/10)

#-----------------------------------------------------------------------------------------------
# Calculating the resistance from the inductance values
# Outputs : resistance
def calculate_resistance_inductor(ind,fo,Q):
    res=2*np.pi*fo*ind/Q
    return res

#-----------------------------------------------------------------------------------------------
# Calculating Ld
# Outputs : Ld
def calculate_Ld():
    Ld=10e-9;
    return Ld

#-----------------------------------------------------------------------------------------------
# Calculating Rl
# Outputs : Rl
def calculate_Rl(Rd):
    Rl=max(50, 20*Rd)
    return Rl

#-----------------------------------------------------------------------------------------------
# Calculating gm
# Outputs : gm
def calculate_gm(gain,Rl):
    return gain/Rl

#-----------------------------------------------------------------------------------------------
# Calculating Io
# Outputs : Io
def calculate_Io(Vdd,Rd,Rl):
    return Vdd/(Rl+2*Rd)

#-----------------------------------------------------------------------------------------------
# Calculating W
# Outputs : W
def calculate_W(gm,Lmin,Id,un,Cox):
    return (gm**2)*Lmin/(2*Id*un*Cox)

#-----------------------------------------------------------------------------------------------
# Calculating coupling cap
# Outputs : C
def calculate_Ccoup(fo,R):
    wo=2*np.pi*fo
    return 10/(wo*R)

"""
===========================================================================================================================
-------------------------------------------- Main Functions ---------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function to calculate the Initial Circuit Parameters  
# Outputs : NONE
def calculate_initial_parameters(cir,optimization_input_parameters):
    
    output_conditions=optimization_input_parameters['output_conditions']
    
    # Getting the output conditions
    fo=output_conditions['wo']/(2*np.pi)
    Rin=output_conditions['Rin']
    gain=db_to_normal(output_conditions['gain_db'])
    Lmin=cir.mos_parameters['Lmin']
    Cox=cir.mos_parameters['cox']
    un=cir.mos_parameters['un']
    Vdd=cir.mos_parameters['Vdd']

    # Calculating the circuit parameters
    circuit_parameters={}
    circuit_parameters['Ld']=calculate_Ld()
    circuit_parameters['Rd']=calculate_resistance_inductor(circuit_parameters['Ld'],fo,50)
    circuit_parameters['Rl']=calculate_Rl(circuit_parameters['Rd'])
    circuit_parameters['Io']=calculate_Io(Vdd,circuit_parameters['Rd'],circuit_parameters['Rl'])
    gm=calculate_gm(gain,circuit_parameters['Rl'])
    circuit_parameters['W']=calculate_W(gm,Lmin,circuit_parameters['Io'],un,Cox)
    circuit_parameters['Rbias']=5000
    circuit_parameters['Ccoup_in']=calculate_Ccoup(fo,circuit_parameters['Rbias'])
    circuit_parameters['Ccoup_out']=calculate_Ccoup(fo,circuit_parameters['Rl'])

    # Running the circuit
    cir.update_circuit(circuit_parameters)

#---------------------------------------------------------------------------------------------------------------------------
# Function to update the Initial Circuit Parameters after calculating the new value of vt
# Inputs  : circuit_parameters, mos_parameters, extracted_parameters, optimization_input_parameters
# Outputs : circuit_parameters, dc_outputs, mos_parameters, extracted_parameters
"""
def update_initial_parameters(cir,optimization_input_parameters):

    i=0
    while i<5 and cir.extracted_parameters['s11_db']>-15.0:

        # Printing the iteration number
        i+=1
        print('----- Iteration ',i,' -----')

        # Updating the values
        fo=optimization_input_parameters['output_conditions']['wo']/(2*np.pi)
        cir.circuit_parameters['Ls']=cir.circuit_parameters['Ls']*50/cir.extracted_parameters['Zin_R']
        cir.circuit_parameters['Lg']=cir.circuit_parameters['Lg']-cir.extracted_parameters['Zin_I']/(2*np.pi*fo)
        
        # Running the circuit
        cir.run_circuit()
"""


"""
===========================================================================================================================
-------------------------------------------- Output Functions -------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function to perform complete pre-optimization with Method-1   
# Outputs : NONE
def automatic_initial_parameters(cir,optimization_input_parameters,optimization_results):
    
    
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Automatic Operating Point Selection 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


    #======================================================== Step 1 =============================================================================================================
    print('\n\n--------------------------------- Operating Point Calculations ------------------------------------')

    # Calculating the Values of Circuit Parameters
    calculate_initial_parameters(cir,optimization_input_parameters)

    # Storing the Circuit and Extracted Parameters
    optimization_results['auto_hc']={}
    optimization_results['auto_hc']['circuit_parameters']=cir.circuit_parameters.copy()
    optimization_results['auto_hc']['extracted_parameters']=cir.extracted_parameters.copy()

    # Printing the values
    cff.print_circuit_parameters(cir.circuit_parameters)
    cff.print_extracted_outputs(cir.extracted_parameters)
    """
    #======================================================== Step 2 =======================================================
    print('\n\n--------------------------------- Operating Point Updations ------------------------------------')

    # Calculating the Values of Circuit Parameters
    update_initial_parameters(cir,optimization_input_parameters)

    # Storing the Circuit and Extracted Parameters
    optimization_results['hc_update']={}
    optimization_results['hc_update']['circuit_parameters']=cir.circuit_parameters.copy()
    optimization_results['hc_update']['extracted_parameters']=cir.extracted_parameters.copy()

    # Printing the values
    cff.print_circuit_parameters(cir.circuit_parameters)
    cff.print_extracted_outputs(cir.extracted_parameters)
    """