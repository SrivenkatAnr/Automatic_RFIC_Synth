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
# Calculating coupling cap
# Outputs : C
def calculate_Ccoup(fo,R,thresh):
    wo=2*np.pi*fo
    return thresh/(wo*R)

#-----------------------------------------------------------------------------------------------
# Calculating drain inductance
# Outputs : Ld
def calculate_Ld():
    Ld_max=10e-9;
    return Ld_max

#-----------------------------------------------------------------------------------------------
# Calculating max op_swing
# Outputs : Vout_max
def calculate_op_swing(Vdd,gain):
    Vout_max=Vdd/(1+(2/gain))
    return Vout_max

#-----------------------------------------------------------------------------------------------
# Calculating load resistance
# Outputs : Rl
def calculate_Rl(Vout_max,Pout):
    Rl=(Vout_max**2)/(2*Pout)
    return Rl

#-----------------------------------------------------------------------------------------------
# Calculating bias resistance
# Outputs : Rb
def calculate_Rb(W,wo):
    cap_per_wid=1.2e-9
    gate_cap=W*cap_per_wid
    gate_imp=1/(wo*gate_cap)
    Rb=max(5000,50*gate_imp)
    return Rb

#-----------------------------------------------------------------------------------------------
# Calculating MOSFET width, bias current
# Outputs : W,Io
def calculate_W_Io(Vout_max,gain,Lmin,Rl,un,cox):
    gm=gain/Rl 
    Vdsat=Vout_max/gain
    W=(gm*Lmin)/(un*cox*Vdsat)
    I=0.5*un*cox*(W/Lmin)*(Vdsat**2)
    return (W,I)

#-----------------------------------------------------------------------------------------------
# Calculating Matching Network params
# Outputs : Lsrc,Lload,Cmn
def calculate_MN_params(Rl_src):
    return(Lsrc,Lload,Cmn)

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
    gain=db_to_normal(output_conditions['gain_db']/2)/2
    Pout=db_to_normal(output_conditions['op1db']-3)*1e-3
    Lmin=cir.mos_parameters['Lmin']
    Cox=cir.mos_parameters['cox']
    un=cir.mos_parameters['un']
    Vdd=cir.mos_parameters['Vdd']

    # Calculating the circuit parameters
    circuit_parameters={}
    circuit_parameters['Ld']=calculate_Ld()
    Vout_max=calculate_op_swing(Vdd,gain)
    #circuit_parameters['Rl']=2*calculate_Rl(Vout_max,Pout)
    Rl_single_exp=calculate_Rl(Vout_max,Pout)
    circuit_parameters['Lsrc'],circuit_parameters['Lload'],circuit_parameters['Cmn']=calculate_MN_params(Rl_single_exp)
    circuit_parameters['W'],circuit_parameters['Io']=calculate_W_Io(Vout_max,gain,Lmin,circuit_parameters['Rl'],un,Cox)
    circuit_parameters['Rb']=calculate_Rb(circuit_parameters['W'],output_conditions['wo'])
    circuit_parameters['Rin']=Rin
    #C1_thresh=cir.circuit_initialization_parameters['simulation']['standard_parameters']['C1_threshold']
    #C2_thresh=cir.circuit_initialization_parameters['simulation']['standard_parameters']['C2_threshold']

    #circuit_parameters['Rd']=calculate_resistance_inductor(circuit_parameters['Ld'],fo,50)
    #circuit_parameters['C1']=calculate_Ccoup(fo,circuit_parameters['Rb'],C1_thresh)
    #circuit_parameters['C2']=calculate_Ccoup(fo,circuit_parameters['Rl'],C2_thresh)

    # Running the circuit
    cir.update_circuit(circuit_parameters)

#---------------------------------------------------------------------------------------------------------------------------
# Function to update the Initial Circuit Parameters after calculating the new value of vt
# Inputs  : circuit_parameters, mos_parameters, extracted_parameters, optimization_input_parameters
# Outputs : circuit_parameters, dc_outputs, mos_parameters, extracted_parameters

def update_initial_parameters(cir,optimization_input_parameters):

    #setting Rin to 50
    i=0
    Rin_exp=cir.circuit_parameters['Rin']
    Rin_ext=0
    while i<3 and abs(Rin_ext-Rin_exp)/Rin_exp > 0.01:

        # Printing the iteration number
        i+=1
        print('----- Iteration ',i,' -----')

        # Updating the values
        Rin_ext = cir.extracted_parameters['Rin_ext']
        cir.circuit_parameters['Rin']=cir.circuit_parameters['Rin']*(Rin_exp/Rin_ext)
        # Running the circuit
        cir.run_circuit()

    #increasing W to match specs
    i=0
    gain_exp=db_to_normal(optimization_input_parameters['output_conditions']['gain_db']/2)
    Rl=cir.circuit_parameters['Rl']
    gm_exp=1.2*gain_exp/Rl
    while i<5 and cir.extracted_parameters['op1db_man']<optimization_input_parameters['output_conditions']['op1db']:

        # Printing the iteration number
        i+=1
        print('----- Iteration ',i,' -----')

        # Updating the values
        gm_ext = cir.extracted_parameters['gm']
        cir.circuit_parameters['W']=cir.circuit_parameters['W']*((gm_exp/gm_ext)**2)
        # Running the circuit
        cir.run_circuit()



"""
===========================================================================================================================
-------------------------------------------- Output Functions -------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function to perform complete pre-optimization with Method-1   
# Outputs : NONE
def automatic_initial_parameters(cir,optimization_input_parameters,optimization_results):
    
    """
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Automatic Operating Point Selection 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    """

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
    cff.print_extracted_parameters(cir.extracted_parameters)
    
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
    cff.print_extracted_parameters(cir.extracted_parameters)
   
