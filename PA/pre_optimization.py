#===========================================================================================================================
"""
Name                : Srivenkat A
Roll Number         : EE18B038
File Name           : pre_optimisation.py
File Description    : This file will perform pre optimization and calculate an initial point to be used 
                      at the start of the gradient descent algorithm

Functions structure in this file:
    --> pre_optimization
        --> save_input_results_pre_optimization
        --> save_mos_results
        --> manual_circuit_parameters
        --> calculate_mos_parameters
        --> save_output_results_pre_optimization

"""

#===========================================================================================================================
import datetime
import common_functions as cff     # type: ignore
import PA.hand_calc as hc1 # type: ignore
from collections import OrderedDict
"""
===========================================================================================================================
----------------------------------- File Writing Functions ----------------------------------------------------------------
"""


#-----------------------------------------------------------------
# Function that stores input data of the simulation
# Inputs  : optimization_input_parameters
# Outputs : NONE
def save_input_results_pre_optimization(optimization_input_parameters):

    # Opening the file
    filename=optimization_input_parameters['filename']['output']+str('/input_data.txt')
    f=open(filename,'a')

    # Storing the results
    f.write('\n\n********************************************************************************\n')
    f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pre Optimization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    f.write('\nPre_Opt_Type    :'+str(optimization_input_parameters['pre_optimization']['type']))
    
    f.close()

#-----------------------------------------------------------------
# Function that stores output data of the pre optimization
# Inputs  : optimization_results, optimization_input_parameters
# Outputs : NONE
def save_output_results_pre_optimization(optimization_results,optimization_input_parameters):
    
    # Opening the file
    filename=optimization_input_parameters['filename']['output']+str('/output_data.txt')
    f=open(filename,'a')

    # Storing the results
    f.write('\n\n********************************************************************************\n')
    f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pre Optimization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    if 'manual_hc' in optimization_results:
        f.write('\n\n--------------------- Manual Hand Calculations ---------------------------------')
        f.write('\n\n---------------- Circuit Parameters Complete ---------------')
        cff.print_output_parameters_complete(f,optimization_results['manual_hc']['circuit_parameters'])
        f.write('\n\n---------------- Circuit Parameters ------------------------')
        cff.print_output_parameters(f,optimization_results['manual_hc']['circuit_parameters'])
        f.write('\n\n---------------- Extracted Parameters ------------------------')
        cff.print_output_parameters(f,optimization_results['manual_hc']['extracted_parameters'])

    if 'auto_hc' in optimization_results:
        f.write('\n\n--------------------- Automatic Hand Calculations ---------------------------------')
        f.write('\n\n---------------- Circuit Parameters Complete ---------------')
        cff.print_output_parameters_complete(f,optimization_results['auto_hc']['circuit_parameters'])
        f.write('\n\n---------------- Circuit Parameters ------------------------')
        cff.print_output_parameters(f,optimization_results['auto_hc']['circuit_parameters'])
        f.write('\n\n---------------- Extracted Parameters ------------------------')
        cff.print_output_parameters(f,optimization_results['auto_hc']['extracted_parameters'])

    if 'hc_update' in optimization_results:
        f.write('\n\n--------------------- Hand Calculations Update ---------------------------------')
        f.write('\n\n---------------- Circuit Parameters Complete ---------------')
        cff.print_output_parameters_complete(f,optimization_results['hc_update']['circuit_parameters'])
        f.write('\n\n---------------- Circuit Parameters ------------------------')
        cff.print_output_parameters(f,optimization_results['hc_update']['circuit_parameters'])
        f.write('\n\n---------------- Extracted Parameters ------------------------')
        cff.print_output_parameters(f,optimization_results['hc_update']['extracted_parameters'])
    
    f.close()


"""
===========================================================================================================================
----------------------------------- Manual Hand Calculations --------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function to manually choose the Initial Circuit Parameters
# Inputs  : optimization_input_parameters
# Outputs : circuit_parameters, extracted_parameters
def manual_initial_parameters(cir,optimization_input_parameters):

    # Running Eldo
    cir.update_circuit(optimization_input_parameters['pre_optimization']['manual_circuit_parameters'].copy(),'basic')



"""
===========================================================================================================================
------------------------------------------- Output Functions --------------------------------------------------------------
"""

#--------------------------------------------------------------------------------------------------------------------------
# Function to perform pre-optimization
# Inputs  : mos_parameters, optimization_input_parameters, timing_results
# Outputs : circuit_parameters, extracted_parameters
def pre_optimization(cir,optimization_input_parameters,timing_results):

    # Opening the Run_Status File
    f=open(optimization_input_parameters['filename']['run_status'],'a')
    f.write('Pre Optimization Start\n Time : '+str(datetime.datetime.now())+'\n\n')
    f.close()

    # Storing the starting time
    timing_results['pre_optimization']=OrderedDict()
    timing_results['pre_optimization']['start']=datetime.datetime.now()

    print('************************************************************************************************************')
    print('*********************************** Pre Optimization *******************************************************')

    save_input_results_pre_optimization(optimization_input_parameters)

    cir.update_simulation_parameters(optimization_input_parameters['pre_optimization']['simulation'])

    optimization_results=OrderedDict()
    
    #======================================================== Manual Initial Points =============================================================================================================

    if optimization_input_parameters['pre_optimization']['type']=='manual':
        
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Manual Operating Point Selection ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #--------------------Initial Point Calculations-------------------------

        # Calculating the Values of Circuit Parameters
        manual_initial_parameters(cir,optimization_input_parameters)

        # Storing the Circuit and Extracted Parameters
        optimization_results['manual_hc']=OrderedDict()
        optimization_results['manual_hc']['circuit_parameters']=cir.circuit_parameters.copy()
        optimization_results['manual_hc']['extracted_parameters']=cir.extracted_parameters.copy()

        # Printing the values
        cff.print_circuit_parameters(cir.circuit_parameters)
        cff.print_extracted_parameters(cir.extracted_parameters)

        #cff.wait_key()

    
    #======================================================== Automatic Initial Points =============================================================================================================

    
    if optimization_input_parameters['pre_optimization']['type']==1:
        
        print('************************************************************************************************************')
        print('***********************************  Automatic Operating Point Selection 1 *********************************')

        # Extracting the MOSFET Parameters from the MOS file
        hc1.automatic_initial_parameters(cir,optimization_input_parameters,optimization_results)
    """
    if optimization_input_parameters['pre_optimization']['type']==2:

        print('************************************************************************************************************')
        print('***********************************  Automatic Operating Point Selection 2 *********************************')

        # Extracting the MOSFET Parameters from the MOS file
        hc2.automatic_initial_parameters(cir,optimization_input_parameters,optimization_results)
    
    if optimization_input_parameters['pre_optimization']['type']==3:

        print('************************************************************************************************************')
        print('***********************************  Automatic Operating Point Selection 3 *********************************')

        # Extracting the MOSFET Parameters from the MOS file
        hc3.automatic_initial_parameters(cir,optimization_input_parameters,optimization_results)
    """
    # Printing the values
    cff.print_circuit_parameters(cir.circuit_parameters)
    cff.print_extracted_parameters(cir.extracted_parameters)

    # Storing the results
    save_output_results_pre_optimization(optimization_results,optimization_input_parameters)

    # Storing the finishing time
    timing_results['pre_optimization']['stop']=datetime.datetime.now()

    # Opening the Run_Status File
    f=open(optimization_input_parameters['filename']['run_status'],'a')
    f.write('Pre Optimization End\n Time : '+str(datetime.datetime.now())+'\n\n')
    f.close()
    

#===========================================================================================================================
