#===========================================================================================================================
"""
Name                : Srivenkat A
Roll Number         : EE18B038
File Name           : ckt_optimizer.py
File Description    : This file will initialize the optimization_input_parameters and run the optimization file

Functions structure in this file:
    --> get_mos_parameters
    --> get_output_conditions
    --> get_simulation_conditions
    --> get_pre_optimization_parameters
    --> get_optimization_parameters
    --> get_temperature_analysis_parameters
    --> get_process_analysis_parameters
    
"""

#===========================================================================================================================
import numpy as np
import Ckt_Optimization.optimization as opt

#===========================================================================================================================
#------------------------------------ Other Functions ----------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the MOSFET parameters to the circuit_initialization_parameters dictionary
def get_mos_parameters(circuit_initialization_parameters,process_name):
    
    circuit_initialization_parameters['MOS']={}
    circuit_initialization_parameters['MOS']['Process']=process_name
    circuit_initialization_parameters['MOS']['filename']={}
    
    f=open('/home/ee18b038/Auto_Ckt_Synth_Codes/Automatic_RFIC_Synth/MOS_Files/'+process_name+'.txt')
    lines=f.readlines()
    f.close()

    # Extracting values from the MOS File
    for i in range(len(lines)):
        line=lines[i][:-1]
        if line=='Vdd':
            circuit_initialization_parameters['MOS']['Vdd']=float(lines[i+1][:-1])
        elif line=='Lmin':
            circuit_initialization_parameters['MOS']['Lmin']=float(lines[i+1][:-1])
        elif line=='u0':
            circuit_initialization_parameters['MOS']['un']=float(lines[i+1][:-1])*1e-4
        elif line=='tox':
            circuit_initialization_parameters['MOS']['tox']=float(lines[i+1][:-1])
        elif line=='vth0':
            circuit_initialization_parameters['MOS']['vt']=float(lines[i+1][:-1])
        elif line=='tt_file':
            circuit_initialization_parameters['MOS']['filename']['tt']=''
            j=i+1
            while lines[j][:-1]!='':
                circuit_initialization_parameters['MOS']['filename']['tt']+=lines[j]
                j+=1
        elif line=='ff_file':
            circuit_initialization_parameters['MOS']['filename']['ff']=''
            j=i+1
            while lines[j][:-1]!='':
                circuit_initialization_parameters['MOS']['filename']['ff']+=lines[j]
                j+=1
        elif line=='ss_file':
            circuit_initialization_parameters['MOS']['filename']['ss']=''
            j=i+1
            while lines[j][:-1]!='':
                circuit_initialization_parameters['MOS']['filename']['ss']+=lines[j]
                j+=1
                
    # Calculating Cox
    eo=8.85*1e-12
    er=3.9
    circuit_initialization_parameters['MOS']['cox']=eo*er/circuit_initialization_parameters['MOS']['tox']

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the output conditions to the optimization_input_parameters dictionary
def get_output_conditions(optimization_input_parameters,fo):
    
    optimization_input_parameters['output_conditions']={
        'gain_db':6.0,
        'op1db':10,
        'am-pm-dev':5,
        'wo':2.0*np.pi*fo,
        'delta_v':0.1,
        'Rin':50
    }

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the simulation conditions to the circuit_initialization_parameters dictionary
def get_simulation_conditions(circuit_initialization_parameters,fo):
    
    circuit_initialization_parameters['simulation']={}
    circuit_initialization_parameters['simulation']['standard_parameters']={}

    #file names
    circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']='/home/ee18b038/cadence_project/PA_tri3/'
    circuit_initialization_parameters['simulation']['standard_parameters']['run_directory']='/home/ee18b038/Auto_Ckt_Synth_Codes/Automatic_RFIC_Synth/'
    circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
    #circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
    circuit_initialization_parameters['simulation']['standard_parameters']['tcsh']=circuit_initialization_parameters['simulation']['standard_parameters']['run_directory']+'Spectre_Run/'
    #circuit_initialization_parameters['simulation']['standard_parameters']['iip3_type']='basic'     # 'basic' or 'advanced' 

    circuit_initialization_parameters['simulation']['standard_parameters']['std_temp']=27
    #circuit_initialization_parameters['simulation']['standard_parameters']['pin_fixed']=-65
    #circuit_initialization_parameters['simulation']['standard_parameters']['pin_start']=-70
    #circuit_initialization_parameters['simulation']['standard_parameters']['pin_stop']=-40
    #circuit_initialization_parameters['simulation']['standard_parameters']['pin_points']=6
    #circuit_initialization_parameters['simulation']['standard_parameters']['iip3_calc_points']=3
    circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']='tt'
    circuit_initialization_parameters['simulation']['standard_parameters']['conservative']='NO'
    circuit_initialization_parameters['simulation']['standard_parameters']['w_finger_max']=2e-6
    circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']=fo
    circuit_initialization_parameters['simulation']['standard_parameters']['f_range']=1e8

    circuit_initialization_parameters['simulation']['netlist_parameters']={
        'pin_start':-60,
        'pin_stop':20,
        'pin_step':1,
        'cir_temp':27,
        'n_harm':10
    }

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the pre_optimization parameters to the optimization_input_parameters dictionary
def get_pre_optimization_parameters(optimization_input_parameters,fo):

    optimization_input_parameters['pre_optimization']={}
    optimization_input_parameters['pre_optimization']['type']=1

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Manual Hand Calculations
    optimization_input_parameters['pre_optimization']['manual_circuit_parameters']={
        'Rin':50,
        'Rb':5000,
        'Rl':12.5,
        'Ld':10e-9,
        'C1':3.2e-12,
        'C2':1.25e-9,
        'W':420e-6,
        'Io':20e-3
    }
    
    optimization_input_parameters['pre_optimization']['Step1b_Limit']=5
    optimization_input_parameters['pre_optimization']['Step2_Limit']=5
    optimization_input_parameters['pre_optimization']['vdsat_reqd']=0.07

    optimization_input_parameters['pre_optimization']['gmrs_threshold']=0.2
    optimization_input_parameters['pre_optimization']['vdsat_threshold']=0.02
    optimization_input_parameters['pre_optimization']['C1_threshold']=200
    optimization_input_parameters['pre_optimization']['C2_threshold']=200
    optimization_input_parameters['pre_optimization']['Rbias_threshold']=100
    optimization_input_parameters['pre_optimization']['Rbias_minimum']=1000

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Pre Optimization Simulation Parameters
    optimization_input_parameters['pre_optimization']['simulation']={}
    optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']={}

    optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
    #optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
    #optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['iip3_type']='basic'
    optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['std_temp']=27
    #optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['pin_fixed']=-65
    #optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['pin_start']=-70
    #optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['pin_stop']=-40
    #optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['pin_points']=6
    #optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['iip3_calc_points']=3
    optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['process_corner']='tt'
    optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['conservative']='NO'

    optimization_input_parameters['pre_optimization']['simulation']['netlist_parameters']={
        'pin_start':-60,
        'pin_stop':20,
        'pin_step':1,
        'cir_temp':27,
        'n_harm':10
    }

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the optimization parameters to the optimization_input_parameters dictionary
def get_optimization_parameters(optimization_input_parameters,fo,optimization_name):

    optimization_input_parameters['optimization']={}

    optimization_input_parameters['optimization']['run']='YES'
    optimization_input_parameters['optimization']['n_runs']=1

    if optimization_name=='LOSS':
        optimization_input_parameters['optimization']['optimization_name']='loss1'
    #else:
    #    optimization_input_parameters['optimization']['optimization_name']='fom1'

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Parameters for Optimization Run 1
    optimization_input_parameters['optimization'][1]={}
    
    optimization_input_parameters['optimization'][1]['max_iteration']=300
    optimization_input_parameters['optimization'][1]['alpha_min']=-1
    optimization_input_parameters['optimization'][1]['consec_iter']=-1

    optimization_input_parameters['optimization'][1]['delta_threshold']=0.001
    optimization_input_parameters['optimization'][1]['alpha_mult']=1
    optimization_input_parameters['optimization'][1]['loss_type']=0
    optimization_input_parameters['optimization'][1]['update_check']=0

    optimization_input_parameters['optimization'][1]['optimizing_parameters']=['Rb','Rl','Ld','Io','W','C1','C2']
    optimization_input_parameters['optimization'][1]['output_parameters_list']=['Isup_hb','Ids_dc','Ids_hb','gain_db','p_source','Psup_dc','op1db_auto','am-pm-dev','ip1db_auto','gm','Voutdc','vdsat','vg','vd','gds']
    
    optimization_input_parameters['optimization'][1]['optimization_type']=0

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Assigning values to the loss weights
    loss_weights={}
    loss_weights['gain_db']=1/6.0  
    loss_weights['op1db']=1/10.0  
    loss_weights['am-pm-dev']=1/5.0   
    loss_weights['Isup']=1000 
    optimization_input_parameters['optimization'][1]['loss_weights']=loss_weights


    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Assigning Values of Alpha
    alpha_parameters={}
    alpha_parameters['common']=0.2
    alpha_parameters['Rb']=1
    alpha_parameters['Ld']=1
    alpha_parameters['Rl']=1
    alpha_parameters['W']=1
    alpha_parameters['Io']=1
    alpha_parameters['C1']=1
    alpha_parameters['C2']=1
    optimization_input_parameters['optimization'][1]['alpha']={}
    optimization_input_parameters['optimization'][1]['alpha']['values']=alpha_parameters

    optimization_input_parameters['optimization'][1]['alpha']['type']='Normal'
    optimization_input_parameters['optimization'][1]['alpha']['start']=0.8
    optimization_input_parameters['optimization'][1]['alpha']['end']=0.05

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Optimization Simulation Parameters
    optimization_input_parameters['optimization']['simulation']={}
    optimization_input_parameters['optimization']['simulation'][1]={}
    optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']={}

    optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
    #optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
    #optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['iip3_type']='basic'
    optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['std_temp']=27
    #optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['pin_fixed']=-65
    #optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['pin_start']=-70
    #optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['pin_stop']=-40
    #optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['pin_points']=16
    #optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['iip3_calc_points']=5
    optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['process_corner']='tt'
    optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['conservative']='NO'

    optimization_input_parameters['optimization']['simulation'][1]['netlist_parameters']={
        'pin_start':-60,
        'pin_stop':20,
        'pin_step':1,
        'cir_temp':27,
        'n_harm':10
    }

    """
    optimization_input_parameters['optimization']['simulation'][2]={}
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']={}

    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65'
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65'
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['iip3_type']='basic'
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['std_temp']=27
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['pin_fixed']=-65
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['pin_start']=-70
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['pin_stop']=-40
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['pin_points']=16
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['iip3_calc_points']=5
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['process_corner']='tt'
    optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['conservative']='YES'

    optimization_input_parameters['optimization']['simulation'][2]['netlist_parameters']={
        'pin':-65,
        'fund_2':fo+1e6,
        'fund_1':fo,
        'cir_temp':27,
        'n_harm':15
    }
    """
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Conditions for acceptable solution 
    """
    if optimization_name=='FOM':
        optimization_input_parameters['acceptable_solution']={}
        optimization_input_parameters['acceptable_solution']['s11_db']=-15
        optimization_input_parameters['acceptable_solution']['gain_db']=6
        optimization_input_parameters['acceptable_solution']['iip3_dbm']=-15
        optimization_input_parameters['acceptable_solution']['nf_db']=10
        optimization_input_parameters['acceptable_solution']['p_source']=10e-3
    """
        
#---------------------------------------------------------------------------------------------------------------------------
"""
# Function that sets the temperature analysis parameters to the optimization_input_parameters dictionary

def get_temperature_analysis_parameters(optimization_input_parameters,fo):

    optimization_input_parameters['temperature_analysis']={}

    optimization_input_parameters['temperature_analysis']['run']='YES'

    optimization_input_parameters['temperature_analysis']['start_temp']=-40
    optimization_input_parameters['temperature_analysis']['stop_temp']=120
    optimization_input_parameters['temperature_analysis']['n_temp']=5

    optimization_input_parameters['temperature_analysis']['start_current']=0.1
    optimization_input_parameters['temperature_analysis']['stop_current']=10
    optimization_input_parameters['temperature_analysis']['n_current']=1


    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Temperature Analysis Simulation Parameters
    optimization_input_parameters['temperature_analysis']['simulation']={}
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']={}

    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['iip3_type']='basic'
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['std_temp']=27
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['pin_fixed']=-65
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['pin_start']=-70
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['pin_stop']=-40
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['pin_points']=16
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['iip3_calc_points']=5
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['process_corner']='tt'
    optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['conservative']='YES'

    optimization_input_parameters['temperature_analysis']['simulation']['netlist_parameters']={
        'pin':-65,
        'fund_2':fo+1e6,
        'fund_1':fo,
        'cir_temp':27,
        'n_harm':15
    }

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the process analysis parameters to the optimization_input_parameters dictionary
def get_process_analysis_parameters(optimization_input_parameters,fo):

    optimization_input_parameters['process_analysis']={}
    optimization_input_parameters['process_analysis']['run']='YES'

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Temperature Analysis Simulation Parameters
    optimization_input_parameters['process_analysis']['simulation']={}
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']={}

    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['iip3_type']='basic'
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['std_temp']=27
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['pin_fixed']=-65
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['pin_start']=-70
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['pin_stop']=-40
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['pin_points']=16
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['iip3_calc_points']=5
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['process_corner']='tt'
    optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['conservative']='YES'

    optimization_input_parameters['process_analysis']['simulation']['netlist_parameters']={
        'pin':-65,
        'fund_2':fo+1e6,
        'fund_1':fo,
        'cir_temp':27,
        'n_harm':15
    }

"""
#===========================================================================================================================
#------------------------------------Main Program Code----------------------------------------------------------------------

# Creating a dictionary with the optimization parameters
circuit_initialization_parameters={}
optimization_input_parameters={}
optimization_name='LOSS'

# ---------- MOSFET Parameters ----------
#get_mos_parameters(circuit_initialization_parameters,'TSMC180')
get_mos_parameters(circuit_initialization_parameters,'TSMC65')
#get_mos_parameters(circuit_initialization_parameters,'TSMC65_2')
#get_mos_parameters(circuit_initialization_parameters,'IBM130')

# ---------- Output Conditions ----------
fo=2e9
get_output_conditions(optimization_input_parameters,fo)

# ---------- Simulation Conditions ----------
get_simulation_conditions(circuit_initialization_parameters,fo)

# ---------- Pre Optimization Parameters ----------
get_pre_optimization_parameters(optimization_input_parameters,fo)

# ---------- Optimization Parameters ----------
get_optimization_parameters(optimization_input_parameters,fo,optimization_name)

# ---------- Temperature Analysis Parameters ----------
#get_temperature_analysis_parameters(optimization_input_parameters,fo)

# ---------- Temperature Analysis Parameters ----------
#get_process_analysis_parameters(optimization_input_parameters,fo)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------- FILE NAMES ------------------------------------------

optimization_input_parameters['filename']={}
optimization_input_parameters['filename']['run_status']='/home/ee18b038/Auto_Ckt_Synth_Codes/Simulation_Results/run_status.txt'

f_directory='/home/ee18b038/Auto_Ckt_Synth_Codes/Simulation_Results/'+str(optimization_name)+'/'


file_choose='S' # 'S' to run a single time; 'M' to run multiple times


if file_choose=='S':

    # ------- Set Any Additional Parameters Here --------
    filename=f_directory+'Test_Simulation_Parameters'                       # SET THE FILENAME HERE
    # ------- Set Any Additional Parameters Here --------
    

    # ------- DON'T CHANGE THESE LINES -------------
    optimization_input_parameters['filename']['output']=filename
    opt.complete_optimization(circuit_initialization_parameters,optimization_input_parameters,'PA')          
    # ------- DON'T CHANGE THESE LINES -------------        


if file_choose=='M':
    for i in range(0,5):    # SET NUMBER OF ITERATIONS HERE

        # ------- Set Any Additional Parameters Here --------
        filename=f_directory+'tsmc_str_1'+str(i)            # SET THE FILENAME HERE
        fo=1e9
        wo=2*np.pi*fo
        optimization_input_parameters['output_conditions']['wo']=wo
        iip_mtd=['basic_hb','basic_pss','advanced_hb','advanced_pss','hb_sweep']
        optimization_input_parameters['simulation']['iip3_method']=iip_mtd[i]
        
        # ------- Set Any Additional Parameters Here --------


        # ------- DON'T CHANGE THESE LINES -------------
        optimization_input_parameters['filename']['output']=filename
        opt.complete_optimization(optimization_input_parameters)
        # ------- DON'T CHANGE THESE LINES -------------

#===========================================================================================================================
