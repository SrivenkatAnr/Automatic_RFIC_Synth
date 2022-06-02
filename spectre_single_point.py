#===========================================================================================================================
"""
Name				: Srivenkat A
Roll Number			: EE18B038
File Name			: spectre_single_point.py
File Description 	: This file will run spectre and extract the parameters for a single point

Functions structure in this file:
	--> get_mos_parameters
	--> get_simulation_conditions
	
"""

#===========================================================================================================================

import PA.spectre as sp
import time

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
# Function that sets the simulation conditions to the optimization_input_parameters dictionary
def get_simulation_conditions_PA(circuit_initialization_parameters,fo):
	
	circuit_initialization_parameters['simulation']={}
	circuit_initialization_parameters['simulation']['standard_parameters']={}

	circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']='/home/ee18b038/cadence_project/PA_single_pt_2/'
	circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_tsmc_65_rcm'
	circuit_initialization_parameters['simulation']['standard_parameters']['run_directory']='/home/ee18b038/Auto_Ckt_Synth_Codes/Automatic_RFIC_Synth/'
	circuit_initialization_parameters['simulation']['standard_parameters']['tcsh']=circuit_initialization_parameters['simulation']['standard_parameters']['run_directory']+'Spectre_Run/'
	circuit_initialization_parameters['simulation']['standard_parameters']['std_temp']=27
	circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']='tt'
	circuit_initialization_parameters['simulation']['standard_parameters']['conservative']='NO'
	circuit_initialization_parameters['simulation']['standard_parameters']['w_finger_max']=2e-6
	circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']=fo
	circuit_initialization_parameters['simulation']['standard_parameters']['f_range']=1e8

	circuit_initialization_parameters['simulation']['standard_parameters']['C1_threshold']=2000
	circuit_initialization_parameters['simulation']['standard_parameters']['C2_threshold']=100
	circuit_initialization_parameters['simulation']['standard_parameters']['Rbias_threshold']=100
	circuit_initialization_parameters['simulation']['standard_parameters']['Rbias_minimum']=1000

	circuit_initialization_parameters['simulation']['netlist_parameters']={
        'pin_start':-25,
        'pin_stop':15,
        'pin_step':0.5,
        'cir_temp':27,
        'n_harm':10
	}


#===========================================================================================================================
#------------------------------------Main Program Code----------------------------------------------------------------------

t_start=time.time()

# Creating a dictionary with the optimization parameters
circuit_initialization_parameters={}

# ---------- MOSFET Parameters ----------
#get_mos_parameters(circuit_initialization_parameters,'TSMC65')
get_mos_parameters(circuit_initialization_parameters,'TSMC65_2')
#get_mos_parameters(circuit_initialization_parameters,'TSMC65_3')
#get_mos_parameters(circuit_initialization_parameters,'TSMC65_4')

# ---------- Simulation Conditions ----------
fo=2e9
get_simulation_conditions_PA(circuit_initialization_parameters,fo)

circuit_parameters={
	'Rin':50,
	'Rb':2359,
	'Rl':100,
    'Rl_expected':100,
	'Ld':5.99e-9,
	'Lsrc':2.274e-9,
	'Lload':2.6526e-9,
	'Cmn':2.92e-12,
	#'C1':89.7e-12,
	#'C2':35.5e-9,
	'W':1.38e-3,
	'Io':17.8e-3
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------- FILE RUN --------------------------------------------
cir=sp.Circuit(circuit_initialization_parameters)
cir.update_circuit(circuit_parameters,'advanced')

print ('____________________________________________________________________')
print ('------------------------Circuit Parameters------------------------\n')
for param_name in cir.circuit_parameters:
	print(param_name,' : ',cir.circuit_parameters[param_name])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------- OUTPUT PRINT --------------------------------------------
print ('____________________________________________________________________')
print ('------------------------Extracted Parameters------------------------\n')
for param_name in cir.extracted_parameters:
    if ('comb_' not in param_name):
	    print(param_name,' : ',cir.extracted_parameters[param_name])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------- SAVE PLOTS --------------------------------------------
cir.plot_ckt_trends("Output_Trends/")

t_end = time.time()
print("\n time taken is {} seconds\n".format(t_end-t_start))

#===========================================================================================================================
