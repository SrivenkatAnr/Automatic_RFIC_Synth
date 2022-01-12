#===========================================================================================================================
"""
Name                : Srivenkat A
Roll Number         : EE18B038
File Name           : spectre.py
File Description    : This file will contain the functions to write, run, and read from the spectre files

Functions structure in this file:
    --> valueName_to_value
    --> valueE_to_value
    --> extract_file
    --> extract_basic_parameters
        --> extract_dc_param
        --> extract_ac_param
        --> extract_sp_param
        --> extract_noise_param
    --> print_param
    --> dict_convert
    --> write_circuit_parameters
    --> write_MOS_parameters
    --> write_simulation_parameters
    --> write_tcsh_file

    --> write_extract
        --> write_extract_basic
        --> write_extract_iip3
        --> run_file


    
"""
#===========================================================================================================================
import numpy as np
import fileinput
import os
import multiprocessing as mp
import common_functions as cff # type: ignore
import copy

"""
====================================================================================================================================================================================
------------------------------------------------------------ CIRCUIT CLASS ---------------------------------------------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Creating a class for the circuit 
class Circuit():
    def __init__(self,circuit_initialization_parameters):
        self.circuit_parameters={}
        self.extracted_parameters={}
        self.simulation_parameters={}
        self.circuit_initialization_parameters=circuit_initialization_parameters
        self.mos_parameters=calculate_mos_parameters(self.circuit_initialization_parameters)
    
    def run_circuit(self):
        #if self.circuit_parameters['W']>900e-6:
        #   self.circuit_parameters['W']=850e-6
        self.extracted_parameters=write_extract(self.circuit_parameters,self.circuit_initialization_parameters)

    def update_circuit(self,circuit_parameters):
        self.circuit_parameters=circuit_parameters
        #if self.circuit_parameters['W']>900e-6:
        #   self.circuit_parameters['W']=850e-6
        self.extracted_parameters=write_extract(circuit_parameters,self.circuit_initialization_parameters)
    
    def update_circuit_parameters(self,circuit_parameters):
        self.circuit_parameters=circuit_parameters

    def update_simulation_parameters(self,simulation_parameters):
        if 'netlist_parameters' in simulation_parameters:
            for param_name in simulation_parameters['netlist_parameters']:
                self.circuit_initialization_parameters['simulation']['netlist_parameters'][param_name]=simulation_parameters['netlist_parameters'][param_name]
        
        if 'standard_parameters' in simulation_parameters:
            for param_name in simulation_parameters['standard_parameters']:
                self.circuit_initialization_parameters['simulation']['standard_parameters'][param_name]=simulation_parameters['standard_parameters'][param_name]
    
    def write_simulation_parameters(self):
        write_simulation_parameters(self.circuit_initialization_parameters)
    
    def update_temp(self,temp):
        self.circuit_initialization_parameters['simulation']['netlist_parameters']['cir_temp']=temp
    
    def reset_temp(self):
        self.circuit_initialization_parameters['simulation']['netlist_parameters']['cir_temp']=self.circuit_initialization_parameters['simulation']['standard_parameters']['std_temp']  



#===========================================================================================================================
#------------------------------------ MOSFET EXTRACTION --------------------------------------------------------------------

#-----------------------------------------------------------------
# Function that extracts the MOSFET File Parameeters
# Inputs  : Optimization Input Parameters
# Outputs : MOS_Parameters
def calculate_mos_parameters(circuit_initialization_parameters):
    
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


"""
====================================================================================================================================================================================
------------------------------------------------------------ EXTRACTION FUNCTION ---------------------------------------------------------------------------------------------------
"""

#===========================================================================================================================================================
#------------------------------------------------------ Character to Real Number Functions -----------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Changing the values extracted as a string to a floating point value 
# Input: Value of the number in string format   
# Output: Value of the number in float
def valueName_to_value(value_name):

    # Checking if the last character of array is a string
    if value_name[-1].isalpha()==0:
        val=float(value_name)
        return val
    
    # Checking if the last character of array is a string
    if (value_name[-1]=='G' and value_name[-2]=='E') or (value_name[-1]=='g' and value_name[-2]=='e'):
        val=float(value_name[:-3])*1e6
        return val
        
    # Extracting the numerical part of the number 
    val=float(value_name[:-1])
    
    # Extracting the character that denotes the units ( i.e, millt, micro, nano, etc)
    mult_name=value_name[-1]
    mult=1.0
    
    # Calculating the value of the unit
    if mult_name=='M' or mult_name=='m':
        mult=1e-3
    elif mult_name=='U' or mult_name=='u':
        mult=1e-6
    elif mult_name=='N' or mult_name=='n':
        mult=1e-9
    elif mult_name=='P' or mult_name=='p':
        mult=1e-12
    elif mult_name=='F' or mult_name=='f':
        mult=1e-15
    elif mult_name=='G' or mult_name=='g':
        mult=1e9
    else:
        mult=1.0
        
    val=val*mult
    return val
    
#---------------------------------------------------------------------------------------------------------------------------
# Changing the values extracted as 10e1, 1.5e-2 to a floating point value 
# Input: Value of the number in string format   
# Output: Value of the number in float
def valueE_to_value(value_name):
    
    # Extracting the number before and after e
    if 'e' in value_name:
        num1=float(value_name.split('e')[0])
        num2=float(value_name.split('e')[1])
        
        # Calculating the final number
        num=num1*(10**num2)
    
    else:
        num=float(value_name)
    
    return num


#===========================================================================================================================================================
#--------------------------------------------------------- Other File Extraction Functions -----------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Extracting the files as an array of lines
# Inputs: file name
# Output: array of lines
def extract_file(file_name):
    f=open(file_name)
    lines=f.readlines()
    f.close()
    return lines


#===========================================================================================================================================================
#------------------------------------------------------ Basic File Extraction Functions --------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------    
# Extracting the DC from the file
# Inputs: circuit_initialization_parameters
# Output: Dictionary with all the parameters
def extract_dc_param(circuit_initialization_parameters):

    # Getting the filename
    file_name=circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/dc.out'
    lines=extract_file(file_name)
    
    extracted_parameters={}

    # Skipping the first few lines
    lines=lines[7:]
    lines=lines[0].split()

    # Extracting the values from the required line
    extracted_parameters['vg1']=valueE_to_value(lines[5])
    extracted_parameters['vd1']=valueE_to_value(lines[6])
    
    extracted_parameters['i_source']=np.absolute(valueE_to_value(lines[4]))
    extracted_parameters['v_source']=np.absolute(valueE_to_value(lines[3]))
    extracted_parameters['p_source']=extracted_parameters['i_source']*extracted_parameters['v_source']
    extracted_parameters['Voutdc']=valueE_to_value(lines[1])
    extracted_parameters['Isup']=valueE_to_value(lines[2])
    extracted_parameters['Vsup']=valueE_to_value(lines[16])
    extracted_parameters['Psup']=extracted_parameters['Isup']*extracted_parameters['Vsup']
    extracted_parameters['gm1']=valueE_to_value(lines[8])
    extracted_parameters['gds1']=valueE_to_value(lines[9])
    extracted_parameters['vth1']=valueE_to_value(lines[10])
    extracted_parameters['vdsat1']=valueE_to_value(lines[12])
    extracted_parameters['Idsat1']=valueE_to_value(lines[7])

    extracted_parameters['cgs1']=np.absolute(valueE_to_value(lines[13]))
    extracted_parameters['cgd1']=np.absolute(valueE_to_value(lines[14]))
    extracted_parameters['region']=np.absolute(valueE_to_value(lines[15]))

    
    return extracted_parameters

#---------------------------------------------------------------------------------------------------------------------------    
# Extracting the AC from the file
# Inputs: circuit_initialization_parameters
# Output: Dictionary with all the parameters
def extract_ac_param(circuit_initialization_parameters):

    # Getting the filename
    file_name=circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/ac.out'
    lines=extract_file(file_name)
    
    extracted_parameters={}

    # Skipping the first few lines
    lines=lines[7:]
    lines=lines[0].split()

    # Extracting the values frim the required line
    extracted_parameters['freq']=valueE_to_value(lines[0])
    vout_re=valueE_to_value(lines[1])
    vout_im=valueE_to_value(lines[2])
    vin_re=valueE_to_value(lines[3])
    vin_im=valueE_to_value(lines[4])
    extracted_parameters['gain_db']=calculate_gain_db(vout_re,vout_im,vin_re,vin_im)
    extracted_parameters['gain_phase']=calculate_gain_phase(vout_re,vout_im,vin_re,vin_im)
    
    return extracted_parameters

#---------------------------------------------------------------------------------------------------------------------------    
# Extracting the AC from the file
# Inputs: circuit_initialization_parameters
# Output: Dictionary with all the parameters
def extract_xdb_auto_param(circuit_initialization_parameters):

    # Getting the filename
    file_name=circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/circ.raw/gcomp_test.xdb.pss_hb'
    lines=extract_file(file_name)

    extracted_parameters={}
    ip1db=0
    op1db=0

    for line in lines:
        if '"Pin"' in line:
            ip1db = line.split()[1]
        if '"Power"' in line:
            op1db = line.split()[1]

    extracted_parameters['ip1db_auto']=valueE_to_value(ip1db)
    extracted_parameters['op1db_auto']=valueE_to_value(op1db)

    return extracted_parameters

#---------------------------------------------------------------------------------------------------------------------------    
# Extracting the AC from the file
# Inputs: circuit_initialization_parameters
# Output: Dictionary with all the parameters
def extract_comp_param(circuit_initialization_parameters,extracted_parameters_xdb):

    # Getting the filename
    fname_template=circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/circ.raw/swp-{}_phdev_test.fd.pss_hb'
    
    pin_start=circuit_initialization_parameters['simulation']['netlist_parameters']['pin_start']
    pin_stop=circuit_initialization_parameters['simulation']['netlist_parameters']['pin_stop']
    pin_step=circuit_initialization_parameters['simulation']['netlist_parameters']['pin_step']
    npin=int((pin_stop-pin_start)/pin_step)
    op_freq=circuit_initialization_parameters['simulation']['netlist_parameters']['fund_1']

    vin_re_arr = np.zeros(npin,dtype=float)
    vin_im_arr = np.zeros(npin,dtype=float)
    vout_re_arr = np.zeros(npin,dtype=float)
    vout_im_arr = np.zeros(npin,dtype=float)
    ph_arr = np.zeros(npin,dtype=float)
    gdb_arr = np.zeros(npin,dtype=float)
    extracted_parameters={}

    freq_flag=0

    for i in range(npin):
        if (i==0):
            filename=fname_template.replace("{}","000")  
        elif (i<=9):
            filename=str(fname_template.format("00"+str(i)))
        elif (i<=99):
            filename=str(fname_template.format("0"+str(i)))
        else:
            filename=str(fname_template.format(i))
        for line in extract_file(filename):
            if ('"freq"' in line) and ('"sweep"' not in line):
                freq=line.split()[1]
                freq=valueE_to_value(freq)
                if(freq==op_freq):
                    freq_flag=1
            if ('"Vin"' in line) and ('"V"' not in line) and (freq_flag==1):
                vin_re_arr[i],vin_im_arr[i] = extract_voltage(line)
            if ('"Vout"' in line) and ('"V"' not in line) and (freq_flag==1):
                vout_re_arr[i],vout_im_arr[i] = extract_voltage(line)
                freq_flag=0
        ph_arr[i] = calculate_gain_phase(vout_re_arr[i], vout_im_arr[i], vin_re_arr[i], vin_im_arr[i])
        gdb_arr[i] = calculate_gain_db(vout_re_arr[i], vout_im_arr[i], vin_re_arr[i], vin_im_arr[i])

    ip1db = extracted_parameters_xdb['ip1db_auto']
    ip_index=int((ip1db-pin_start)/pin_step)
    ph_min=min(ph_arr[:ip_index+1])
    ph_max=max(ph_arr[:ip_index+1]) 
    extracted_parameters["am-pm-dev"]=(ph_max-ph_min)

    gdb_ss=max(gdb_arr)
    ip_man_index=np.argmin(abs(gdb_arr-(gdb_ss-1)))
    extracted_parameters["ip1db_man"]=pin_start+(ip_man_index*pin_step)

    return extracted_parameters    

#---------------------------------------------------------------------------------------------------------------------------    
# Extracts Vout_magnitude from hb,pss file line
# Inputs: Line
# Output: Vout_Magnitude
def extract_voltage(lines):
    
    # Extracting Vout Magnitude
    lines=lines.split()
    char_r=lines[1].split('(')[1]
    char_i=lines[2].split(')')[0]

    # Converting string to floating point value
    vout_r=valueE_to_value(char_r)
    vout_i=valueE_to_value(char_i)

    return (vout_r,vout_i)

#---------------------------------------------------------------------------------------------------------------------------    
# Calculating the gain and angle from the vout and vin values
# Inputs: vout and vin
# Output: gain_db and phase
def calculate_gain_db(vout_re,vout_im,vin_re,vin_im):
    
    # Calculating gain_dB
    gain=(vout_re**2+vout_im**2)/(vin_re**2+vin_im**2)
    gain_db=10*np.log10(gain)

    return gain_db

#---------------------------------------------------------------------------------------------------------------------------    
# Calculating the gain and angle from the vout and vin values
# Inputs: vout and vin
# Output: gain_db and phase
def calculate_gain_phase(vout_re,vout_im,vin_re,vin_im):
    
    vout=np.complex(vout_re,vout_im)
    vin=np.complex(vin_re,vin_im)
    phase=np.angle(vout/vin)*180/np.pi

    return phase
#---------------------------------------------------------------------------------------------------------------------------
# Extracting all the output parameters from chi file
# Inputs: optimization_input parameters
# Outputs: output parameters dictionary 
def extract_basic_parameters(circuit_initialization_parameters):
    
    # Extracting the outputs 
    extracted_parameters_dc=extract_dc_param(circuit_initialization_parameters)
    extracted_parameters_ac=extract_ac_param(circuit_initialization_parameters)
    extracted_parameters_xdb=extract_xdb_auto_param(circuit_initialization_parameters)
    extracted_parameters_comp=extract_comp_param(circuit_initialization_parameters,extracted_parameters_xdb)
    
    # Storing the outputs in a single dictionary
    extracted_parameters={}
    
    for param_name in extracted_parameters_dc:
        extracted_parameters[param_name]=extracted_parameters_dc[param_name]
    for param_name in extracted_parameters_ac:
        extracted_parameters[param_name]=extracted_parameters_ac[param_name]
    for param_name in extracted_parameters_xdb:
        extracted_parameters[param_name]=extracted_parameters_xdb[param_name]
    for param_name in extracted_parameters_comp:
        extracted_parameters[param_name]=extracted_parameters_comp[param_name]

    return extracted_parameters


#===========================================================================================================================================================
#------------------------------------------- AM-AM, AM-PM Deviation Extraction Functions --------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------    

"""
#---------------------------------------------------------------------------------------------------------------------------    
# Calculating the IIP3 after extraction of Vout data
# Inputs: circuit_initialization_parameters, Vout_fund, Vout_im3, pin
# Output: IIP3
def calculate_iip3_multiple_points(circuit_initialization_parameters,vout_fund_mag,vout_im3_mag,pin):

    # Calculating values in log scale
    vout_fund_log=20*np.log10(vout_fund_mag)
    vout_im3_log=20*np.log10(vout_im3_mag)

    n_pin=circuit_initialization_parameters['simulation']['standard_parameters']['pin_points']
    n_points=circuit_initialization_parameters['simulation']['standard_parameters']['iip3_calc_points']
    
    # Creating arrays for slopes and y-intercepts of fundamental and im3 components
    fund_slope=np.zeros(n_pin+1-n_points,dtype=float)
    fund_intercept=np.zeros(n_pin+1-n_points,dtype=float)
    im3_slope=np.zeros(n_pin+1-n_points,dtype=float)
    im3_intercept=np.zeros(n_pin+1-n_points,dtype=float)

    # Calculating the slopes and y-intercepts
    for i in range(n_pin+1-n_points):
        fund_slope[i],fund_intercept[i]=calculate_slope(pin[i:i+n_points-1],vout_fund_log[i:i+n_points-1])
        im3_slope[i],im3_intercept[i]=calculate_slope(pin[i:i+n_points-1],vout_im3_log[i:i+n_points-1])
    
    # Finding the best points for iip3 calculation
    best_point=calculate_best_iip3_point(fund_slope,im3_slope)
    
    # Calculating the iip3 given the slope and y-intercept of fundamental and im3
    iip3=(im3_intercept[best_point]-fund_intercept[best_point])/(fund_slope[best_point]-im3_slope[best_point])

    return iip3

#---------------------------------------------------------------------------------------------------------------------------    
# Calculating the slope and y-intercept
# Inputs: x and y coordinates of the points
# Output: slope, y-intercept
def calculate_slope(x,y):
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    return m,c

#---------------------------------------------------------------------------------------------------------------------------    
# Calculating the point with slope closest to 1dB/dB for fund and 3dB/dB for im3
# Inputs: Slope of fundamental and im3
# Output: Location of the best point
def calculate_best_iip3_point(fund_slope,im3_slope):
    
    # Getting the length of the array
    l=len(fund_slope)

    # Choosing the best point as the first point
    best_point=0
    best_error=(fund_slope[0]-1)**2+(im3_slope[0]-3)**2

    for i in range(1,l):
        error=(fund_slope[i]-1)**2+(im3_slope[i]-3)**2
        if error<best_error:    # Checking if the current point is better than the best point
            best_point=i
            best_error=error
    return best_point

#---------------------------------------------------------------------------------------------------------------------------    
# Checks if the frequency is within range ( within (target-error,target+error) )
# Inputs: Test Frequency, Target Frequency, Error
# Output: 1 if Yes and 0 if No
def check_freq(f_test,f_target,f_error):
    if f_test<f_target+f_error and f_test>f_target-f_error:
        return 1
    else:
        return 0

#---------------------------------------------------------------------------------------------------------------------------    
# Extracting Vout magnitude of fundamental and im3 from file ( for hb_sweep )
# Inputs: Filename, Optimization Input Parameters
# Output: Magnitude of Vout at fundamental and im3
def extract_vout_magnitude(file_name,circuit_initialization_parameters):

    lines=extract_file(file_name)
    
    fund_1=circuit_initialization_parameters['simulation']['netlist_parameters']['fund_1']
    fund_2=circuit_initialization_parameters['simulation']['netlist_parameters']['fund_2']

    f_im3=2*fund_2-fund_1
    f_error=(fund_2-fund_1)/100

    flag=0
    flag_fun=0
    flag_im3=0
    flag_test=0

    while 1:
        if len(lines[0].split())<2:
            lines=lines[1:]
        
        elif 'freq' in lines[0].split()[0] and flag==0:
            flag=1
            lines=lines[1:]
        
        elif 'freq' in lines[0].split()[0] and flag==1:
            if flag_fun==0 and check_freq(float(lines[0].split()[1]),fund_2,f_error)==1 :
                
                #Extracting Vout for fundamental
                flag_test=1
                while flag_test==1:
                    if 'Vout' in lines[0].split()[0]:
                        flag_test=0
                        vout_fund=extract_vout(lines[0])
                    else:
                        lines=lines[1:]
                flag_fun=1
            
            elif flag_im3==0 and check_freq(float(lines[0].split()[1]),f_im3,f_error)==1 :
                
                #Extracting Vout for im3
                flag_test=1
                while flag_test==1:
                    if 'Vout' in lines[0].split()[0]:
                        flag_test=0
                        vout_im3=extract_vout(lines[0])
                    else:
                        lines=lines[1:]
                flag_im3=1
            lines=lines[1:]
            
            if flag_fun==1 and flag_im3==1:
                break
        
        else:
            lines=lines[1:]

    return vout_fund,vout_im3

#---------------------------------------------------------------------------------------------------------------------------    
# Extracts Vout_magnitude from hb,pss file line
# Inputs: Line
# Output: Vout_Magnitude
def extract_vout(lines):
    
    # Extracting Vout Magnitude
    lines=lines.split()
    char_r=lines[1].split('(')[1]
    char_i=lines[2].split(')')[0]

    # Converting string to floating point value
    vout_r=valueE_to_value(char_r)
    vout_i=valueE_to_value(char_i)
    
    # Calculating the magnitude of the output
    vout_mag=np.sqrt(vout_r*vout_r+vout_i*vout_i)

    return vout_mag


#===========================================================================================================================
"""

"""
====================================================================================================================================================================================
------------------------------------------------------------ FILE WRITE FUNCTIONS --------------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Command that returns the string that has to be printed in the .scs file
# Inputs  : Name of the parameter, Value of the parameter
# Outputs : String to be printed
def print_param(param_var,val):
    return "parameters "+param_var+'='+str(val)+'\n'

#-----------------------------------------------------------------      
# Function that converts input parameter dictionary to writing dictionary
# Inputs  : Circuit Parameters Dictionary, Optimization Input Parameters
# Outputs : The dictionary containing the parameters to be written to the .scs file
def dict_convert(circuit_parameters,circuit_initialization_parameters):
    write_dict={}
    
    # param_names in write_dict will contain the name of the parameters as it is written in the .scs file
    cir_writing_dict={
        'wid_bias':'W',
        'cur0':'Io',
        'res_bias':'Rb',
        'res_in':'Rin',
        'res_drain':'Rd',
        'res_load':'Rl',
        'ind_drain':'Ld',
        'cap_coup_in':'C1',
        'cap_coup_out':'C2',
    }
    for param_name in cir_writing_dict:
        try:
            write_dict[param_name]=circuit_parameters[cir_writing_dict[param_name]]
        except:
            pass
        
    # Getting the value of resistances
    write_dict['res_drain']=circuit_parameters['Ld']*circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']*2*np.pi/50
    
    # Calculating the number of fingers
    n_finger=int(circuit_parameters['W']/circuit_initialization_parameters['simulation']['standard_parameters']['w_finger_max'])+1
    write_dict['n_finger']=n_finger

    return write_dict
            
#-----------------------------------------------------------------
# Function that modifies the .scs file
# Inputs  : circuit_parameters, optimization input parameters
# Outputs : NONE
def write_circuit_parameters(circuit_parameters,circuit_initialization_parameters):
    
    # We will convert the circuit parameters to write_dict
    write_dict=dict_convert(circuit_parameters,circuit_initialization_parameters)
    
    # Getting the filenames
    filename1=circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/circ.scs'

    # We will write the new values to the Basic Circuit
    f=open(filename1,'r+') 
    s=''
    for line in fileinput.input(filename1):
        for param_name in write_dict:
            if "parameters "+param_name+'=' in line:    # Checking for a particular parameter in the .scs file
                line=line.replace(line,print_param(param_name,write_dict[param_name]))  # Replacing the parameter in the .scs file
        s=s+line

    f.truncate(0)
    f.write(s)
    f.close()

#-----------------------------------------------------------------
# Function that adds Simulation Parameters
# Inputs  : Optimization Input Parameters
# Outputs : NONE
def write_simulation_parameters(circuit_initialization_parameters):
    
    # Adding simulation_parameters to write_dict
    write_dict={}
    for param_name in circuit_initialization_parameters['simulation']['netlist_parameters']:
        write_dict[param_name]=circuit_initialization_parameters['simulation']['netlist_parameters'][param_name]
    process_corner=circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']
    
    write_dict['len']=circuit_initialization_parameters['MOS']['Lmin']
    write_dict['v_dd']=circuit_initialization_parameters['MOS']['Vdd']

    # Getting the filenames
    filename1=circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/circ.scs'

    # Writing the simulation parameters to Basic File
    f=open(filename1,'r+')
    s=''
    write_check=1
    include_check=0
    
    # Replacing the lines of .scs file
    for line in fileinput.input(filename1):
        if "include " in line:  # This line is used to include the MOS file in the .scs file
            include_check=0
            write_check=1

        elif "include" not in line and include_check==1:
            s=s+circuit_initialization_parameters['MOS']['filename'][process_corner]
            include_check=0
            write_check=1
        
        for param_name in write_dict:   # This line is used to replace the MOS parameters and simulation_parameters
            if "parameters "+param_name+'=' in line:
                line=line.replace(line,print_param(param_name,write_dict[param_name]))
        
        if write_check==1:
            s=s+line

    f.truncate(0)
    f.write(s)
    f.close()

#-----------------------------------------------------------------
# Function that modifies tcsh file
# Inputs  : circuit_initialization_parameters
# Outputs : NONE
def write_tcsh_file(circuit_initialization_parameters,optimiztion_type):
    
    filename=circuit_initialization_parameters['simulation']['standard_parameters']['tcsh']
    f=open(filename,'r+')
    s=''
    s='#tcsh\n'
    s=s+'source ~/.cshrc\n'
    
    if optimiztion_type=='basic':
        s=s+'cd '+circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'\n'
    """
    else:
        s=s+'cd '+circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']+circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']+'\n'
    """
    #s=s+'cp '+circuit_initialization_parameters['simulation']['standard_parameters']['run_directory']+'Netlists_Ref/circ.scs'+' ./\n'
    s=s+'spectre circ.scs =log circ_log.txt\n'
    s=s+'exit'
    
    f.truncate(0)
    f.write(s)
    f.close()

#-----------------------------------------------------------------
# Function that modifies tcsh file
# Inputs  : circuit_initialization_parameters
# Outputs : NONE
def copy_netlists(circuit_initialization_parameters,circuit_initialization_parameters_run):
    
    src=circuit_initialization_parameters['simulation']['standard_parameters']['run_directory']+'Netlists_Ref/circ.scs'
    for i in circuit_initialization_parameters_run.keys():
        #print(circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['sim_directory'],circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['basic_circuit'])
        dst_dir=circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['sim_directory']+circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['basic_circuit']
        s='cp '+src+' '+dst_dir+'/ '
        os.system(s)

#===========================================================================================================================


"""
====================================================================================================================================================================================
------------------------------------------------------------ SPECTRE RUNNING FUNCTIONS ---------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------    
# This function will run the shell commands to run Spectre
# Inputs  : Optimization Input Parameters
# Outputs : NONE
def run_file(circuit_initialization_parameters):
    os.system('cd /home/ee18b038/cadence_project')
    os.system('tcsh '+circuit_initialization_parameters['simulation']['standard_parameters']['tcsh'])   # This is the command to run the spectre file
    

#-----------------------------------------------------------------------------------------------
# This function will perform simulation for Basic Parameters
# Inputs  : Circuit_Parameters, circuit_initialization_parameters
# Outputs : Extracted_Parameters
def write_extract_basic(circuit_initialization_parameters):
    
    # Writing the tcsh file for Basic Analysis
    write_tcsh_file(circuit_initialization_parameters,'basic')
    
    # Writing the simulation parameters
    write_simulation_parameters(circuit_initialization_parameters)
    
    # Running netlist file
    run_file(circuit_initialization_parameters)
    
    # Extracting the Basic Parameters
    basic_extracted_parameters=extract_basic_parameters(circuit_initialization_parameters)
    
    return basic_extracted_parameters
"""
#-----------------------------------------------------------------------------------------------
# This function will perform simulation for Advanced Parameters
# Inputs  : Circuit_Parameters, circuit_initialization_parameters
# Outputs : Extracted_Parameters
def write_extract_advanced(circuit_initialization_parameters):

    pin_start=circuit_initialization_parameters['simulation']['standard_parameters']['pin_start']
    pin_stop=circuit_initialization_parameters['simulation']['standard_parameters']['pin_stop']
    pin_points=circuit_initialization_parameters['simulation']['standard_parameters']['pin_points']

    pin=np.linspace(pin_start,pin_stop,pin_points)
            
    pow_ph_params=np.zeros(pin_points,dtype=float)

    for i in range(pin_points): 
            
        circuit_initialization_parameters['simulation']['netlist_parameters']['pin']=pin[i]
                
        # Writing the simulation parameters
        write_simulation_parameters(circuit_initialization_parameters)

        # Writing the tcsh file for Basic Analysis
        write_tcsh_file(circuit_initialization_parameters,'basic')

        # Running netlist file
        run_file(circuit_initialization_parameters)

        # Extracting Vout Magnitude
        file_name=circuit_initialization_parameters['simulation']['standard_parameters']['sim_directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/circ.scs'
        pow_ph_params[i]=extract_dev_param(file_name,circuit_initialization_parameters)

    advanced_extracted_parameters=extract_advanced_parameters(pow_ph_params)
    
    return advanced_extracted_parameters
"""
#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
# Inputs  : Circuit_Parameters, circuit_initialization_parameters
# Outputs : Extracted_Parameters
def write_extract_single(i,circuit_parameters,circuit_initialization_parameters):
     
    # Writing to netlist file
    write_circuit_parameters(circuit_parameters,circuit_initialization_parameters)
    
    # Extracting the Basic Parameters
    basic_extracted_parameters=write_extract_basic(circuit_initialization_parameters)
    
    # Extracting the Advanced Parameters
    #advanced_extracted_parameters=write_extract_advanced(circuit_initialization_parameters)
    
    # Extracting Parameters from output files
    extracted_parameters=basic_extracted_parameters.copy()
    #for param_name in advanced_extracted_parameters:
    #    extracted_parameters[param_name]=advanced_extracted_parameters[param_name]

    return (i,extracted_parameters)

#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
# Inputs  : Circuit_Parameters, circuit_initialization_parameters
# Outputs : Extracted_Parameters
def get_final_extracted_parameters(extracted_parameters_combined):
    
    final_extracted_parameters={}

    extracted_parameters_select={
        'vg1':'mid',
        'vd1':'mid',
        'i_source':'mid',
        'v_source':'mid',
        'p_source':'mid',
        'Isup':'max',
        'Vsup':'max',
        'Psup':'max',
        'gm1':'mid',
        'gds1':'mid',
        'vth1':'mid',
        'vdsat1':'mid',
        'Idsat1':'max',
        'cgs1':'mid',
        'cgd1':'mid',
        'freq':'mid',
        'region':'mid',
        'am-pm-dev':'max',
        'ip1db_man':'min',
        'Voutdc':'mid',
        #'s12_db':'max',
        #'s21_db':'max',
        #'s22_db':'max',
        #'k':'min',
        #'nf_db':'max',
        #'iip3_dbm':'max'
    }

    for param in extracted_parameters_combined[0]:
        for i in range(3):
            final_extracted_parameters['comb_'+str(i)+'_'+param]=extracted_parameters_combined[i][param]
    
    for param in extracted_parameters_select:
        if extracted_parameters_select[param]=='mid':
            final_extracted_parameters[param]=extracted_parameters_combined[1][param]
        elif extracted_parameters_select[param]=='min':
            param_array=[]
            for i in extracted_parameters_combined:
                param_array.append(extracted_parameters_combined[i][param])
            final_extracted_parameters[param]=min(param_array)
        else:
            param_array=[]
            for i in extracted_parameters_combined:
                param_array.append(extracted_parameters_combined[i][param])
            final_extracted_parameters[param]=max(param_array)
    
    # Calculating the value of gain
    gain_array=[]
    gain_phase_array=[]
    for i in extracted_parameters_combined:
        gain_array.append(extracted_parameters_combined[i]['gain_db'])
        gain_phase_array.append(extracted_parameters_combined[i]['gain_phase'])
    gain_min=min(gain_array)
    gain_index=gain_array.index(gain_min)
    final_extracted_parameters['gain_db']=gain_min
    final_extracted_parameters['gain_phase']=extracted_parameters_combined[gain_index]['gain_phase']

    # Calculating the value of 1dB compression point
    ip1db_array=[]
    op1db_array=[]
    for i in extracted_parameters_combined:
        ip1db_array.append(extracted_parameters_combined[i]['ip1db_auto'])
        op1db_array.append(extracted_parameters_combined[i]['op1db_auto'])
    ip1db_min=min(ip1db_array)
    ip1db_index=ip1db_array.index(ip1db_min)
    final_extracted_parameters['ip1db_auto']=ip1db_min
    final_extracted_parameters['op1db_auto']=extracted_parameters_combined[ip1db_index]['op1db_auto']

    return final_extracted_parameters


#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
# Inputs  : Circuit_Parameters, circuit_initialization_parameters
# Outputs : Extracted_Parameters
def write_extract(circuit_parameters,circuit_initialization_parameters):
    
    pool=mp.Pool()

    # Creating new circuit parameter files
    circuit_parameters_run={}
    circuit_parameters_run[0]=circuit_parameters.copy()
    circuit_parameters_run[1]=circuit_parameters.copy()
    circuit_parameters_run[2]=circuit_parameters.copy()
    circuit_initialization_parameters_run={}
    

    # Getting the values of frequency and range
    f_operating=circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']
    f_range=circuit_initialization_parameters['simulation']['standard_parameters']['f_range']

    # Creating new circuit initialization parameters
    circuit_initialization_parameters_run[0]={}
    circuit_initialization_parameters_run[0]=copy.deepcopy(circuit_initialization_parameters)
    circuit_initialization_parameters_run[0]['simulation']['standard_parameters']['sim_directory']=circuit_initialization_parameters_run[0]['simulation']['standard_parameters']['sim_directory']+'T1/'
    circuit_initialization_parameters_run[0]['simulation']['standard_parameters']['tcsh']=circuit_initialization_parameters_run[0]['simulation']['standard_parameters']['tcsh']+'T1/spectre_run.tcsh'
    circuit_initialization_parameters_run[0]['simulation']['netlist_parameters']['fund_1']=f_operating-f_range

    circuit_initialization_parameters_run[1]={}
    circuit_initialization_parameters_run[1]=copy.deepcopy(circuit_initialization_parameters)
    circuit_initialization_parameters_run[1]['simulation']['standard_parameters']['sim_directory']=circuit_initialization_parameters_run[1]['simulation']['standard_parameters']['sim_directory']+'T2/'
    circuit_initialization_parameters_run[1]['simulation']['standard_parameters']['tcsh']=circuit_initialization_parameters_run[1]['simulation']['standard_parameters']['tcsh']+'T2/spectre_run.tcsh'
    circuit_initialization_parameters_run[1]['simulation']['netlist_parameters']['fund_1']=f_operating
   
    circuit_initialization_parameters_run[2]={}
    circuit_initialization_parameters_run[2]=copy.deepcopy(circuit_initialization_parameters)
    circuit_initialization_parameters_run[2]['simulation']['standard_parameters']['sim_directory']=circuit_initialization_parameters_run[2]['simulation']['standard_parameters']['sim_directory']+'T3/'
    circuit_initialization_parameters_run[2]['simulation']['standard_parameters']['tcsh']=circuit_initialization_parameters_run[2]['simulation']['standard_parameters']['tcsh']+'T3/spectre_run.tcsh'
    circuit_initialization_parameters_run[2]['simulation']['netlist_parameters']['fund_1']=f_operating+f_range

    #copying netlist files
    copy_netlists(circuit_initialization_parameters,circuit_initialization_parameters_run)
        
    # Creating processes
    results_async=[pool.apply_async(write_extract_single,args=(i,circuit_parameters_run[i],circuit_initialization_parameters_run[i])) for i in range(3)]
    
    extracted_parameters_combined={}
    for r in results_async:
        try:
            (i,extracted_parameters)=r.get()
        except mp.TimeoutError:
            print("timeout error")

        extracted_parameters_combined[i]=extracted_parameters
        
    final_extracted_parameters=get_final_extracted_parameters(extracted_parameters_combined)
    
    pool.close()
    pool.join()

    return final_extracted_parameters



#===========================================================================================================================

