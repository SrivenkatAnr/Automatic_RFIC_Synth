#===========================================================================================================================
"""
Name                : Srivenkat A
Roll Number         : EE18B038
File Description    : This file will perform the process analysis by changing the process and measuring the values
"""

#===========================================================================================================================
import numpy as np
import os
import common_functions as cf # type: ignore
from matplotlib import pylab
from pylab import *
#===========================================================================================================================


"""
============================================================================================================================
------------------------------------- Defining the functions ---------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of circuit_parameters to a txt file
def write_circuit_parameters(circuit_parameters,optimization_input_parameters):
    
    filename=optimization_input_parameters['filename']['output']    # Getting the filename
    
    newpath =filename+'/Process_Analysis/Results'   # Creating the folder if it is not present
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    filename=filename+'/Process_Analysis/Results/circuit_parameters.txt'
    
    f=open(filename,'w')
    for param_name in circuit_parameters:
        f.write(str(param_name)+'\t'+str(circuit_parameters[param_name])+'\n')  # Writing the values in the file
    f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the header row for extracted parameters to a csv file
def write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters,process_corner):
    
    newpath=optimization_input_parameters['filename']['output']+'/Process_Analysis/Results/'+str(process_corner)+'/'    # Creating the folder if it is not present
    if not os.path.exists(newpath):
        os.makedirs(newpath)
        
    filename=optimization_input_parameters['filename']['output']+'/Process_Analysis/Results/'+str(process_corner)+'/extracted_parameters.csv'   # Getting the filename

    f=open(filename,'w')
    f.write('Temperature,')
    for param_name in extracted_parameters:
        f.write(param_name+',') # Writing the names in the file
    f.write('\n')
    f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each process iteration to a csv file
def update_extracted_parameters(extracted_parameters,optimization_input_parameters,process_corner,temp):
    filename=optimization_input_parameters['filename']['output']+'/Process_Analysis/Results/'+str(process_corner)+'/extracted_parameters.csv'   # Getting the filename
    
    f=open(filename,'a')
    f.write(str(temp)+',')
    for param_name in extracted_parameters:
        f.write(str(extracted_parameters[param_name])+',')  # Writing the values to the file
    f.write('\n')
    f.close()
    

""" 
============================================================================================================================
------------------------------------- Output Functions ---------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function that will perform the temperature analysis
def process_analysis(cir,optimization_input_parameters,timing_results):
    
    if optimization_input_parameters['process_analysis']['run']=='NO':
        return

    # Opening the Run_Status File
    f=open(optimization_input_parameters['filename']['run_status'],'a')
    f.write('Process Analysis Start\n Time : '+str(datetime.datetime.now())+'\n\n')
    f.close()

    # Storing the starting time
    timing_results['process_analysis']={}
    timing_results['process_analysis']['start']=datetime.datetime.now()

    print('************************************************************************************************************')
    print('************************************* Process Analysis *****************************************************')
    
    cir.update_simulation_parameters(optimization_input_parameters['process_analysis']['simulation'])
    
    initial_extracted_parameters=cir.extracted_parameters.copy()
    initial_circuit_parameters=cir.circuit_parameters.copy()

    # Creating Dictionaries to Store Values
    extracted_parameters_iter={}
    
    # Writing the values to output files
    write_circuit_parameters(cir.circuit_parameters,optimization_input_parameters)
    
    # Creating the list of process corners
    process_corner_list=['tt','ff','ss']

    # Creating an array for temperature analysis
    start_temp=optimization_input_parameters['process_analysis']['start_temp']
    stop_temp=optimization_input_parameters['process_analysis']['stop_temp']
    n_temp=optimization_input_parameters['process_analysis']['n_temp']
    if n_temp==1:
        temp_array=np.array(['27'])
    else:
        temp_array=np.linspace(start_temp,stop_temp,n_temp)

    # Performing the analysis
    for process_corner in process_corner_list:
        extracted_parameters_iter[process_corner]={}
        cir.circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']=process_corner
        write_extracted_parameters_initial(cir.extracted_parameters,optimization_input_parameters,process_corner)
        cir.write_simulation_parameters()
        for temp in temp_array:
            cir.update_temp(temp)
            cir.run_circuit()
            update_extracted_parameters(cir.extracted_parameters,optimization_input_parameters,process_corner,temp)     # Writing the values to the output file
            extracted_parameters_iter[process_corner][temp]=cir.extracted_parameters.copy()
        
    # Restoring the value of initial extracted and circuit parameters
    cir.reset_temp()
    cir.update_circuit(initial_circuit_parameters.copy())
    
    # Plotting the graphs
    file_directory=optimization_input_parameters['filename']['output']
    plot_process_analysis(extracted_parameters_iter,file_directory)

    # Storing the stopping time
    timing_results['process_analysis']['stop']=datetime.datetime.now()

    # Opening the Run_Status File
    f=open(optimization_input_parameters['filename']['run_status'],'a')
    f.write('Process Analysis End\n Time : '+str(datetime.datetime.now())+'\n\n')
    f.close()
    

"""
============================================================================================================================
------------------------------------- Plotting Functions -------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# Plotting results ( Parameters vs Io at different temperatures )
def plot_process_analysis(extracted_parameters_iter,file_directory):
    
    # Creating a folder to store the results
    pathname=file_directory+'/Process_Analysis/Plots/'
    if not os.path.exists(pathname):
        os.makedirs(pathname)
    
    # Creating output conditions dictionary
    colour_dict_process={0:'red',1:'green',2:'blue',3:'cyan',4:'lime'}
    
    # Extracting the values
    extracted_matrix,process_array,temp_array,param_array=extract_process_analysis(extracted_parameters_iter)

    # Getting the lengths of the arrays
    n_process=len(process_array)
    n_temp=len(temp_array)
    n_param=len(param_array)

    # Plotting parameters vs temp for different process corners
    for k in range(n_param):
        figure()
        for i in range(n_process):
            plot(temp_array,extracted_matrix[i,:,k],color=colour_dict_process[i],label='Process : '+str(process_array[i]))
        xlabel('Temperature')
        ylabel(param_array[k])
        legend()
        grid()
        savefig(pathname+str(param_array[k])+'.pdf')
        close()

#-----------------------------------------------------------------------------------------------
# Function to extract the data from extracted_parameters_iter dictionary and store it in the form of a matrix
def extract_process_analysis(extracted_parameters_iter):
    
    # Assigning the array to store the loss variables
    process_array=[]
    temp_array=[]
    param_array=[]
    
    # Creating variables for the counting no of lines and variables
    n_process=0
    n_temp=0
    n_param=0
    
    # Storing the value of temp, current, and param arrays
    for process in extracted_parameters_iter:
        process_array.append(process)
        n_process+=1
        final_process=process

    for temp in extracted_parameters_iter[final_process]:
        temp_array.append(temp)
        n_temp+=1
        final_temp=temp

    for param in extracted_parameters_iter[final_process][final_temp]:
        param_array.append(param)
        n_param+=1
        
    # Creating array to store the values 
    extracted_matrix=np.zeros((n_process,n_temp,n_param),dtype=float)
    
    # Extracting the values of the variables
    i=0
    for process in extracted_parameters_iter:
        j=0
        for temp in extracted_parameters_iter[process]:
            k=0
            for param in extracted_parameters_iter[process][temp]:
                extracted_matrix[i,j,k]=extracted_parameters_iter[process][temp][param]
                k+=1
            j+=1
        i+=1
    
    return extracted_matrix,process_array,temp_array,param_array

#===========================================================================================================================