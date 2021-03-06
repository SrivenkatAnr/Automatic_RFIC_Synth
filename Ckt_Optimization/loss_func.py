#===========================================================================================================================
"""
Name                : Srivenkat A
Roll Number         : EE18B038
File Name           : loss_func.py
File Description    : This file will contain the functions for Loss Optimization

Functions structure in this file:
    --> ramp_func
    --> calc_fom_1
    --> update_circuit_parameters
    --> calc_check_loss
    --> check_best_solution
    
"""
#===========================================================================================================================

from collections import OrderedDict
#===========================================================================================================================
#------------------------------------Defining the functions -----------------------------------------

#-----------------------------------------------------------------------------------------------
# This is the ramp function
# Inputs  : x
# Outputs : r(x)
def ramp_func(x):
    if x>0:
        return x
    else:
        return 0
    
#===========================================================================================================================
#--------------------------------------------Output Functions---------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
# This function calculates the loss for Io Optimization
# Inputs  : extracted_parameters,output_conditions,loss_weights
# Outputs : loss_dict
def calc_loss(extracted_parameters,output_conditions,loss_weights):
    
    # Extracted Values
    gain=extracted_parameters['gain_db']
    op1db=extracted_parameters['op1db_man']
    am_pm_dev=extracted_parameters['am-pm-dev']
    Io=extracted_parameters['Isup_hb']
    p_harm_ratio=extracted_parameters['p_harm_ratio']
    gain_phase=extracted_parameters['gain_phase']
    
    # Reference Values
    gain_ref=output_conditions['gain_db']
    op1db_ref=output_conditions['op1db']
    am_pm_dev_ref=output_conditions['am-pm-dev']
    p_harm_ratio_ref=output_conditions['p-harm-ratio']
    gain_phase_dev=output_conditions['gain-phase-dev']
    
    #Defining the weights to calculate Loss
    A1=loss_weights['gain_db']  # Weight for gain
    A2=loss_weights['op1db'] # Weight for output 1dB power
    A3=loss_weights['am-pm-dev']   # Weight for am-pm-deviation
    A4=loss_weights['Isup']   # Weight for Io
    A5=loss_weights['p-harm-ratio']   # Weight for Io
    A6=loss_weights['gain-phase'] #Weight for gain-phase
    
    # Calculating Loss
    loss_gain=A1*ramp_func(gain_ref-gain)
    loss_op1db=A2*ramp_func(op1db_ref-op1db)
    loss_am_pm_dev=A3*ramp_func(am_pm_dev-am_pm_dev_ref)
    loss_Io=A4*Io
    loss_p_harm=A5*ramp_func(p_harm_ratio_ref-p_harm_ratio)
    loss_gain_phase=A6*ramp_func(min(abs(gain_phase+180),abs(gain_phase-180))-gain_phase_dev)
    loss=loss_gain+loss_op1db+loss_am_pm_dev+loss_Io+loss_p_harm+loss_gain_phase
    #loss=loss_gain+loss_s11+loss_nf+loss_Io
    loss_dict={'loss':loss,'loss_gain':loss_gain,'loss_op1db':loss_op1db,'loss_am_pm_dev':loss_am_pm_dev,'loss_Io':loss_Io,'loss_p_harm':loss_p_harm,'loss_gain_phase':loss_gain_phase}
    #loss_dict={'loss':loss,'loss_gain':loss_gain,'loss_s11':loss_s11,'loss_nf':loss_nf,'loss_Io':loss_Io}
    
    return loss_dict

#-----------------------------------------------------------------------------------------------
# This function updates the values of circuit parameters by trying to minimize loss
# Inputs  : circuit_parameters,circuit_parameters_slope,check_loss,optimization_input_parameters
# Outputs : circuit_parameters
def update_circuit_parameters(cir,circuit_parameters_slope,check_loss,optimization_input_parameters,run_number):

    alpha_parameters=optimization_input_parameters['optimization'][run_number]['alpha']['values']

    # Calculating the value to update each parameter with
    for param_name in circuit_parameters_slope:
        
        # Calculating the Increment Value
        if check_loss==-1:
            change=circuit_parameters_slope[param_name]['loss']*(cir.circuit_parameters[param_name]**2)*alpha_parameters['common']*alpha_parameters[param_name]
        elif check_loss==1:
            change=circuit_parameters_slope[param_name]['loss_Io']*(cir.circuit_parameters[param_name]**2)*alpha_parameters['common']*alpha_parameters[param_name]
        else:
            change=(circuit_parameters_slope[param_name]['loss']-circuit_parameters_slope[param_name]['loss_Io'])
            change=change*(cir.circuit_parameters[param_name]**2)*alpha_parameters['common']*alpha_parameters[param_name]

        #print(change,param_name,'hi1')
    
    
        # Checking if the parameter is updated by a large value
        change_limit=0.25 # If the incremented value is more than +- change_limit*parameter_name, then we will limit the change
        if change>change_limit*cir.circuit_parameters[param_name]:
            change=change_limit*cir.circuit_parameters[param_name]
        elif change<-1*change_limit*cir.circuit_parameters[param_name]:
            change=-1*change_limit*cir.circuit_parameters[param_name]

        #print(change,param_name,'hi2')
        
        #if (param_name=='Rl') and ((cir.circuit_parameters[param_name]-change)>50):
        #    change=cir.circuit_parameters[param_name]-50  
        if (param_name=='Ld') and ((cir.circuit_parameters[param_name]-change)>10e-9):
            change=cir.circuit_parameters[param_name]-10e-9
        
        # Updating circuit_parameters
        cir.circuit_parameters[param_name]=cir.circuit_parameters[param_name]-change
        #print(cir.circuit_parameters[param_name],param_name)
        
    
#-----------------------------------------------------------------------------------------------
# This function will check the loss of gain, iip3, nf, and s11
# Inputs  : loss_iter,i,loss_type
# Outputs : check_loss ( -1 if s11 loss is 0 )
def calc_check_loss(loss_iter,i,loss_type):

    if loss_type==0:
        if loss_iter[i-1]['loss']==loss_iter[i-1]['loss_Io']:
            check_loss=1
        else:
            check_loss=0
                
    elif loss_type==1:
        check_loss=-1
        
    return check_loss
    
#---------------------------------------------------------------------------------------------------------------------------
# Function to check the best solution
# Inputs  : optimization_results,loss_max
# Outputs : opt_dict
def check_best_solution(optimization_results,loss_max):

    # Defining some values
    n_iter=optimization_results['n_iter']
    iter_min=0
    loss_Io_min=optimization_results['loss_iter'][0]['loss_Io']

    if (optimization_results['loss_iter'][0]['loss']-optimization_results['loss_iter'][0]['loss_Io'])>loss_max:
        flag=0
    else:
        flag=1
    
    for i in range(1,n_iter):
        if (optimization_results['loss_iter'][i]['loss']-optimization_results['loss_iter'][i]['loss_Io'])>loss_max:
            continue

        if flag==0 or (flag==1 and optimization_results['loss_iter'][i]['loss_Io']<loss_Io_min):
            iter_min=i
            loss_Io_min=optimization_results['loss_iter'][i]['loss_Io']
            flag=1

    # Creating output dictionary
    opt_dict=OrderedDict()
    opt_dict['loss_max']=loss_max
    if flag==1:
        opt_dict['perfect_point']='Yes'
    else:
        opt_dict['perfect_point']='No'
    opt_dict['iter_number']=iter_min+1
    opt_dict['Io_loss']=loss_Io_min
    
    return opt_dict
    

#===========================================================================================================================


        
    


