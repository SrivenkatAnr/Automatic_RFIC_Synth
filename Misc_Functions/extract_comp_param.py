import fileinput
import numpy as np
import matplotlib.pyplot as plt

fname="circ.raw/swp-{}_phdev_test.fd.pss_hb"

pin_start=-60
pin_stop=20
npin=pin_stop-pin_start+1

rvin_arr = np.zeros(npin,dtype=float)
ivin_arr = np.zeros(npin,dtype=float)
rvout_arr = np.zeros(npin,dtype=float)
ivout_arr = np.zeros(npin,dtype=float)
ph_arr = np.zeros(npin,dtype=float)
gdb_arr = np.zeros(npin,dtype=float)

def extract_file(file_name):
    f=open(file_name)
    lines=f.readlines()
    f.close()
    return lines

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

def extract_volt(lines):
    lines=lines.split()
    char_r=lines[1].split('(')[1]
    char_i=lines[2].split(')')[0]

	# Converting string to floating point value
    vout_r=valueE_to_value(char_r)
    vout_i=valueE_to_value(char_i)
    return (vout_r, vout_i)

def calculate_phase_old(vout_re,vout_im,vin_re,vin_im):
	
	# Calculating the phase of vout and vin
	if vout_re>0:
		vout_phase=np.arctan(vout_im/vout_re)*180/np.pi
	else:
		vout_phase=180+np.arctan(vout_im/vout_re)*180/np.pi
	if vin_re>0:
		vin_phase=np.arctan(vin_im/vin_re)*180/np.pi
	else:
		vin_phase=180+np.arctan(vin_im/vin_re)*180/np.pi
	
	# Calculating the phase of the gain
	phase=vout_phase-vin_phase
	while phase<-180:
		phase+=180
	while phase>180:
		phase-=180

	return phase

def calculate_phase(vout_re,vout_im,vin_re,vin_im):
    
    vout=np.complex(vout_re,vout_im)
    vin=np.complex(vin_re,vin_im)
    phase=np.angle(vout/vin)*180/np.pi
    while phase<-180:
        phase+=360
    while phase>180:
        phase-=360

    return phase

def calculate_gain(vout_re,vout_im,vin_re,vin_im):
    gain=(vout_re**2+vout_im**2)/(vin_re**2+vin_im**2)
    gain_db=10*np.log10(gain)

    return gain_db

freq_flag=0
op_freq=1.9e9
for i in range(npin):
    if (i==0):
        filename=fname.replace("{}","000")  
    elif (i<=9):
        filename=str(fname.format("00"+str(i)))
    elif (i<=99):
        filename=str(fname.format("0"+str(i)))
    else:
        filename=str(fname.format(i)) 
    for line in extract_file(filename):
        if ('"freq"' in line) and ('"sweep"' not in line):
            freq=line.split()[1]
            freq=valueE_to_value(freq)
            if(freq==op_freq):
                freq_flag=1
        if ('"Vin"' in line) and ('"V"' not in line) and freq_flag==1:
            rvin_arr[i],ivin_arr[i] = extract_volt(line)
        if ('"Vout"' in line) and ('"V"' not in line) and freq_flag==1:
            rvout_arr[i],ivout_arr[i] = extract_volt(line)
            freq_flag=0
    ph_arr[i] = calculate_phase(rvout_arr[i], ivout_arr[i], rvin_arr[i], ivin_arr[i]) 
    gdb_arr[i] = calculate_gain(rvout_arr[i], ivout_arr[i], rvin_arr[i], ivin_arr[i])

plt.plot(np.linspace(pin_start,pin_stop,npin),ph_arr)
plt.title("ph_dev")
plt.show()

plt.plot(np.linspace(pin_start,pin_stop,npin),gdb_arr)
plt.title("gain_comp")
plt.show()



