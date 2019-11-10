import subprocess
import numpy as np
import math
import psrchive
#First File=====================================================
cmd = "psrcat -db_file /home/kbansal/psrcat_tar/psrcat.db B2217+47"
p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
stdout,stderr = p.communicate()
period = 1/float(stdout.split('\n')[4].split()[3])
df1 = 9.536891E-15 # took off the negative sign


filename = raw_input('Enter First Filename:')
a = psrchive.Archive_load(filename)#reads in the raw data. 
data = a.get_data() #returns raw data as numpy array
profile = a.get_Profile(0,0,0) # This always requires four inputs: subintegrations, polarization, channel, and profile bin. 
d1 = data[0,0,0,:]
d = d1.ravel()#reshape the array
peak1 = np.argmax(d)

mjd1 = subprocess.Popen(["vap -c stt_imjd,stt_smjd,stt_offs %s"%filename],stdout=subprocess.PIPE,shell=True).communicate()[0].split()#[1]
print mjd1
start_time1 = float(mjd1[-3]) + (float(mjd1[-2])+float(mjd1[-1]))/(24.0*3600.0)

nbin = subprocess.Popen(["vap -c nbin %s"%filename],stdout=subprocess.PIPE,shell=True).communicate()[0].split()
#=================File2============================================

for i in range(0,10):
	filename = raw_input('Enter Second Filename:')
	a = psrchive.Archive_load(filename)#reads in the raw data. 
	data = a.get_data() #returns raw data as numpy array
	profile = a.get_Profile(0,0,0) # This always requires four inputs: subintegrations, polarization, channel, and profile bin. 
	d1 = data[0,0,0,:]
	d = d1.ravel()#reshape the array
	peak2 = np.argmax(d)

	mjd2 = subprocess.Popen(["vap -c stt_imjd,stt_smjd,stt_offs %s"%filename],stdout=subprocess.PIPE,shell=True).communicate()[0].split()#[1]
	print mjd2
	start_time2 = float(mjd2[-3]) + (float(mjd2[-2])+float(mjd2[-1]))/(24.0*3600.0)
	time_off = (start_time2 - start_time1)*24*3600.0

	phase_shift = (time_off/period - int(time_off/period))*1024#tperiod
	print time_off, time_off/period, int(time_off/period)
	print 1024 - phase_shift
	print peak1,peak2
	print period,int(time_off/period)
	i=i+1
	pulseTime = 0

'''
for i in range (0,(int(time_off/period)-1)):
	period=period + math.pow(round(period,10),3) * df1 * i
	pulseTime = pulseTime+period
print period
binoff =1024 - 1024* (time_off - pulseTime)/(period)#+period**2 * df1 * time_off)
print binoff,phase_shift

MJDPar = 46599.00
nPulse = (start_time1 - MJDPar)*24*3600.0/ period
nPulse2 = (start_time2 - MJDPar)*24*3600.0/ period
nbin = (nPulse - int(nPulse))*1024
nbin2 = (nPulse2 - int(nPulse2))*1024
print 'nbin from par MJD',nbin,nbin2
'''

