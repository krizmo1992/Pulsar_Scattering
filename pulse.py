#from frequency_func import freqFit
from template_high import addGauss
from template_high import off_pulse
from smooth import smooth
from boxCar import boxCar
import numpy as np
import psrchive
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy.optimize as opt
from scipy import stats
import subprocess
import operator
import sys
from pulseLib import template,chi2,filter
 #=======================Lower Frequencies=========================#
# This program needs two average profile files. The first input file is at lower frequency (lwa DATA) and the second profile is at a higher frequency or expected model for a lower frequency. The aim of this program is to obtain scattering coefficient of the interstellar medium using the pulsar data. --kb,July,2017

SCREEN_WIDTH = 80
centered = operator.methodcaller('center', SCREEN_WIDTH)
print(centered("#########START Pulse FITTING###########"))
print "'''''''''''''''''''''''''''''''''''''''''''''"
print "Input Low-Frequency filename (Average profile)" 

# Read in the archive files
orig_data = sys.argv[1] #takes in data file_name
print 'File Name:',orig_data
#Pulsar Name from file
psrname = orig_data.split('_')[1]
#MJD from file
MJD = orig_data.split('_')[0]

#Pulsar period from psrcat command
cmd = "psrcat -db_file /home/kbansal/psrcat_tar/psrcat.db %s" % psrname
p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
stdout,stderr = p.communicate()
print stdout,stderr
period = 1/float(stdout.split('\n')[4].split()[3])
DM = float(stdout.split('\n')[4].split()[6])
print "'''''''''''''''''''''''''''''''''''''''''''''"
print 'Pulsar Period:', period
print 'Pulsar DM:', DM
print "'''''''''''''''''''''''''''''''''''''''''''''"

# estimate the average of pulse such that it excludes the on-pulse region lower and upper bin range to separate the on-pulse region from on-pulse region
# need to round off the average such that it matches the data points. To find the location of FWHM, some range will be useful as it may be difficult to find exact number. 

#===============================================================================================================
'''
This file should be averaged over all the subintegrations and contains only 4 channels. We are not using this since it will provide very poor SNR.
'''
orig_data = open(orig_data,'r')
freqMat = []
dMat = np.empty((0,256))
#dMat = np.empty((0,512))
for files in orig_data:
	if '#' in files:
		pass
	else:
		nu = float(files.strip().split('_')[2][0:4])
		print 'Current Frequency:',nu
		a = psrchive.Archive_load(files.strip())#reads in the raw data. 
		data = a.get_data() #returns raw data as numpy array
		#print data.shape #array dimensions
		profile = a.get_Profile(0,0,0) # This always requires four inputs: subintegrations, polarization, channel, and profile bin. 
		#print 'Size of Profile:', np.size(profile), np.shape(profile)
		#nbin = profile.get_nbin() #Maximum bin number in the average profile
		d = data[0,0,:,:]
		print 'Size and Shape of Data Array:',np.size(d),np.shape(d)
		nchan,nbin = np.shape(d)
		bins  = np.linspace(0,1,nbin)#bin array
		if nchan==2:
			nu1 = nu - 4.9
			nu2 = nu + 4.9
			freqMat.extend([nu1,nu2])
		elif nchan==4: # Using only top three frequencies
			nu1 = nu - 9.8 + 2.45
			nu2 = nu - 4.9 + 2.45
			nu3 = nu + 2.45
			freqMat.extend([nu1,nu2,nu3])
		else:
			freqMat.extend([nu])
		dMat = np.vstack((dMat,d))	
print 'Frequencies',freqMat
print 'Shape of data matrix',np.shape(dMat)

#========================== Reading in the Pulse Model Parameters===============================================
filename = psrname + '.mod'
f = open(filename,'r')
par = []# Pulse Components 
for line in f:
	if '#' in line:
		pass
	else:
		par.append(line)
		#print par
		#print map(float,pars.split())
f.close()

# Gaussian components
pars = map(float,par[0].strip().strip('[]').split(','))

# component-width frequency dependence
pnWidth = map(float,par[1].replace('[','').replace(']','').split(','))
print 'Model Width Parameters',pnWidth, 'Model With par size',np.size(pnWidth)

# component-Amplitude Ratio frequency dependence
pnAmpr = map(float,par[2].replace('[','').replace(']','').split(','))
print 'Model Amplitude Ratio Parameters', pnAmpr

ncomp = np.size(pnWidth)/3 # Number of guassian componenets in model

#pnSep = map(float,par[2].replace('[','').replace(']','').split(','))

#================================Main Program==========================
plt.figure(12)
plt.figure(figsize=(7,9)) # height and width of plot in inches
plt.ylabel('Amplitude',fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(np.arange(-0.1,1.0,0.5),fontsize=8)
plt.title('Pulse fitting at %s MHz' %(nu))
tauMat = []
tauErr = []

x = bins # X-axis Variable
count=-1
for nu in freqMat:
	print ''
	print '----------------------------------------------------------------------------------'
	print '===========================Frequency:',nu,'Fitting==============================='		

	#Calculating the model
	if ncomp==1:
		alpha_wd = pnWidth[1]
		kappa_wd = pnWidth[2]
		pars[2] = alpha_wd*(nu**(kappa_wd))# Width of pulse
	elif ncomp>1:
		#pars[1] = pars[4] + alpha_sep*nu**(kappa_sep)# + c_sep #separation of pulse 
		for j in range(0,ncomp):
			alpha_wd = pnWidth[j*3+1]
			kappa_wd = pnWidth[j*3+2]
			pars[3*j+2] = alpha_wd*	nu**(kappa_wd) #+ c_wd #Width of pulse 
		if ncomp==2:
			'''
			alpha_sep = pnSep[1]
			kappa_sep = pnSep[2]
			pars[4] = pars[1] +alpha_sep*nu**(kappa_sep) #Width of pulse 
			'''
			alpha_ampr = pnAmpr[1]
			kappa_ampr = pnAmpr[2]
			pars[3] = pars[0] /(alpha_ampr*(nu**(kappa_ampr))) #Amplitude of pulse
		if ncomp>=3:
			for j in range(0,ncomp-1):
				'''
				alpha_sep = pnSep[j*3+1]
				kappa_sep = pnSep[j*3+2]
				pars[3*(j+1)+1] = pars[1] -alpha_sep*nu**(kappa_sep)
				print 3*j+2, pars,nu,alpha_wd,kappa_wd
				'''
				alpha_ampr = pnAmpr[j*3+1]
				kappa_ampr = pnAmpr[j*3+2]
				pars[3*(j+1)] = pars[0]/(alpha_ampr*nu**(kappa_ampr))
				print pars[0],alpha_ampr,nu,kappa_ampr,alpha_ampr*nu**(kappa_ampr) 
	else:
		print 'No GAUSSIAN Components are available.'
	count = count+1
	model = addGauss(x,*pars)
	model = model/sum(model)
	print 'Model Parameters:',pars

#============================Simulating Data======================================================
	'''
	noise = np.random.normal(0,1,nbin)
	d = np.convolve(model,np.exp(-5*(count+1)*(x)),'same') /sum(np.exp(-5*(count+1)*(x))) #+ 0.1*noise
	d = d/max(d) + 0.1*noise
	'''
#====================================================================================
	d = dMat[count,:]
	
	nbin = np.size(d)
	print 'size(d):',np.size(d)
	'''
	plt.figure(20)
	plt.plot(x,d,'-')
	plt.figure(21)
	plt.plot(x,model,'-')
	plt.show()
	'''
	# shifting the peak to center
	peak = np.argmax(smooth(d,5))
	d = np.roll(d, (int(nbin*0.1) -peak))
	peak = np.argmax(smooth(d,5))

	#ON_Off Pulse regions
	#j3 = peak - int(0.25*nbin)#/((count+1)**0.5))
	j2 = peak + int(0.25*nbin)#/((count+1)**0.5))
	j3 = 1
	print 'Bin Selection Range',j2,j3

	# REMOVING Pedestal
	mean0 = np.median(d)
	d = d - mean0
	d = d/np.max(d)
	#==================Off and On Pulse region=============================='
	offpulse2 = d[j2:(nbin-1)]
	offpulse1 = d[0:j3]
	offpulse = np.concatenate((offpulse1,offpulse2),axis=0)

	sigma0 = np.average((offpulse)**2)**0.5
	sigma0 = np.std(offpulse)
	print '===================Data Parameters=========================='
	print 'peak:',peak,'SNR of Pulse', np.max(d)/sigma0
	print 'Pedestal:', mean0,'Sigma:',sigma0#, np.size(yy)

	'''
	plt.figure(3)
	plt.plot(xx,modelCut)
	plt.figure(4)
	plt.plot(xx,yy)
	plt.show()
	'''

#=====================================================================================================		
	#Fourier tansform of Signal and probable exponential

	#r1 = 0.05;r2=0.20#B0823+26
	r1 = 0.15;r2=0.4#B1842+14
	#r1 = 0.15;r2=0.30#B1822-09
	#r1 = 0.1;r2=0.15#B1842+14
	#r1 = 0.05;r2=0.1

	#plt.figure(11)
	#plt.plot(xx,yy)

	#modelCut = model[j3:j2] #Selected model
	#xx = np.linspace(0,1,len(modelCut)) # Selected x-axis
	#j3 = peak - int(0.1*nbin)#/((count+1)**0.5))
	#j2 = peak + int(0.4*nbin)#/((count+1)**0.5))
	j2 = peak + int(0.5*nbin)#/((count+1)**0.5))
	j3 = 1
	yy = d[j3:j2] #Selected data
	xx = x[j3:j2]
	dd = filter(d,x,r1,r2)

	scale_f = max(smooth(d,5))/max(model)
	print '============Fitting Initial Pars==============='
	print 'peak_value_Scale:',scale_f
	print 'Error bar on points:',sigma0
	print '================================================'

	boundRange = ([0,1./256], [np.inf,1.0])
	ytest = np.ones(np.size(yy))
	#Baseline offset has already been taken care of.We fit for the amplitude scaling and tau value
	#pt,pot = opt.curve_fit(template,(xx,modelCut,dd),yy,initialPar,sigma0,maxfev=2000)
	lowestChi2 = 10000
	p1 = 1000 
	for k in range(0,3):
		initialPar = np.array([scale_f, 10.0**(-k)])#/(count+1.0)])#,mean0])
		#print 'Intial test Parameters',initialPar
		#ptTest,pot = opt.curve_fit(template,(x,model,dd,j3,j2),yy,sigma=sigma0*ytest,p0=initialPar,bounds = boundRange,maxfev=2000)
		ptTest,pot = opt.curve_fit(template,(x,model,dd,j3,j2),yy,sigma=sigma0*ytest,p0=initialPar,maxfev=2000)
		print 'Scattering Fitting Function Values Set1:',ptTest#,pot
		print 'Reduced Chi-square Value:', chi2(ptTest,(x,model,dd,j3,j2),yy,sigma0)
		if lowestChi2>=chi2(ptTest,(x,model,dd,j3,j2),yy,sigma0):
			lowestChi2 = chi2(ptTest,(x,model,dd,j3,j2),yy,sigma0)
			pt = ptTest
		'''
		if p1>=(ptTest[0]):
			p1 = ptTest[0]
			pt = ptTest
		'''
	#Revised Initial Parameters
	initialPar = np.array([pt[0],pt[1]])#,mean0])

	ymin = min(yy)
	ymax = max(yy)
	'''
	plt.figure(2)
	plt.subplot(211)
	plt.ylabel('Amplitude',fontsize=12)
	plt.title('Pulse fitting at %s MHz' %(nu))
	plt.ylim([ymin,ymax])
	plt.plot(xx,yy,'r.',label='Observed Profile')
	plt.plot(xx,template((x,model,dd,j3,j2),*pt),'b',label='Template')#.mean(0),extent=(0,1,f_lo,f_hi))
	'''

	#r1 = 0.1;r2=0.35
	#dd = filter(d,x,r1,r2)
	ytest = np.ones(np.size(d))
	pt,pot = opt.curve_fit(template,(x,model,dd,0,nbin),d,initialPar,sigma0*ytest,maxfev=2000)
	print 'Scattering Fitting Function Values Set2:',pt,pot

	#pt1 = opt.fmin(residual,initialPar,args=(x,d))
	
	#pt,pot = opt.least_squares(template,initialPar,sigma0,maxfev=2000)
	#=================Scattering Matrix===================================		
	tauc = abs(pt[1])
	print np.diag(pot)
	tauErr = np.append(tauErr,np.sqrt(np.diag(pot))[1])
	tauMat = np.append(tauMat,tauc)
	
	'''
	plt.subplot(212)
	#plt.subplot(np.size(freqMat),1,count+1)
	plt.plot(x,d,'r',label='Observed Profile at %s MHz' %(nu))
	plt.plot(x,template((x,model,dd,0,nbin),*pt),'b',label='Template')
	plt.legend(loc=2, prop={'size': 1})
	plt.show()	

	'''
	plt.subplot(np.size(freqMat),1,count+1)
	plt.plot(x,d,'r-',label='Observed Profile at %s MHz' %(nu))
	plt.plot(x,template((x,model,dd,0,nbin),*pt),'b',label='Template')
	plt.legend(loc=2, prop={'size': 1})
	#plt.show()	
	#y-axis limits
	ymin = min(d)
	ymax = max(d)
	'''
	plt.figure(2)
	plt.subplot(211)
	plt.ylabel('Amplitude',fontsize=12)
	plt.title('Pulse fitting at %s MHz' %(nu))
	plt.ylim([ymin,ymax])
	plt.plot(x,d,'r',label='Observed Profile')
	plt.plot(x,template((x,model,dd),*pt),'b',label='Template')#.mean(0),extent=(0,1,f_lo,f_hi))
	plt.legend(loc=2)

	plt.subplot(212)
	resid = d - template((x,model,dd),*pt)
	chi2 = (resid**2).sum()/((nbin-np.size(pt))*sigma0**2)
	print 'chi-square value', chi2
	plt.title('Residual')
	plt.xlabel('Pulse Phase',fontsize=12)
	plt.ylim([ymin,ymax])
	plt.plot(x,resid,label='Residual')
	#plt.show()
	'''
	print '-----------------------------------xxxxxxx-----------------------------------------------'

print'============================================'
print 'TAU_MAT:',tauMat
print 'TAU_MATERR:',tauErr
print'============================================'
a1 = 1;a2 = len(tauMat)#-2
#tauErr = tauErr[a1:a2]*(1.0/(tauMat[a1:a2]**2)) * period#*1.0e3
tauErr = tauErr[a1:a2] * period#*1.0e3
tauMat = tauMat[a1:a2] * period #*1.0e3
plt.savefig('Frequncy_fit.pdf')
logF = np.log(freqMat[a1:a2])
tauErrl = tauErr/tauMat
logTau = np.log(tauMat)
def fitt(x,*par):
	alpha = par[0]
	b = par[1]
	return b*(x**alpha)
def fittlog(x,*par):
	alpha = par[0]
	b = par[1]
	return b + alpha * x
pi = [.1,np.log(DM**2)]
par,pcov = opt.curve_fit(fittlog,logF,logTau,pi,tauErrl)

par = [par[0],np.exp(par[1])]
print "Tau spectral Index:", par
tauIndex = par[0]
tauEindex = np.sqrt(np.diag(pcov))[0]
print "Spectral Index error",tauEindex
#print "Spectral Index error",((np.diag(pcov))**0.5)[0]
plt.figure(5)
plt.title('Frequency dependence of Scattering time')
plt.xlabel('Frequency in MHz')
plt.ylabel('Tau in seconds')
plt.errorbar(freqMat[a1:a2],tauMat,yerr=tauErr,label = 'Tau values for each frequency')
x = np.linspace(freqMat[a1],freqMat[np.size(freqMat[a1:a2])-1],20)
plt.plot(x,fitt(x,*par),'r-',label = 'Fitted Model')
plt.legend(loc=1)
#plt.show()
plt.savefig('TauFitplt')
filename = MJD +'_'+psrname+'.alpha'
np.savetxt(filename,(tauIndex,tauEindex))
filename = MJD +'_'+psrname+'.tau'
np.savetxt(filename,(tauMat,tauErr,freqMat[a1:a2]))
