import numpy as np
import re
from scipy import signal
import psrchive
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy import stats
import pdb
from smooth import smooth
#returns template profile, off_pulse bin range


'off_pulse function takes in a pulse and returns peak-offset-pulse, with onpulse bounds'
def off_pulse(pulsef,nbin):
	peak = np.argmax(pulsef) #Peak position
	offset = nbin/2 - peak
	pulsef = np.roll(pulsef,offset) #move the peak to the center
	peak = np.argmax(pulsef) #Revised Peak position
	print 'Peak:', peak
	j1 = peak - int(0.2*nbin)
	j2 = peak + int(0.2*nbin)
	offpulse_mean = np.median(pulsef)
	#print 'off Pulse_mean:',offpulse_mean
	pulse = pulsef - offpulse_mean
	sigma = np.sqrt((np.mean(pulse[0:j1]**2) + np.mean(pulse[j2:nbin]**2))/2.0)
	print 'Sigma:', sigma
#	print 'SNR0,peak',max(pulse)/sigma,peak
	#Use signalTOnoise to identify onPulse bounds
	for i in range(0,peak):
		snr1 = pulse[i]/sigma
		#print 'SNR:', snr1,i
		if snr1>=3:
			j1 = i
			break
		i += 1
	for i in range(0,peak):
		j = nbin-i-1
		snr2 = pulse[j]/sigma
		if snr2>=3:
			j2 = j
			break
		i += 1
	return j1,j2,pulsef

def normalize(ydata):
	ydata = ydata - np.median(ydata)
	ymax = np.max(ydata)
	ydata = ydata/ymax
	return ydata


def f_test(chi1,chi2,dof1,dof2):
	#dof1>dof2
	F = ( (chi1-chi2)/(dof1-dof2) ) / (chi2/dof2)
	prob = stats.f.cdf(F,dof1-dof2,dof2)
	SF = stats.f.sf(1,dof1-dof2,dof2)
	ISF = stats.f.isf(0.000001,dof1-dof2,dof2)	
	print 'Inverse F-value for probability 1e-6:',ISF
	return F, 1-prob

def fwhm(pulse,nbin):
	'''
	width of pulse
	n1,n2 start and end regions of pulse
	'''
	n1 = 0
	n2 = np.size(pulse)
	peak = np.argmax(pulse)
	maxima = pulse[peak]
	i1=i2=0
	#print 'n1,n2,peak:',n1,n2,peak
	for i in range(n1,peak):
		t = pulse[i]
		if t>=maxima/2:
			i1 = i	
	for i in range(peak,n2):
		t = pulse[i]
		if t>=maxima/2:
			i2 = i
	return float(i2-i1)/(nbin)

def addGauss(x,*param):
	n= len(param)/3
	#offset = param[3*n]
	f = np.zeros(np.size(x))# + offset
	for i in range(0,n):
		amp = param[3*i]
		mu = param[3*i+1]
		sd = param[3*i+2]
		f = f+(amp * np.exp(-0.5*(((x-mu)/sd)**2)))#+amp2*np.exp(-((x-mu2)**2/2*(sig2**2))))
	return f

def addlorentz(x,*param):
	n= len(param)/3
	#offset = param[3*n]
	f = np.zeros(np.size(x))# + offset
	for i in range(0,n):
		amp = param[3*i]
		mu = param[3*i+1]
		sd = param[3*i+2]
		f = f+(amp *sd**2 /((x-mu)**2 + sd**2))#+amp2*np.exp(-((x-mu2)**2/2*(sig2**2))))
	return f

def gauss(x,amp,x0,sd):
	return amp * np.exp(-((x-x0)**2)/(2*(sd**2)))

def onpulse(pulse,nbin,x):
	peak = np.argmax(pulse) #Peak position
	offset = nbin/2 - peak
	pulse = np.roll(pulse,offset) #move the peak to the center
	peak = np.argmax(pulse) #Revised Peak position
	j1 = peak - int(0.25*nbin)
	j2 = peak + int(0.25*nbin)
	print j1,j2 #bound

	#Remove mean
	pulse = pulse - np.median(pulse)
	print 'Highest Peak:',np.max(pulse)

	#Sigma Calculation
	noise = np.sqrt((np.mean(pulse[0:j1]**2) + np.mean(pulse[j2:nbin]**2))/2.0)
	print 'Off-Pulse Sigma:',noise

	onpulse = pulse[0:nbin]
	xx = x[0:nbin]
	#print 'Variable List', varList

	#Considering only above 50% peak value
	limit = max(onpulse)*0.80
	#Revised time and pulse
	ty = onpulse[onpulse[:]>=limit]
	tx = x[onpulse[:]>=limit]
	#print 'Onpluse Matrix Size', ty.size

	#Initializing parameters before fitting
	mu = float(peak)/nbin
	initialPar = np.array([pulse[peak],mu,fwhm(pulse,nbin)]) #Highly depends on initial values!
	#print 'Initial Input Parameters for pulse fittin',initialPar
	print '**********************-----------------****************************'

	plt.figure(11)
	plt.title('50% pulse')
	plt.plot(tx,ty,'r')
	#first model fitting
	pt,popt = opt.curve_fit(addGauss,tx,ty,initialPar,maxfev=2000)

	#Testing the fit only for On_Pulse
	model = addGauss(xx,*pt)
	resid = onpulse - model

	#initial pars for f_test
	chi2 = ((resid**2).sum()/(noise**2))
	dof2 = np.size(xx) - np.size(pt)
	print 'Reduced chi-square value for fit1', chi2/dof2
	plt.figure(1)
	plt.plot(xx,model,'r')
	plt.figure(2)
	plt.subplot(211)
	plt.plot(xx,onpulse,'r.',xx,model,'b')#.mean(0),extent=(0,1,f_lo,f_hi))
	plt.subplot(212)
	plt.plot(xx,resid,'b.')#.mean(0),extent=(0,1,f_lo,f_hi))
	plt.show()
	Fval = 100 #for initializing the loop
	j=1
	#while (max(resid)/noise>=2.5 or chi2>1):
	colors = ["b.", "y.", "g.", "k.",'r.']

	while (Fval>15):
		ptFinal = pt

		#selecting points above 50-60% of residual peak, limit may vary with each pulsar
		limit = max(resid)*0.2
		'''
		print 'DATA Selected Bins:',x[onpulse[:]>=limit]
		print 'Residual Selected Bins:',x[resid[:]>=limit], x[np.argmax(resid)],np.argmax(resid)
		#find the jump between two peak regions		
		resBin = x[resid[:]>=limit]
		for i in range(1,np.size(resBin)):
			if int(1024*abs(resBin[i-1]-resBin[i]))>=10:
				print resBin[i-1],resBin[i]
				i1 = int(resBin[i-1]*1024)
				i2 = int(resBin[i]*1024)
		print i1,i2
		residP = np.argmax(resid)
		if abs(residP - i1)> abs(residP - i2):
			ty = onpulse[onpulse[i1:]>=limit]
			tx = x[onpulse[i1:]>=limit]
		else:
			print 'second loop'
			ty = onpulse[onpulse[:i2]>=limit]
			tx = x[onpulse[:i2]>=limit]	
		'''

		ty = onpulse[onpulse[:]>=limit]
		tx = x[onpulse[:]>=limit]	

		plt.figure(11)
		plt.title('50% pulse')
		plt.plot(tx,ty,'r')
		#print tx
		print '============Comparing Residuals:======='
		print 'Maximum(residual), 5% of pulse peak', max(resid), 0.05*max(pulse)
		print'Maximum:residual/noise', max(resid)/noise
		print ''
		
		# Comparing peaks of residual with main pulse
		if (max(resid))<(0.2*max(pulse)):
			print 'Break Reason: Residual Peak is smaller than 5% of main pulse peak'
			break
		# Comparing peaks of residual rms from low frequency with pulse residual
		'''
		if len(varList)>0:
			print 'Residual Noise at 79.2 Mhz:', varList[0]
			print 'Residual Noise at higher frequency:', (np.var(resid))**0.5
			if (varList[0]/((np.var(resid))**0.5))>2:
				print 'Break Reason: Residual have noise of the order 79.2 MHz'
				break
		'''	
		chiold = chi2
		dof1 = dof2
		
		#re-initialize parameters	
		mu = np.argmax(resid)/float(nbin)
		maxima = max(resid)
		wd = fwhm(resid,nbin)
		print 'Old parameter for add Gauss:',pt
		initialPar = np.append(pt,[maxima,mu,wd])
		#print 'Initial Parameter for add Gauss:',initialPar
		print ''

		#fitting
		try:
			pt,popt = opt.curve_fit(addGauss,tx,ty,initialPar,noise)
		except RuntimeError:
			print 'Break Reason: There is runtimeError'
			break
		print 'Parameters from add Gauss:',pt
		model = addGauss(xx,*pt)
		resid = onpulse - model
		print 'fwhm(residual):',fwhm(pulse,nbin),'peak:',peak

		print '=========Chi-square tests============='
		chi2 = ((resid**2).sum()/(noise**2))
		dof2 = np.size(xx) - np.size(pt)
		print 'Reduced chi-square value', chi2/dof2
		print chiold,chi2,dof1,dof2

		#comparison
		test1 = f_test(chiold,chi2,dof1,dof2)
		prob = test1[1]
		Fval = test1[0]
		print 'F-test parameters:F-value,p-value',test1
		j+=1 #counter
		print 'Number of components',j

		plt.figure(1)
		plt.plot(xx,model,'r')
		for i in range(0,j):	
			plt.plot(xx,gauss(xx,pt[3*i],pt[3*i+1],pt[3*i+2]),colors[i])	
		plt.figure(2)
		plt.subplot(211)
		plt.plot(xx,onpulse,'r.',xx,model,'b')#.mean(0),extent=(0,1,f_lo,f_hi))
		plt.subplot(212)
		plt.plot(xx,resid,'b.')#.mean(0),extent=(0,1,f_lo,f_hi))
		plt.show()
		
		#Break points:
		if (chi2/dof2)<1: #if reduced chi2 is less than 1
			ptFinal = pt
			print 'Break Reason: reduced chi2 is less than 1'
			break	

		if (chi2/dof2)>(chiold/dof1): # If new reduced chi2 is larger than previous, stop!
			print 'Break Reason: reduced chi2 is larger than previous'
			break
	
		if (Fval)<15: # if Fval reaches below 15 stop and assign previous pt value as final
			print 'Break Reason: Fval reaches below 15 '
			break

		if (max(resid)/noise<=3):
			ptFinal = pt
			print 'Break Reason: SNR for Residual<=3'
			break

		if j>=5:
			print 'Break Reason: Number of components is more than 5'
			break
	
	# Parameters for width and separations
	n = np.size(ptFinal)/3
	pars = []
	for i in range(0,n):			
		pars = np.append(pars,[ptFinal[3*i+1],ptFinal[3*i+2]])
	return ptFinal#,(np.mean(resid**2))**0.5

def profilefit(pulse,x,par,peak):
	# probably needs normalization with peak
	# Assuming both pulse and par are for normalized profiles
	nbin = x.size
	j1,j2,pulse = off_pulse(pulse,nbin)
	pulseMax = np.max(pulse)
	'''	
	offset = nbin/2 - np.argmax(pulse)
	pulse = np.roll(pulse,offset) #move the peak to the center
	print 'before fitting',par
	'''	
	xx = x#[j1:j2]	
	pulse = normalize(pulse)#[j1:j2]	
	y=pulse
	#Normalization
	ngauss = np.size(par)/3
	for i in range(0,ngauss):
		par[3*i] = pulseMax*par[3*i] / peak
	#pn,popt = opt.leastsq(addGauss,x,pulse,par)

	pn,popt = opt.curve_fit(addGauss,xx,y,par)
	print 'After fitting',pn
	resid = pulse - addGauss(x,*pn)
	#chi-square
	noise = np.median((pulse - np.median(pulse))**2) **0.5
	chi2 = (((pulse-addGauss(x,*pn))**2).sum()/(noise**2))
	dof2 = np.size(x) - np.size(pn)
	print 'postfit Reduced Chi2', chi2/dof2
	#Plotting Figures
	colors = ["b.", "y.", "g.", "k.",'r.']
	plt.figure(1)
	for i in range(0,ngauss):	
		plt.plot(x,gauss(x,pn[3*i],pn[3*i+1],pn[3*i+2]),colors[i])	
	plt.figure(2)
	plt.subplot(211)
	plt.plot(x,addGauss(x,*pn),'r',x,pulse,'b.')
	plt.subplot(212)
	plt.plot(x,resid,'b.')
	plt.show()
	return pn
	
def template(filename):
#	print 'Enter the name of high frequency profile file (a text file)'
	print filename
	if ".txt" in filename:
		data_high = np.loadtxt(filename.strip())
		pulse_h = data_high[:,3] #the third column contains intensities along the time axis 
		#print pulse_h
	elif ".avg" in filename:
		print filename
		a = psrchive.Archive_load(filename.strip())#reads in the raw data
		d1 = a.get_data()[0,0,0,:] #returns raw data as numpy array
		pulse_h = d1.ravel()
	elif '.asc' in filename:
		data = open(filename.strip(),'r')
		pulse_h = []
		for line in data:
		    t = line.strip()
		    y = t.split()
			#print len(y)
		    if y[1]=='(s)' or y[1]=='Data...':
				pass
		    elif y[1]=='Stokes' or y[1]=='I' or y[1]=='Q':
				pass
		    elif y[1]=='Power' or y[1]=='only':
				pass
		    elif y[1]=='Time':
				pass
		    else:
				pulse_h = np.append(pulse_h,float(y[1]))
		data.close()
	else:
		pulse_h = filename
	nbin = pulse_h.size
	# Defining the x axis
	xAxis  = np.linspace(0,1,nbin)
	return pulse_h, xAxis,nbin

if __name__ == '__main__':
	filename = raw_input('Enter LOFAR frequency file: ')
	flist  = open(filename,'r')
	pulse,x,nbin = template(filename)	
	parLofar = onpulse(pulse,nbin,x)
	print parLofar
	#parLofar = [0.90504737, 0.50092745, 0.00355561, 0.11579287, 0.53353479, 0.00688251, 0.13855569, 0.49343699, 0.0096567,  0.0580849, 0.4545365, 0.00434808] #B0329+54
	#parLofar = [0.81053122, 0.50096634, 0.00673544, 0.19579593, 0.49808568, 0.01625497, 0.03744238, 0.58914707, 0.010614]#B0823+26
	#parLofar = [0.61751921,0.49924811,0.00386558,0.21667341,0.48038841,0.01465707,0.36833315,0.506996,0.01199936]
	
	lofarPeak = max(parLofar)
	fname = raw_input('Enter FileLIst name: ')
	psrname = fname.split('.')[0]
	print psrname
	parFile = psrname + '.par'
	flist  = open(fname,'r')
	fpar = np.empty([np.size(parLofar)+1,0])
	print fpar
	for line in flist:
		if '#' in line:#avoid comments
			pass
		else:
			if ".avg" in line:
				nu = float(re.findall("\d+\.\d+", line.split('_')[2])[0]) # identifying frequency
			else:
				nu = float(line.split('.')[0].split('_')[1]) # identifying frequency

			pulse,x,nbin = template(line.strip())
			pars = profilefit(pulse,x,parLofar,lofarPeak)
			pars = np.append(pars,nu)
			fpar = np.column_stack((fpar,pars))
	np.savetxt(parFile,fpar)
	#parFile.close()
			#Save them to a file
	print 'testing this programming procedure'	
