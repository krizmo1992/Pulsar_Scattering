from template_high import template
import re
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import sys
import subprocess #subprocess is a way to run a command form another program
import argparse

# This file takes high frequency profiles and generates an intrinsic pulse model
#===============================================
parser = argparse.ArgumentParser()
parser.add_argument('Usage: python freq_func.py High_frequency_profile_list PSRNAME',nargs = '+')

if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)

filename = sys.argv[1] #High frequency pulse file list
#print "input Low-Frequency filename (Average profile)"
psrname = filename.split('.')[0]
print 'PULSAR NAME', psrname

def func_ini(x,a,k): #test linear fit for initial values
	return a + (k*x)
def freq_func(x,alpha,kappa): #frequency dependence function
	return alpha*(x**kappa)

#===================================================================================

def freqFit(filename,psrname):
	flist = open(filename,'r')

	cmd = "psrcat -db_file /home/kbansal/psrcat_tar/psrcat.db -c \"p0\" -o short %s" % psrname
	p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	stdout,stderr = p.communicate()
	print stdout,stderr
	period = float(stdout.split('\n')[4].split()[1])
	print 'PULSAR PERIOD:',period

	freq = []
	compSep = []
	widthErr= []
	widthmat = []
	widthmat = []
	varList = []
	for line in flist:
		if '#' in line:#avoid comments
			pass
		else:
			if ".avg" in line:
				nu = float(re.findall("\d+\.\d+", line.split('_')[2])[0]) # identifying frequency
			else:
				nu = float(line.split('.')[0].split('_')[1]) # identifying frequency

			high_freq = line.strip()
			pars,j1,j2,nbin,var = template(high_freq,varList)
			varList.append(var)
			#print j1,j2,nu

			#Number of fitted components
			ncomp = np.size(pars)/3 #Three elements: mu,sig,amp
			widthErr = np.append(widthErr,1.0/nbin)
			freq= np.append(freq,nu)
			print freq
			if ncomp>1:
				width = []
				peakPos = []
				for i in range(0,ncomp):
					peakPos.append(pars[3*i+1])
					width.append(pars[3*i+2])
					print "Number of components:", ncomp
				print width
				widthmat.append(width)
				if ncomp==2:
					sep1 = peakPos[0]-peakPos[1]
					compSep.append(sep1)#*period)
				elif ncomp==3:
					sep1 = (peakPos[0]-peakPos[1])
					sep2 = (peakPos[1]-peakPos[2])
					sep3 = (peakPos[0]-peakPos[2])
					compSep.append([sep1,sep2,sep3])#*period)
			elif ncomp==1:
				widthmat.append(width)
			else:
				print "Pulse consists of zero components"
	flist.close()

	widthmat = np.asarray(widthmat)
	compSep = np.asarray(compSep)
	print 'width of main components:'
	print widthmat
	print 'Separation between components:'
	print compSep

#def nudepe():	
	if np.size(freq)>1:
		freqWid = [] # Width Frequency function
		freqSep = [] # Separation Frequency function
		x = np.log(freq)
		for j in range(0,ncomp):
			ywid = []
			ycompsep = [] # This is not equal to ncomp but ncompC2
			for i in range(0,np.size(freq)):
				ywid = np.append(ywid,widthmat[i][j])
				ycompsep = np.append(ycompsep,compSep[i][j])
			print 'width matrix :', np.size(ywid), ywid
			# Main fit
			pn,popt = fitting(freq,ywid,widthErr)
			if np.size(pn)<1:
				pass
			else:
				freqWid.append([j,pn])
				print 'Pulse width frequency parameters:',pn
				plt.figure(1)
				plt.title('Pulse width Evolution')
				plt.plot(freq,ywid,'r.',freq,freq_func(freq,pn[0],pn[1]),'b.')

			#Component separation fit
			print 'Separation between components:', compSep
			pn,popt = fitting(freq,ycompsep,widthErr)
			if np.size(pn)<1:
				pass
			else:
				freqSep.append([j,pn])
				print 'Pulse component separation frequency parameters:',pn
				plt.figure(2)
				plt.title('Component Separation Evolution')
				plt.plot(freq,ycompsep,'r.',freq,freq_func(freq,pn[0],pn[1]),'b.')
				plt.show()
	else:
		print pars
		freqWid =[]
		freqSep =[]
	return pars, freqWid, freqSep

def fitting(x,y,err):
	# Log test
	pn,popt = opt.curve_fit(func_ini,np.log(x),np.log(y))
	pn = [np.exp(pn[0]),pn[1]]
	return pn,popt	
#Return the new model for higher frequency which includes frequency dependence. 
#====================================================================================

pars,freqWid, freqSep = freqFit(filename,psrname)
filename = psrname+'.mod'
filename = open(filename,'w')
filename.write('%s\n' % '#Gauss_parmeters')
filename.write('%s\n' % str(pars))
filename.write('%s\n' % '#Width Frequency_parmeters')
filename.write('%s\n' % str(freqWid))
filename.write('%s\n' % '#Separation Frequency_parmeters')
filename.write('%s\n' % str(freqSep))
filename.close()

#subprocess is a way to run a command form another program
# scan psrcat to find the pulse period
