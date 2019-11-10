import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from astropy.io import fits

N = 20
m = 10
t = np.linspace(1,N,N)

x = np.zeros(N)

for i in range(N/2 -int(.2*N),N/2+int(.2*N)):
	x[i] = 1.0

y = np.zeros(N)
for i in range(0,m/2):
	y[i] = i
for i in range(m/2,m):
	y[i] = m - i

print y, np.size(y)
z = np.append(np.ones(m),np.zeros(N-m))

print z, np.size(z)
conv = np.zeros(N)
for i in range(0,N):
	for k in range(0,N):	
		if k-i>=0:
			conv[i]  = conv[i] + y[k] * x[i-k]
print np.size(conv)

test = signal.convolve(x,y,mode='same')

plt.figure(1)

plt.plot(y,'r',x,'b')
plt.figure(2)
plt.plot(t,conv)

plt.figure(3)
plt.plot(t,test)

print'enter the file name'
orig_data = raw_input("")


#with open(orig_data) as f:
#	lines = f.readlines()

data = np.loadtxt(orig_data)

print len(data)
plt.figure(4)
plt.plot(data[:,2],data[:,3])
plt.show()
