import numpy as np

def boxCar(seq,ns):
		if ns%2==0:
			ns = raw_input('Enter an odd number:')
		nsize = np.size(seq)
		tempSeq = np.zeros(nsize)
		nsum = ns/2+1
		for i in range(0,nsize):
				temp = seq[i]
				for j in range (1,nsum):
						if i+j < 20:    
								temp = temp + seq[i-j]+seq[i+j]
						#        print i-j,i+j
						else:
								#print (seq[i-nsum:20]), seq[0:i+nsum-19]
								temp = temp + seq[i-j]+seq[i+j-20]
			   
				tempSeq[i] = temp/float(ns)
		return tempSeq

if __name__ == "__main__":
        a = np.linspace(1,10,20)
        boxCar(a,3)
