#This code uses a 4th order Runge-Kutta method (see Iserles pg.41) to solve
#Burgers equation: 
#u_t + b uu_x = 0
################################################################################
import numpy as np
import scipy as sp
import scipy.fftpack as spf
################################################################################
#In terms of the global variables in main we'll usually have:
#uf = ufourier
#M = Mtime
#N = 2*Nfourier +2*Kpad + 1
#X = Xspatial
#b = b
#A = Adap_param
#Tf = Tfin

#This is the function for the Euler / RK4 method for the gCLM Equation
#n = current fourier coefficient
def gCLMfourier(uf,n,N,b):

	z=0
	l = n
	l = float(l)
	if n > 0:
		for k in range(n-N,N+1): 
			z = z - (1/l)*(1j*(n-k)*abs(n-k)*uf[k]*uf[n-k]) \
			- (1/l)*b*(1j*(k)*abs(n-k)*(uf[k]*uf[n-k]))
			# z = z - (1j*k*abs(k)*uf[k]*uf[n-k]) \
			# - b*(1j*(n-k)*abs(k)*(uf[k]*uf[n-k]))
			#z = z - (1j*k*uf[k]*uf[n-k]) - b*(1j*k*(uf[k]*uf[n-k]))
	elif n < 0:
		for k in range(-N,N+n+1):
			z = z - (1/(-l))*(1j*(n-k)*abs(n-k)*uf[k]*uf[n-k]) \
			- (1/(-l))*b*(1j*(k)*abs(n-k)*(uf[k]*uf[n-k]))
			#z = z - (1j*k*uf[k]*uf[n-k]) - b*(1j*k*(uf[k]*uf[n-k]))

	return z

def gCLM_RK4_FFT(M,N,K,X,Tstep,Xstep,uf,b,A,Tf):
	u = []
	ufourier = uf
	#.tolist() gives the transpose of what we want
	ufourier = np.transpose(ufourier).tolist()
	ufourier.append(np.zeros(2*N+2*K+1).tolist())
	adap_bound = A
	t = Tstep
	Tst = Tstep
	m=0
	#The matrix of intermediate RK values
	xi = np.zeros((2*(N+K)+1,3),dtype = complex)
	#Implements the RK4 method:
	#for m in range(M-1):
	#Change the original for loop to a while loop for adaptive code.
	print("The current Runge-Kutta time step:")
	while t< (Tf):
		udx = []
		#This will stop the evaluation if the values of ufourier get too large.
		print(m, " out of ", len(ufourier))
		for i in ufourier[m]:
			if (abs(np.real(i))+abs(np.imag(i)))>100:
				return ufourier
		#######################################################################
		#The Runge-Kutta Method
		
		for n in range(-N,N + 1):
			xi[n,0]=ufourier[m][n]+Tst*(.5)*gCLMfourier(ufourier[m][:],n,N,b)

		for n in range(-N,N + 1):
			xi[n,1]=ufourier[m][n]+Tst*(.5)*(gCLMfourier(xi[:,0],n,N,b))

		for n in range(-N,N + 1):
			xi[n,2]=ufourier[m][n]+Tst*(1)*(gCLMfourier(xi[:,1],n,N,b))

		for n in range(-N,N + 1):
			ufourier[m+1][n] = ufourier[m][n] \
			+ Tst*(1.0/6.0)*(gCLMfourier(ufourier[m][:],n,N,b)) \
			+ Tst*(1.0/3.0)*(gCLMfourier(xi[:,0],n,N,b))

			ufourier[m+1][n] = ufourier[m+1][n] \
			+ Tst*(1.0/3.0)*(gCLMfourier(xi[:,1],n,N,b)) \
			+ Tst*(1.0/6.0)*(gCLMfourier(xi[:,2],n,N,b))
		#######################################################################
		#Compute u and udx
		
		thisthing = np.conj(np.array(ufourier[m]))
		u.append(np.real(spf.fft(np.delete(thisthing,X/2,0))))

		for x in range(X-1):
			udx.append((u[m][x+1] - u[m][x-1])/(2*Xstep))
		#udx[m][X-1] = (u[m][0] - u[m][X-2])/(2*Xstep)
		udx = np.real(udx)
		#######################################################################
		#Decide whether or not to lower the step size.
		#Each time adap_bound is breached, refine the mesh. Stop if udx
		#gets too large.
		if (np.amax(np.absolute(udx))> adap_bound)&(np.amax(np.absolute(udx))<10):
			#change Tst 
			print("\n")
			print("udx = %f"%np.amax(np.absolute(udx)))
			print("Inserting new time mesh points")
			adap_bound = adap_bound*1.5
			temp = (Tf-t)/Tst
			Tst = Tst*(.5)
			print("Appending %d rows of zeros to ufourier"%int(np.round((Tf - t)/Tst-temp)))
			for i in range(int(np.round((Tf - t)/Tst-temp))):
				ufourier.append(np.zeros(2*N+2*K+1).tolist())
			print("The new total number of time steps is: %d"%len(ufourier))
			print("\n")
		t+=Tst
		m+=1
	print("The new total number of time steps is: %d"%len(ufourier))
	Mtime = m


	#u.append(spf.fft(np.delete(np.conj(np.transpose(ufourier[:,M-1])),X/2,0)))
	return [ufourier,u,Mtime]
