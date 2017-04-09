#This code uses a 4th order Runge-Kutta method (see Iserles pg.41) to solve
#Burgers equation: 
#u_t + b uu_x = 0
################################################################################
import numpy as np
################################################################################
#In terms of the global variables in main we'll usually have:
#uf = ufourier
#M = Mtime
#N = 2*Nfourier +2*Kpad + 1
#X = Xspatial
#b = b

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

def gCLM_RK4(M,N,K,X,Tstep,Xstep,uf,b):
	ufourier = np.zeros((2*(N+K)+1,M),dtype = complex)
	ufourier[:,0] = uf[:,0]
	#The matrix of intermediate RK values
	xi = np.zeros((2*(N+K)+1,3),dtype = complex)
	#Implements the RK4 method:
	for m in range(M-1):
		#This will stop the evaluation if the values of ufourier get too large.
		for i in ufourier[:,m]:
			if (abs(np.real(i))+abs(np.imag(i)))>100:
				return ufourier

		for n in range(-N,N + 1):
			xi[n,0]=ufourier[n,m]+Tstep*(.5)*gCLMfourier(ufourier[:,m],n,N,b)

		for n in range(-N,N + 1):
			xi[n,1]=ufourier[n,m]+Tstep*(.5)*(gCLMfourier(xi[:,0],n,N,b))

		for n in range(-N,N + 1):
			xi[n,2]=ufourier[n,m]+Tstep*(1)*(gCLMfourier(xi[:,1],n,N,b))

		for n in range(-N,N + 1):
			ufourier[n,m+1] = ufourier[n,m] \
			+ Tstep*(1.0/6.0)*(gCLMfourier(ufourier[:,m],n,N,b)) \
			+ Tstep*(1.0/3.0)*(gCLMfourier(xi[:,0],n,N,b))

			ufourier[n,m+1] = ufourier[n,m+1] \
			+ Tstep*(1.0/3.0)*(gCLMfourier(xi[:,1],n,N,b)) \
			+ Tstep*(1.0/6.0)*(gCLMfourier(xi[:,2],n,N,b))

	return ufourier
