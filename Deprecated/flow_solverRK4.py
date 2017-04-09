#Solves for the Lagrangian flow of u.
#In terms of the global variables in main we'll usually have:
#u = u
#M = Mtime
#N = 2*Nfourier +2*Kpad + 1
#X = Xspatial
#b = b
#############################################################################
import numpy as np
import math 
from scipy.interpolate import interp1d
#############################################################################
def flow_RK4(M,X,Tstep,Xstep,u):
	u = np.real(u)
	eta = np.zeros((X+1,M+1))
	eta[:,0] = np.linspace(0,2*(math.pi),num = X+1)
	#The matrix of intermediate RK values
	xi = np.zeros((X+1,3))
	#Implements the RK4 method:
	for m in range(M):
		f = interp1d(eta[:,0],u[:,m])

		xi[:,0]=eta[:,m]+Tstep*(.5)*f(eta[:,m])

		xi[:,1]=eta[:,m]+Tstep*(.5)*f(xi[:,0])

		xi[:,2]=eta[:,m]+Tstep*(1)*f(xi[:,1])

		eta[:,m+1] = eta[:,m] \
			+ Tstep*(1.0/6.0)*f(eta[:,m]) \
			+ Tstep*(1.0/3.0)*f(xi[:,0])

		eta[:,m+1] = eta[:,m+1] \
			+ Tstep*(1.0/3.0)*f(xi[:,1])\
			+ Tstep*(1.0/6.0)*f(xi[:,2])

	return eta
