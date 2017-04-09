#This is the python version of the matlab program De-Gregorio Solver Suite v01.
#This contains code for studying the solutions to the generalized 
#Constantin-Lax-Majda equation:

#w_t + uw_x + bwu_x = 0, w = Hu_x 
################################################################################
#Imported packages
import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
import math 
import scipy as sp
import scipy.fftpack as spf
from scipy.integrate import odeint
#Other imported files
import gCLM_Solver_RK4_Adap as gCLMS
#import flow_solverRK4 as flows
from scipy.interpolate import interp2d

import time

starttime = time.time()
################################################################################
#Global Variables to be set by user:

#The number of positive fourier coefficients. There will then
#be 2*Nfourier + 1 not necessarily zero fourier coefficients total
Nfourier = 128
#Initial time
Tinit = 0.0 #make sure this is a float
#Final time
Tfin = .8 #make sure this is a float
#Number of starting time steps
Mtime = 512
Mtimearr = list(np.linspace(0,Tfin, num = Mtime))
#INSERT
Adap_param = 3
#Number of spatial steps. Should be no less than Nfourier*4 
#to prevent aliasing. Also generally best to make a power of 2.
Xspatial = 512
#The number of zeros to pad with to prevent aliasing (e.g. 3/2 rule) and make 
#fft output according to Xspatial resolution. 
Kpad = int(float(Xspatial)/2) - Nfourier 
#Initial Fourier coefficients of u, the solution (make sure they're floats)
ufourier = np.zeros((2*Nfourier+2*Kpad+1,Mtime),dtype = complex)
ufourier[1,0] =-.5*1j
ufourier[-1,0] = .5*1j
#The parameter b in the gCLM equation
b = 2.0

################################################################################
#Global variables that will be computed:

#The matrix giving the values of u
#u = np.zeros((Xspatial,Mtime))
#Time step size
Tstep = (Tfin-Tinit)/float(Mtime)
#Spatial step size
Xstep = (2*math.pi)/float(Xspatial)

################################################################################
#Implement the Fourier-Galerkin method.
[ufourier,u,Mtime] = gCLMS.gCLM_RK4_FFT(Mtime,Nfourier,Kpad,Xspatial,Tstep,Xstep,
	ufourier,b,Adap_param,Tfin)
print("Runge-Kutta method and fft completed")
#We transpose so that we fft across the right way. Additionally scipy uses
#sum x*exp(-2pi i n) rather than sum x*exp(2pi i n) for fft as I thought, so 
#I must conjugate first.

u = np.array(u)
u = np.real(u)
u = np.transpose(u)

ufourier = np.transpose(np.array(ufourier))

#The vorticity w:
wfourier = np.zeros((2*Nfourier+2*Kpad+1,Mtime),dtype = complex)

#Compute the fourier coefficients of the vorticity w:
for iindex in range(Mtime):
	for jindex in range(int(Nfourier/2)+1):
		wfourier[jindex,iindex] = float(jindex)*ufourier[jindex,iindex]
		wfourier[-jindex,iindex] = float(jindex)*ufourier[-jindex,iindex]

#Compute w:
w = spf.fft(np.conj(np.transpose(wfourier)))
w = np.real(w)
w = np.transpose(w)

#Compute the energy (H^1/2 norm of u)
energy = np.zeros(Mtime)

for m in range(Mtime):
	for (a,b) in zip(ufourier[:,m],wfourier[:,m]):
		energy[m] += np.real(a*np.conj(b))
################################S################################################
#Solve for the flow eta
eta = np.zeros((Xspatial,Mtime))
#Space partition
xspace = np.linspace(0,2*(math.pi)-Xstep,num = Xspatial)

eta[:,0] = xspace

#Time partition for solving flow eqn
tspace=np.linspace(0,Tfin,num = Mtime)

f = interp2d(xspace,tspace,np.transpose(u))
for x in range(Xspatial):
  	temp= sp.integrate.odeint(f,eta[x,0],tspace)
  	eta[x,:]= np.reshape(temp,(Mtime))
################################################################################
#Take the derivative of eta
etadx = np.zeros((Xspatial,Mtime))
for m in range(Mtime):     
    for x in range(Xspatial-1):
        
        etadx[x,m]= (eta[x+1,m]-eta[x,m])/(Xstep)
        
    etadx[-1,m] =((eta[0,m]+2*math.pi)-eta[-1,m])/float(Xstep)
        #etadx[Xspatial+1,m]=etadx[1,m];
################################################################################
#Generate the angle in the Ermakov-Pinney representation
etaxint = np.zeros((Xspatial,Mtime))

for p in range(Xspatial-1):
	for k in range(Mtime):
		for i in range(k):
			etaxint[p, k] = etaxint[p, k] + (w[p, 0]/(etadx[p, i]**2.0))*Tstep;

		etaxint[p, k]= etaxint[p, k] + 2*math.pi*(float(p)/float(Xspatial));
		etaxint[-1, k] = etaxint[0, k]+2*math.pi;

################################################################################
#Print output to various files. It's easier for me to make the space steps rows 
#rather than columns, hence the transpose.
np.savetxt("Output/ufourier.txt",np.transpose(ufourier),fmt = '%.2f')
np.savetxt("Output/uout.txt",np.transpose(u),fmt = '%.2f')
np.savetxt("Output/eta.txt",np.transpose(eta),fmt = '%.2f')
np.savetxt("Output/etadx.txt",np.transpose(etadx),fmt = '%.2f')
np.savetxt("Output/wfourier.txt",np.transpose(wfourier),fmt = '%.2f')
np.savetxt("Output/wout.txt",np.transpose(w),fmt = '%.2f')
np.savetxt("Output/etaxint.txt",np.transpose(etaxint),fmt = '%.2f')
np.savetxt("Output/energy.txt",np.transpose(energy),fmt = '%.2f')

fintime =time.time()

print("Total Time: %f " %  (fintime-starttime))
