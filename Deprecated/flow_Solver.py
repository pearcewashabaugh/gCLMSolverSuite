#Given u solves for the flow eta.

eta = np.zeros((2*Xspatial+1,Mtime+1))
#Space partition
xspace = np.linspace(0,2*(math.pi),num = 2*Xspatial+1)

eta[:,0] = xspace

#Time partition for solving flow eqn
tspace=np.linspace(0,Tfin,num = Mtime+1)

f = interp2d(xspace,tspace,np.transpose(u))
for x in range(2*Xspatial+1):
  	temp= sp.integrate.odeint(f,eta[x,0],tspace)
  	eta[x,:]= np.reshape(temp,(Mtime+1))