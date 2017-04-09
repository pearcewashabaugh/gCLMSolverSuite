#This code plots the Ermakov-Pinney trajectory
#WARNING: Make sure that frames = Mtime

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import style

style.use('ggplot')

Mtime = 128
fig = plt.figure()
ax1 = plt.axes(xlim = (-2,2), ylim = (-2,2))
line, = ax1.plot([],[],lw = 1.5, color = "#3C4644")

#import radii of trajectory
rad = open('Output/etadx.txt','r').read()
#import angles
theta = open('Output/etaxint.txt','r').read()

radlines = rad.split('\n')
thetalines = theta.split('\n')

#Make lists of radii and angles
rads = []
thetas = []
for lin in radlines[:-1]:
    x  = lin.split(' ')
    rads.append(x)

for lin in thetalines[:-1]:
    y = lin.split(' ')
    thetas.append(y)

polar = []
for i in range(len(rads)):
    polar.append(list(zip(rads[i],thetas[i])))

#initialization function: plot the background of each frame
def init():
    line.set_data([],[])
    return line, 

#animation function. This is called sequentially
def animate(i):
    """temp = (i - Mtime/2.0)/(Mtime/2.0)
    if i>=Mtime/2:
        ax1.set_xlim(-2 + temp*1.9, 2 - temp*1.9)
        ax1.set_ylim(-2 + temp*1.9, 2 - temp*1.9)"""
    xs = []
    ys = []
    for radtheta in polar[int(i)]:
        xs.append(float(radtheta[0])*np.cos(float(radtheta[1])))
        ys.append(float(radtheta[0])*np.sin(float(radtheta[1])))
    line.set_data(xs,ys)
    return line,

#call the animator. 
#WARNING: Make sure that frames = Mtime
ani = animation.FuncAnimation(fig, animate, init_func = init, 
    frames = Mtime, interval = 80)
ani.save('/home/pearce/Desktop/ermpin.gif', writer='imagemagick', fps = 30)
plt.show()

#Partially taken from:

"""
Matplotlib Animation Example:

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
