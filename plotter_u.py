#This code makes 2D plots of the desired variables such as u or w = omega
#WARNING: Make sure that frames = Mtime

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import style
import time

style.use('ggplot')

fig = plt.figure()
ax1 = plt.axes()
ax1 = plt.axes(xlim = (0,2*np.pi), ylim = (-3,3))
line, = ax1.plot([],[],lw = 1.2, color = "#3C4644")

#import radii of trajectory
ufile = open('Output/wout.txt','r').read()

ulines = ufile.split('\n')


#Make lists of radii and angles
u = []

for lin in ulines[:-1]:
    x  = lin.split(' ')
    u.append(x[:-1])

#initialization function: plot the background of each frame
def init():
    line.set_data([],[])
    return line, 

xs = []

for j in range(256):
    xs.append((j/256)*2*np.pi)

#animation function. This is called sequentially

def animate(i):
    line.set_data(xs,u[i])
    return line,

#call the animator. 
#WARNING: Make sure that frames = Mtime
ani = animation.FuncAnimation(fig, animate, init_func = init, 
    frames = 128, interval = 40)
ani.save('/home/pearce/Desktop/wout.gif', writer='imagemagick', fps = 30)
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
