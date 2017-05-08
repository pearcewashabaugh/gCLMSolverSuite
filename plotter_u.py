#This code makes 2D plots of the desired variables such as u or w = omega
#WARNING: Make sure that frames = Mtime

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import style
import time
###############################################################################
# Main variables to be set by user. Note there are some other places to tweak,
# such as the the y limits of the plot.

# Import relevant variable from output of main.py
ufile = open('Output/uout.txt','r').read()

# Adjust Mtime and Xspatial to be the Mtime and Xspatial used in main.py. 
# This will get fixed in the future.
Mtime = 64
Xspatial = 128

###############################################################################

style.use('ggplot')

fig = plt.figure()
ax1 = plt.axes()
ax1 = plt.axes(xlim = (0,2*np.pi), ylim = (-1.5,1.5))
line, = ax1.plot([],[],lw = 1.2, color = "#3C4644")

ulines = ufile.split('\n')

u = []

for lin in ulines[:-1]:
    y = lin.split(' ')
    u.append(y)

#initialization function: plot the background of each frame
def init():
    line.set_data([],[])
    return line, 

xs = []

for j in range(Xspatial):
    xs.append((j/Xspatial)*2*np.pi)

#animation function. This is called sequentially

def animate(i):
    line.set_data(xs,u[i])
    return line,

#call the animator. 

ani = animation.FuncAnimation(fig, animate, init_func = init, 
    frames = Mtime, interval = 40)

#ani.save('put_save_path_here', writer='imagemagick', fps = 30)

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
