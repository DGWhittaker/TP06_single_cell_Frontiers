#This script can be used to plot APs
#Author: Dominic Whittaker
#Date: 29th January 2020
#Compile with python plot_APs.py

import numpy as np
import matplotlib.pyplot as p
from matplotlib import * 

p.rcdefaults()

p.rc('lines', markeredgewidth=2)
p.rc('lines', markersize=7)
p.rc('xtick.major', size=5)
p.rc('ytick.major', size=5) #changes size of y axis ticks on both sides
p.rc('xtick', direction='out')
p.rc('ytick', direction='out')

a = np.loadtxt( 'APs/EPI.txt', unpack=True )
b = np.loadtxt( 'APs/EPI_BASE.txt', unpack=True )

EPI_t = a[0] - 98975
EPI_v = a[1]

EPI_B_t = b[0] - 98975
EPI_B_v = b[1]

fig = p.figure(1, figsize=(3,3)) #(8,6) seems to be the default
fig.set_facecolor('white') #this changes the background colour from the default gray/blue to white

ax1 = fig.add_subplot(111)
ax1.plot( EPI_t, EPI_v, color=(0.11765,0.23529,1.0), linewidth=1.5 )
ax1.plot( EPI_B_t, EPI_B_v, color='lightseagreen', linewidth=1.5 )
ax1.set_xlabel( 'Time (ms)' )
ax1.set_ylabel( 'Voltage (mV)' )
ax1.axis([0, 500, -100, 50])
ax1.spines['right'].set_visible(False) #removes top and right borders
ax1.spines['top'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
p.xticks(np.arange(0,501,250))
p.yticks(np.arange(-100,50+1,50))
p.legend(['APEX', 'BASE'], loc='upper right', fontsize=9, frameon=False) #add frameon=False to remove the legend box
p.tight_layout()
p.show()

p.close()