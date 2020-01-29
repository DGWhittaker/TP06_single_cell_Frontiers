#This script can be used to plot short APs
#Author: Dominic Whittaker
#Date: 29th January 2020
#Compile with python plot_APs_short.py

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

a = np.loadtxt( 'APs/ENDO.txt', unpack=True )
b = np.loadtxt( 'APs/MCELL.txt', unpack=True )
c = np.loadtxt( 'APs/EPI.txt', unpack=True )

a2 = np.loadtxt( 'APs/ENDO_SHORT.txt', unpack=True )
b2 = np.loadtxt( 'APs/MCELL_SHORT.txt', unpack=True )
c2 = np.loadtxt( 'APs/EPI_SHORT.txt', unpack=True )

ENDO_t = a[0] - 98975
ENDO_v = a[1]

MCELL_t = b[0] - 98975
MCELL_v = b[1]

EPI_t = c[0] - 98975
EPI_v = c[1]

ENDO_S_t = a2[0] - 98975
ENDO_S_v = a2[1]

MCELL_S_t = b2[0] - 98975
MCELL_S_v = b2[1]

EPI_S_t = c2[0] - 98975
EPI_S_v = c2[1]

fig = p.figure(1, figsize=(6,3)) #(8,6) seems to be the default
fig.set_facecolor('white') #this changes the background colour from the default gray/blue to white

#*********************NORMAL********************************************#
ax1 = fig.add_subplot(121)
ax1.plot( ENDO_t, ENDO_v, label = 'WT', color='red', linewidth=1.5 )
ax1.plot( MCELL_t, MCELL_v, label = 'V307L', color='forestgreen', linewidth=1.5 )
ax1.plot( EPI_t, EPI_v, label = 'WT', color=(0.11765,0.23529,1.0), linewidth=1.5 )
ax1.set_xlabel( 'Time (ms)' )
ax1.set_ylabel( 'Voltage (mV)' )
ax1.axis([0, 500, -100, 50])
ax1.spines['right'].set_visible(False) #removes top and right borders
ax1.spines['top'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
p.xticks(np.arange(0,501,250))
p.yticks(np.arange(-100,50+1,50))
p.title('Normal') 

#*********************SHORT********************************************#
ax1 = fig.add_subplot(122)
ax1.plot( ENDO_S_t, ENDO_S_v, label = 'WT', color='red', linewidth=1.5 )
ax1.plot( MCELL_S_t, MCELL_S_v, label = 'V307L', color='forestgreen', linewidth=1.5 )
ax1.plot( EPI_S_t, EPI_S_v, label = 'WT', color=(0.11765,0.23529,1.0), linewidth=1.5 )
ax1.set_xlabel( 'Time (ms)' )
ax1.axis([0, 500, -100, 50])
ax1.spines['right'].set_visible(False) #removes top and right borders
ax1.spines['top'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
p.xticks(np.arange(0,501,250))
p.yticks(np.arange(-100,50+1,50))
[label.set_visible(False) for label in ax1.get_yticklabels()]
p.legend(['ENDO', 'MCELL', 'EPI'], loc='upper right', fontsize=9, frameon=False) #add frameon=False to remove the legend box
p.title('Short') 
p.tight_layout()
p.show()

p.close()