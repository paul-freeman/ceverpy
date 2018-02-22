# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 15:28:34 2017

@author: edur409
"""

def onpick(event):
    global coords
    thisline = event.artist
    xdata = thisline.get_xdata()
    ydata = thisline.get_ydata()
    ind = event.ind
    points = np.array([xdata[ind], ydata[ind]])
    print('onpick points:', points)
    coords.append((xdata[ind], ydata[ind]))
    return coords
    
#Import Modules
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
from scipy.signal import hilbert
from matplotlib.ticker import NullFormatter, FormatStrFormatter
import mcerp as mc
import RockPhysics as rp
import glob

#Read the waveforms 
#REORDER THE FILENAMES SO THE FIRST INDEX IS THE HIGHEST PRESSURE
filenames = sorted(glob.glob('*.csv'),key=lambda f: int(filter(str.isdigit, f)),reverse=True)#

Pc=[]
phases=[]
states=[]
cycles=[]    
for i in range(0,len(filenames)):
    name,Pconf,phase,state,cycle=rp.extract_parameters(filenames[i])
    Pc=np.append(Pc,Pconf)
    phases=np.append(phases,phase)
    states=np.append(states,state)
    cycles=np.append(cycles,cycle)

#CHOOSE THE PHASE TO PICK
indices=np.where((phases=='PP')& (cycles=='d2')& (states=='sat3000'))[0]

#Now we're ready to go over the indices and do all the picking
for k in range(0,len(indices)):
    name=filenames[indices[k]]
    #load a file
    timeP, P=rp.load_csv(name,0)
    #plt.plot(timeP,np.convolve(P,np.ones(30)/30,'same'),lw=2)

    print("Please, check the smoothing is satisfactory by zooming in!!!")
    #Plot smoothed waveform
    plt.figure('P-waveform')
    #plt.plot(timeP,P,label='P')
    l=60 #Length of moving average window
    P_smooth=np.convolve(P,np.ones(l)/l,'same')
    plt.title('Check smoothing is ok on '+name+'')
    plt.plot(timeP,P,label='P')
    plt.plot(timeP,P_smooth,lw=2,label='P smoothed')
    plt.show(block=True)
    
    print("Zoom in, deactivate the magnifying glass button, and pick")
    #Pick the times from the graph   
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('click on points')
    #
    coords=[]
    ax.plot(timeP,P_smooth,lw=2,label='P smooth', picker=5)  # 5 points tolerance
    plt.grid('on')
    
    fig.canvas.mpl_connect('pick_event', onpick)
    
    plt.show(block=True) #Blocks exectution of the code until the figure is closed
    
    #Extract the times from the global variable coords and plot their distribution
    tP=rp.P_values(coords) #P-arrival times
    plt.xlabel('time (microsecs.)'),plt.ylabel('Normalized histogram')
    plt.savefig(''+phases[indices[k]]+'_time_Picks_sat3000_'+np.str(np.int(Pc[indices[k]]))+'.pdf',bbox_inches='tight')
    np.savetxt(''+phases[indices[k]]+'_time_Picks_sat3000_'+np.str(np.int(Pc[indices[k]]))+'.out',tP) #save the picked times
    print('The sample '+name+' is ready!!!' )
    
fig = plt.figure('Vp (PP) vs. Pressure')
ax = fig.add_subplot(111)
dP=50*np.ones(len(indices))  
Vp=[]
dVp=[]
for i in range(0,len(indices)):
    tP=np.loadtxt(''+phases[indices[i]]+'_time_Picks_sat3000_'+np.str(np.int(Pc[indices[i]]))+'.out')
    Vp=np.append(Vp,np.mean(tP))    
    dVp=np.append(dVp,np.std(tP))
    ax.scatter(Pc[indices[i]],np.mean(tP))
ax.errorbar(Pc[indices], Vp, yerr=[2*dVp, 2*dVp], xerr=[2*dP, 2*dP], fmt='o')
    