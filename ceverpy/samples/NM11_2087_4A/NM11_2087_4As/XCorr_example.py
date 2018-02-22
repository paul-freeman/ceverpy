# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 12:34:00 2017

@author: leo
"""
#Import Modules
import numpy as np
from scipy.stats import mode
import matplotlib.pyplot as plt

M=np.loadtxt('Waveforms.out',skiprows=1,delimiter=',')
time=M[:,0] #Time 
template=np.convolve(M[:,1],np.ones(5)/5,'same') #P dry
query=np.convolve(M[:,2],np.ones(5)/5,'same') #P sat.

dt=np.mean(np.diff(time)) #sampling interval in microsecs

plt.figure('Waveforms')
plt.plot(time,template,label='P (dry)')
plt.plot(time,query,label='P (sat.)')
plt.xlabel('time ($\mu s$)')
plt.ylabel('A.U.')
plt.legend()

WS=[5]
TAU=[]
for j in range(0,len(WS)):
    ws=WS[j] #Window size in microsecs
    w=np.int(np.round(ws/dt)) #window size in number of samples
    p=np.arange(np.min(time),np.max(time),5)
    pos=np.min(np.where(time>=np.max(p)))
    #slide the correlation window of length ws
    Corr=[]
    Taus=[]
    Th=[] #half time of the sliding window
    for i in range(0,len(time)-w):
        timew=time[i:i+w]
        C=np.correlate(template[i:i+w],query[i:i+w],mode='same') #dry is template, sat is query
        Corr=np.append(Corr,np.max(C))
        lags=np.linspace(-dt*len(C)/2,dt*len(C)/2,len(C))    
        pmax=np.where(C==np.max(C))[0]
        Taus=np.append(Taus,lags[pmax])
        Th=np.append(Th,np.min(timew)+(np.max(timew)-np.min(timew))/2) #mid-point of window
    
    print(('Mode of Tau for window size '+np.str(ws)+': ',sp.stats.mode(Taus).mode[0]))
    C=np.correlate(template,query,mode='same') #dry is template, sat is query
    pmax=np.where(C==np.max(C))[0]
    lag=np.linspace(-dt*len(C)/2,dt*len(C)/2,len(C))
    tau=lag[pmax] #delta time between waveforms
    TAU=np.append(TAU,tau)
    
#Simple Cross-correlation
plt.figure('Cross-correlation')
plt.plot(lag,C)
plt.axvline(tau,ymin=1.1*np.min(C),ymax=1.1*np.max(C),c='r',lw=2,label='$\Delta$t: '+np.str(np.round(tau[0],decimals=2))+' microsecs')
plt.axvline(sp.stats.mode(Taus)[0],ymin=1.1*np.min(C),ymax=1.1*np.max(C),c='r',lw=2,label='$\Delta t_w$: '+np.str(np.round(sp.stats.mode(Taus).mode[0],decimals=2))+' microsecs')
plt.xlabel('lag ($\mu$s)'),plt.ylabel('Correlation Coefficient')
plt.legend()
plt.grid('on')

#Windowed cross-correlation vs. middle point of moving window
plt.figure('Tau with time')
plt.plot(Th,Taus)
plt.xlabel('time ($\mu s$)'),plt.ylabel('tau'),plt.grid('on')
