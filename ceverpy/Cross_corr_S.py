# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 11:32:10 2017

@author: leo
"""

# PICK THE S-ARRIVAL TIMES FOR THE DRY AND SATURATED WAVEFORMS USING DTW
"""
Created on Fri Mar 10 15:54:27 2017

@author: edur409
"""

COORDS = []

#Subroutine to pick values from the active plot
def onpick(event):
    global COORDS
    thisline = event.artist
    xdata = thisline.get_xdata()
    ydata = thisline.get_ydata()
    ind = event.ind
    points = np.array([xdata[ind], ydata[ind]])
    print('onpick points:', points)
    COORDS.append((xdata[ind], ydata[ind]))
    return COORDS

def crosscorr_lags(A,B):
    C=np.correlate(A,B,mode='full')
    lags=np.linspace(-len(A),len(A),len(C))
    return lags,C

#Import Modules
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
from scipy.signal import hilbert
from matplotlib.ticker import NullFormatter, FormatStrFormatter
import mcerp3 as mc
from . import RockPhysics as rp
#from waveletFunctions import wavelet

def cross_corr_s():
    global COORDS
    rpy2.robjects.numpy2ri.activate()

    #Font size for plots
    font= {'size' : 18}
    plt.rc('font',**font)

    # Set up our R namespaces
    R = rpy2.robjects.r
    DTW = importr('dtw') #import the "dtw" function from R

    #Choose the sample's folder
    well=input("Type name of sample (e.g. 'NM11_2087_4A'): \n")

    #Use waveform or envelope for the DTW match
    method=int(input("Use property of signal to math: Press 1 (waveform), Press 2 (envelope), Press 3 (both) \n"))

    #1) READ THE WAVEFORMS
    #1.1) READ THE DRY SAMPLE WAVEFORM FIRST
    timeSd,T,Sd=rp.plot_csv('./'+well+'/'+well+'d/tek0002ALL.csv',0)

    print('Choose the beginning and end that you want to compare')
    #Choose the two extremes to compare
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('click on points')
    COORDS=[]
    ax.plot(timeSd,Sd,'-',picker=5),plt.xlabel('Time ($\mu$s)')
    plt.grid(color='0.5')
    cid = fig.canvas.mpl_connect('pick_event', onpick)
    #Work around for plt.show(block=True) in IPython using Spyder
    plt.show(block=True)#plt.show(block=True)
    #Pick the times
    t_id=np.min(COORDS[0][0])
    t_fd=np.min(COORDS[1][0])

    #Indices to crop the dry sample time vector
    i_d=np.min(np.where(timeSd>=t_id))
    f_d=np.min(np.where(timeSd>=t_fd))

    #1.2) READ THE SATURATED SAMPLE WAVEFORM SECOND
    timeSs,T,Ss=rp.plot_csv('./'+well+'/'+well+'s/tek0002ALL.csv',0)

    print('Choose the beginning and end that you want to compare')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('click on points')
    COORDS=[]
    ax.plot(timeSs,Ss,'-',picker=5),plt.xlabel('Time ($\mu$s)')
    plt.grid(color='0.5')
    fig.canvas.mpl_connect('pick_event', onpick)
    plt.show(block=True)#plt.show(block=True)

    t_is=np.min(COORDS[0][0])
    t_fs=np.min(COORDS[1][0])

    #Indices to crop for the saturated sample time vector
    i_s=np.min(np.where(timeSs>=t_is))
    f_s=np.min(np.where(timeSs>=t_fs))

    #2) DO THE COMPARISON USING THE DTW ALGORITHM
    #THE DYNAMIC TIME WARPING ALGORITHM BEGINS HERE!!!!
    idxt=timeSd[i_d:f_d] #Indices of the template function (dry sample waveform)
    idxq=timeSs[i_s:f_s] #Indices of the query function (saturated sample waveform)

    #2.1) CHOOSE WHAT INFORMATION TO COMPARE (EITHER WAVEFORMS OR ENVELOPES)
    if method==1: #Method 1 compares the waveforms
        template=Sd[i_d:f_d]/np.max(np.abs(Sd[i_d:f_d]))
        query=Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s]))
        suffix='wav' #Suffix for the files
    elif method==2: #Method 2 compares the energy of the envelope of the waveforms
        template=np.abs(hilbert(Sd[i_d:f_d]/np.max(np.abs(Sd[i_d:f_d]))))
        query=np.abs(hilbert(Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s]))))
        suffix='env' #Suffix for the files
    elif method==3: #Method 3 compares envelope and waveforms
        template=Sd[i_d:f_d]/np.max(np.abs(Sd[i_d:f_d]))
        query=Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s]))
        templateh=np.abs(hilbert(Sd[i_d:f_d]/np.max(np.abs(Sd[i_d:f_d]))))
        queryh=np.abs(hilbert(Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s]))))
        suffix='both' #Suffix for the files
    else:
        print("Please, follow the instructions. Try again!")
        exit()
        
    #Cross-correlation
    dtt=np.mean(np.diff(idxt)) #sampling rate of template
    dtq=np.mean(np.diff(idxq)) #sampling rate of query
    if dtt<=dtq:
        dt=dtt
        query=np.interp(idxt,idxq,query)
        time=idxt
    elif dtt>dtq:
        dt=dtq
        template=np.interp(idxq,idxt,template)
        time=idxq
    #Both have the same support (domain)
    ws=15 #microsecs
    w=np.int(np.round(ws/dt))
    p=np.arange(np.min(time),np.max(time),5)
    pos=np.min(np.where(time>=np.max(p)))
    #slide the correlation window 
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
        Th=np.append(Th,np.min(timew)+(np.max(timew)-np.min(timew))/2)

    print('Mode of Tau: ',sp.stats.mode(Taus).mode[0])
    C=np.correlate(template,query,mode='same') #dry is template, sat is query
    pmax=np.where(C==np.max(C))[0]
    lag=np.linspace(-dt*len(C)/2,dt*len(C)/2,len(C))
    tau=lag[pmax] #delta time between waveforms
    #Cross-correlation
    plt.plot(lag,C,lw=2)
    plt.axvline(tau,ymin=1.1*np.min(C),ymax=1.1*np.max(C),c='b',lw=2,label='$\Delta t$: '+np.str(np.round(tau[0],decimals=2))+' $\mu s$')
    plt.axvline(sp.stats.mode(Taus).mode[0],ymin=1.1*np.min(C),ymax=1.1*np.max(C),c='r',lw=2,label='$\Delta t_w$: '+np.str(np.round(sp.stats.mode(Taus).mode[0],decimals=2))+' $\mu s$')
    plt.xlabel('lag ($\mu s$)'),plt.ylabel('Correlation Coefficient')
    plt.legend()
    plt.grid(color='0.5')
    plt.savefig('./'+well+'/'+well+'s/Vs_sat_Xcorr.pdf',bbox_inches='tight')
    plt.show(block=True)

    plt.figure('Windowed waveforms')
    plt.plot(time,template,label='S (dry)',color='r',lw=2)
    plt.plot(time,query,label='S (sat)',color='b',lw=2)
    plt.grid(color='0.5')
    plt.legend()
    plt.xlabel('time ($\mu s$)'),plt.ylabel('A.U.')
    plt.savefig('./'+well+'/'+well+'s/Vs_waveforms_xcorr.pdf',bbox_inches='tight')

    #Import the lengths measured with the calliper
    Lengths=np.loadtxt('./'+well+'/'+well+'_Lengths.txt',delimiter=',')#np.array([5.256, 5.250, 5.254, 5.254, 5.252, 5.252, 5.258, 5.265, 5.255, 5.252])
    #Calculate mean and standard deviation
    L,StdL=rp.Length_sample(Lengths)

    tS=np.loadtxt('./'+well+'/'+well+'d/S_time_Picks_dry_both.out',delimiter=',')
    tSm=np.mean(tS) 
    dtS=np.std(tS)

    if tau>0:
        Vsat=1e4*L/(tSm-np.abs(tau))
        dVsat=Vsat*np.sqrt((StdL/L)**2+(dtS/(tSm-tau))**2)
    elif tau<0:
        Vsat=1e4*L/(tSm+np.abs(tau))
        dVsat=Vsat*np.sqrt((StdL/L)**2+(dtS/(tSm+tau))**2)
    #print values in Vsat
    #Velocity standard way
    VelS,DVelS,Vs_down,Vsmc,Vs_up=rp.Velocity_S(L,tSm,StdL,dtS)
    print('S-Velocity (dry) (2*std): '+np.str(np.round(VelS-2*DVelS,decimals=0))+' < '+np.str(np.round(VelS,decimals=0))+' < '+np.str(np.round(VelS+2*DVelS,decimals=0))+' m/s')
    print('S-Velocity (sat.) (2*std): '+np.str(np.round(Vsat[0]-2*dVsat[0],decimals=0))+' < '+np.str(np.round(Vsat[0],decimals=0))+' < '+np.str(np.round(Vsat[0]+2*dVsat[0],decimals=0))+' m/s')

    #Save the Velocity from the cross-correlation
    np.savetxt('./'+well+'/'+well+'s/Velocities_S_sat_Xcorr.out', [Vsat,dVsat], fmt='%1.2f',delimiter=',',header='Vs (sat), StdVs (sat)')    
        
    #plt.figure('Correlation')
    #plt.plot(time,template),plt.plot(time,query)
    #plt.plot(Th,Corr/np.max(Corr))
