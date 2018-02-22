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
import matplotlib.pyplot as plt
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
from scipy.signal import hilbert
from matplotlib.ticker import NullFormatter, FormatStrFormatter
import mcerp3 as mc
from . import RockPhysics as rp
#from waveletFunctions import wavelet

def pick_vp_dtw():
    rpy2.robjects.numpy2ri.activate()
    global COORDS
    #Font size for plots
    font= {'size' : 18}
    plt.rc('font',**font)

    # Set up our R namespaces
    R = rpy2.robjects.r
    DTW = importr('dtw') #import the "dtw" function from R

    #Choose the sample's folder
    well=input("Type name of sample (e.g. 'NM11_2087_4A'): \n")

    #1) READ THE WAVEFORMS
    #1.1) READ THE DRY SAMPLE WAVEFORM FIRST
    timeSd,T,Sd=rp.plot_csv('./'+well+'/'+well+'d/tek0001ALL.csv',0)

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
    timeSs,T,Ss=rp.plot_csv('./'+well+'/'+well+'s/tek0001ALL.csv',0)

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

    #2.1) CHOOSE INFORMATION TO COMPARE (WAVEFORMS OR ENVELOPES)
    template=Sd[i_d:f_d]/np.max(np.abs(Sd[i_d:f_d]))
    query=Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s]))
    templatea=(np.angle(hilbert(Sd[i_d:f_d]/np.max(np.abs(Sd[i_d:f_d])))))/np.pi
    querya=(np.angle(hilbert(Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s])))))/np.pi
    templateh=np.abs(hilbert(Sd[i_d:f_d]/np.max(np.abs(Sd[i_d:f_d]))))
    queryh=np.abs(hilbert(Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s]))))
    suffix='both' #Suffix for the files



    # Calculate the alignment vector and corresponding distance (DTW)
    alignment = R.dtw(query, template, keep=True)
    alignmenth = R.dtw(queryh, templateh, keep=True)
    alignmenta = R.dtw(querya, templatea, keep=True)
    dist = alignment.rx('distance')[0][0] #The lower the distance of alignment, the better the match
    disth = alignmenth.rx('distance')[0][0]
    dista = alignmenta.rx('distance')[0][0]

    print('Distance of the DTW algorithm (waveform): ',dist)
    print('Distance of the DTW algorithm (phase): ',dista)
    print('Distance of the DTW algorithm (envelope): ',disth)

    #Indices of alignment
    I1=alignment.rx('index1')[0]
    I2=alignment.rx('index2')[0]
    I1h=alignmenth.rx('index1')[0]
    I2h=alignmenth.rx('index2')[0]
    I1a=alignmenta.rx('index1')[0]
    I2a=alignmenta.rx('index2')[0]


    print('Close the figures to continue running the code...')

    #3) PLOT THE OUTPUTS OF THE DTW ALGORITHM
    #PLOTTING THE OUTPUTS OF THE DTW ALGORITHM BEGINS HERE!!!
    #Plots
    #Plot #1: points of match between the waveforms
    plt.figure('Points of match')
    plt.plot(idxt,template,label='Pd',c='b')
    plt.plot(idxq,query,label='Ps',c='g')
    plt.axis('tight')
    plt.legend()
    plt.xlabel('time ($\mu$s)')
    #Plot the matching points
    for i in np.arange(0,len(I1)-1,20):
        plt.plot([idxt[np.int(I2[i]-1)],idxq[np.int(I1[i]-1)]],[template[np.int(I2[i]-1)],query[np.int(I1[i]-1)]],'r-',lw=0.5)

    #Save the figure
    manager= plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.show()
    plt.savefig('./'+well+'/'+well+'d/DTW_Match_Pwaves_'+suffix+'.pdf',bbox_inches='tight')
    plt.savefig('./'+well+'/'+well+'s/DTW_Match_Pwaves_'+suffix+'.pdf',bbox_inches='tight')
    plt.show(block=True)#plt.show(block=True)

    #Plot #2 (mostly unnecessary): Compare a function over the other
    #plt.figure('One over the other')
    qo=[]
    qoa=[]
    qoh=[]
    to=[]
    toa=[]
    toh=[]
    idxto=[]
    idxtoa=[]
    idxtoh=[]
    idxqo=[]
    idxqoa=[]
    idxqoh=[]
    for i in range(0,len(I1)):
        idxto=np.append(idxto,idxt[np.int(I2[i])-1]) #time vector arranged by index
        idxqo=np.append(idxqo,idxq[np.int(I1[i])-1]) #time vector arranged by index
        qo=np.append(qo,query[np.int(I1[i])-1]) #Query function sampled by the index
        to=np.append(to,template[np.int(I2[i])-1]) #Template function sampled by the index
    for i in range(0,len(I1h)):    
        idxtoh=np.append(idxtoh,idxt[np.int(I2h[i])-1]) #time vector arranged by index
        idxqoh=np.append(idxqoh,idxq[np.int(I1h[i])-1]) #time vector arranged by index
        qoh=np.append(qoh,queryh[np.int(I1h[i])-1]) #Query function sampled by the index
        toh=np.append(toh,templateh[np.int(I2h[i])-1]) #Template function sampled by the index
    for i in range(0,len(I1a)):    
        idxtoa=np.append(idxtoa,idxt[np.int(I2a[i])-1]) #time vector arranged by index
        idxqoa=np.append(idxqoa,idxq[np.int(I1a[i])-1]) #time vector arranged by index
        qoa=np.append(qoa,querya[np.int(I1a[i])-1]) #Query function sampled by the index
        toa=np.append(toa,templatea[np.int(I2a[i])-1]) #Template function sampled by the index
    #plt.figure('Superposition')
    #plt.plot(qo)
    #plt.plot(to)

    #Plot #3: A Summary plot with the alignment function linking the two waveforms
    # definitions for the axes
    left, width = 0.12, 0.60
    bottom, height = 0.08, 0.60
    bottom_h =  0.16 + width
    left_h = left + 0.27
    rect_plot = [left_h, bottom, width, height]
    rect_x = [left_h, bottom_h, width, 0.2]
    rect_y = [left, bottom, 0.2, height]

    # start with a rectangular Figure
    fig=plt.figure('Summary', figsize=(10, 10))

    #axplot = fig.add_subplot(111)
    axplot = plt.axes(rect_plot)
    axx = plt.axes(rect_x)
    axy = plt.axes(rect_y)

    # Plot the matrix
    axplot.plot(idxqo,idxto,'k',picker=5,lw=2)
    #axplot.plot(idxqoa,idxtoa,'m',picker=5,lw=2)
    axplot.plot(idxqoh,idxtoh,'--',color='m',picker=5,lw=2)
    axplot.axis([t_id,t_fd,t_id,t_fd]) #Give same time scale as template
    axplot.grid(color='0.5')


    #Indices in P correspond to locations where S-times saturated are larger than S-times dry
    P=np.where(idxqo>=idxto) #Find the locations where time of S sat. is larger than time of S dry 
    Ph=np.where(idxqoh>=idxtoh)
    #axplot.scatter(idxqo[P],idxto[P],c='r',s=10)
    #axplot.scatter(idxqoh[Ph],idxtoh[Ph],c='r',s=10)
    #Define and plot the 1:1 line of match
    x1=np.linspace(0,np.max([timeSd,timeSs]),10)
    y1=x1
    axplot.plot(x1,y1,'g',lw=2) #plot the 1:1 line of match

    axplot.tick_params(axis='both', which='major', labelsize=18)

    # Plot time serie horizontal
    axx.plot(idxq,query,'-', color='b',lw=2)
    axx.plot(idxq,queryh,'--', color='m',lw=2) #The envelope
    #axx.plot(idxq,querya,'-',color='m',lw=2)
    axx.grid(color='0.5')
    axx.axis([t_id,t_fd,1.1*np.min(query),1.1*np.max(queryh)])
    axx.tick_params(axis='both', which='major', labelsize=18)


    # Plot time serie vertical
    axy.plot(template,idxt,'-',color='r',lw=2)
    axy.plot(templateh,idxt,'--',color='m',lw=2)
    #axy.plot(templatea,idxt,'-',color='m',lw=2)
    axy.grid(color='0.5')
    axy.axis([1.1*np.min(template),1.1*np.max(templateh),t_id,t_fd])
    axy.invert_xaxis()
    axy.tick_params(axis='both', which='major', labelsize=18)

    #Limits
    axx.set_xlim(axplot.get_xlim())
    axy.set_ylim(axplot.get_ylim())


    #A plot to pick the S-wave arrival from the alignment function
    #Pick the S-arrival time.

    #Pick the times from the graph
    axplot.set_title('click on points')

    #Global variable for storing the picked coordinates 
    COORDS=[]
    fig.canvas.mpl_connect('pick_event', onpick)

    plt.show(block=True)#plt.show(block=True)

    #END OF THE VS PICKS 

    #ESTIMATION OF THE VELOCITIES FROM THE PICKED TIMES
    ##
    #Load the caliper lengths
    Lengths=np.loadtxt('./'+well+'/'+well+'_Lengths.txt',delimiter=',')
    L,StdL=rp.Length_sample(Lengths)

    #Extract the times from the DTW picking
    tD,tS=rp.S_values(COORDS)
    np.savetxt('./'+well+'/'+well+'d/P_time_Picks_dry_'+suffix+'.out',tD) #save the picked times
    np.savetxt('./'+well+'/'+well+'s/P_time_Picks_sat_'+suffix+'.out',tS) #save the picked times
    #Save the histograms
    rp.histogram(tD,'P-arrival (dry)')
    plt.xlabel('Time ($\mu$s)')
    plt.savefig('./'+well+'/'+well+'d/P_times_hist_dry_'+suffix+'.pdf',bbox_inches='tight')
    plt.show(block=True)
    rp.histogram(tS,'P-arrival (sat.)')
    plt.xlabel('Time ($\mu$s)')
    #save the figure
    manager= plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.show()
    plt.savefig('./'+well+'/'+well+'s/P_times_hist_sat_'+suffix+'.pdf',bbox_inches='tight')
    plt.show(block=True)

    #CALCULATE THE VELOCITIES FROM THE TIMES PICKED FOR THE DRY SAMPLE
    Vd,dVd,Vd_down,Vdmc,Vd_up=rp.Velocity_P(L,np.mean(tD),StdL,np.std(tD)) #Velocity S is more general than Velocity_P
    #Save the velocities and intervals
    np.savetxt('./'+well+'/'+well+'d/Velocities_P_dry_'+suffix+'.out', [Vd,dVd,Vd_down,Vdmc.mean,Vd_up], fmt='%1.2f',delimiter=',',header='Vp (dry), StdVp (dry), V_0025, V_mean, V_0975')
    #Plot the velocities distributions
    plt.figure('Vp (dry)')
    Veld=mc.N(Vd,dVd)
    Veld.plot(label='Vp (dry)',lw=2,color='b')
    Vdmc.plot(hist=True,label='Vp (dry) (MC)',color='g')
    plt.legend()
    plt.xlabel('Velocity (m/s)')
    #Save the figure
    manager= plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.show()
    plt.savefig('./'+well+'/'+well+'d/Vp_hist_dry_'+suffix+'.pdf',bbox_inches='tight')
    plt.show(block=True)#plt.show(block=True)

    #CALCULATE THE VELOCITIES FROM THE TIMES PICKED FOR THE SATURATED SAMPLE
    Vp,dVp,Vp_down,Vpmc,Vp_up=rp.Velocity_P(L,np.mean(tS),StdL,np.std(tS)) #Velocity S is more general than Velocity_P
    #Save the velocities and intervals
    np.savetxt('./'+well+'/'+well+'s/Velocities_P_sat_'+suffix+'.out', [Vp,dVp,Vp_down,Vpmc.mean,Vp_up], fmt='%1.2f',delimiter=',',header='Vp (sat), StdVp (sat), V_0025, V_mean, V_0975')
    #Plot the velocities distributions
    plt.figure('Vp (sat.)')
    Vels=mc.N(Vp,dVp)
    Vels.plot(label='Vp (sat.)',lw=2,color='b')
    Vpmc.plot(hist=True,label='Vp (sat.) (MC)',color='g')
    plt.legend()
    plt.xlabel('Velocity (m/s)')
    #save the figure
    manager= plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.show()
    plt.savefig('./'+well+'/'+well+'s/Vp_hist_sat_'+suffix+'.pdf',bbox_inches='tight')
    plt.show(block=True)#plt.show(block=True)

    #Plot #5: A Summary plot with the picked intervals
    # definitions for the axes
    left, width = 0.12, 0.60
    bottom, height = 0.08, 0.60
    bottom_h =  0.16 + width
    left_h = left + 0.27
    rect_plot = [left_h, bottom, width, height]
    rect_x = [left_h, bottom_h, width, 0.2]
    rect_y = [left, bottom, 0.2, height]

    # start with a rectangular Figure
    fig=plt.figure('Summary', figsize=(10, 10))

    #axplot = fig.add_subplot(111)
    axplot = plt.axes(rect_plot)
    axx = plt.axes(rect_x)
    axy = plt.axes(rect_y)

    # Plot the matrix
    axplot.plot(idxqo,idxto,'k',lw=2)
    #axplot.plot(idxqoa,idxtoa,'m',lw=2)
    axplot.plot(idxqoh,idxtoh,'--',color='m',lw=2)
    axplot.axvspan(np.min(tS),np.max(tS),alpha=0.5,color='blue') #plot range of picked times
    axplot.axvline(np.mean(tS),ymin=t_id,ymax=t_fd,linewidth=2, color='b') #plot mean of the picked times 
    axplot.axhspan(np.min(tD),np.max(tD),alpha=0.5,color='red') #plot range of picked times
    axplot.axhline(np.mean(tD),xmin=t_id,xmax=t_fd,linewidth=2, color='r') #plot mean of the picked times 
    axplot.axis([t_id,t_fd,t_id,t_fd]) #Give same time scale as template
    axplot.grid(color='0.5')


    #Indices in P correspond to locations where S-times saturated are larger than S-times dry
    P=np.where(idxqo>=idxto) #Find the locations where time of S sat. is larger than time of S dry 
    Ph=np.where(idxqoh>=idxtoh)
    #axplot.scatter(idxqo[P],idxto[P],c='r',s=10)
    #axplot.scatter(idxqoh[Ph],idxtoh[Ph],c='r',s=10)
    #Define and plot the 1:1 line of match
    x1=np.linspace(0,np.max([timeSd,timeSs]),10)
    y1=x1
    axplot.plot(x1,y1,'g') #plot the 1:1 line of match

    axplot.tick_params(axis='both', which='major', labelsize=18)

    # Plot time serie horizontal
    axx.plot(idxq,query,'-', color='b',label='P (sat.)',lw=2)
    axx.plot(idxq,queryh,'--', color='m',lw=2) #The envelope
    axx.axvspan(np.min(tS),np.max(tS),alpha=0.5,color='blue') #plot range of picked times
    axx.axvline(np.mean(tS),ymin=-10,ymax=10,linewidth=2, color='b') #plot mean of the picked times 
    axx.grid(color='0.5')
    axx.axis([t_id,t_fd,1.1*np.min(query),1.1*np.max(queryh)])
    axx.tick_params(axis='both', which='major', labelsize=18)


    # Plot time serie vertical
    axy.plot(template,idxt,'-',color='r',label='P (dry)',lw=2)
    axy.plot(templateh,idxt,'--',color='m',lw=2)
    axy.axhspan(np.min(tD),np.max(tD),alpha=0.5,color='red') #plot range of picked times
    axy.axhline(np.mean(tD),xmin=-10,xmax=10,linewidth=2, color='r') #plot mean of the picked times 
    axy.grid(color='0.5')
    axy.axis([1.1*np.min(template),1.1*np.max(templateh),t_id,t_fd])
    axy.invert_xaxis()
    axy.tick_params(axis='both', which='major', labelsize=18)

    #Limits
    axx.set_xlim(axplot.get_xlim())
    axy.set_ylim(axplot.get_ylim())

    #save the figure
    plt.savefig('./'+well+'/'+well+'d/Summary_Match_Pwaves_pick_'+suffix+'.pdf',bbox_inches='tight')
    plt.savefig('./'+well+'/'+well+'s/Summary_Match_Pwaves_pick_'+suffix+'.pdf',bbox_inches='tight')
    plt.show(block=True)#plt.show(block=True)


    #Read the density files
    A=np.loadtxt('./'+well+'/Dry_density.out')
    B=np.loadtxt('./'+well+'/Sat_density.out')
    C=np.loadtxt('./'+well+'/'+well+'d/Shear_dry.out')
    D=np.loadtxt('./'+well+'/'+well+'s/Shear_sat.out')
    rhod=A[0]
    drhod=A[1]
    rhos=B[0]
    drhos=B[1]
    mud=C[0]
    dmud=C[1]
    mus=D[0]
    dmus=D[1]
    #Create the distributions for the MC simulations
    RHOD=mc.N(rhod,drhod)
    RHOS=mc.N(rhos,drhos)

    #Create distributions for Mud and Mus
    MUD=mc.N(mud,dmud)
    MUS=mc.N(mus,dmus)

    #Bulk modulus
    KD=1e-6*RHOD*Veld**2-(4/3)*MUD
    KS=1e-6*RHOS*Vels**2-(4/3)*MUS

    #Young Modulus 
    PD=(3*KD-2*MUD)/(2*(3*KD+MUD))
    PS=(3*KS-2*MUS)/(2*(3*KS+MUS))
    #Poisson's ratio
    ED=9*KD*MUD/(3*KD+MUD)
    ES=9*KS*MUS/(3*KS+MUS)

    ##Mu
    ##Read the density files
    #A=np.loadtxt('./'+well+'/Dry_density.out')
    #B=np.loadtxt('./'+well+'/Sat_density.out')
    #rhod=A[0]
    #drhod=A[1]
    #rhos=B[0]
    #drhos=B[1]
    ##Create the distributions for the MC simulations
    #RHOD=mc.N(rhod,drhod)
    #RHOS=mc.N(rhos,drhos)
    #
    #MUD=1e-6*RHOD*Veld**2 #Values in GPa
    #MUS=1e-6*RHOD*Vels**2 #Values in GPa
    #
    #plt.figure('Rho (sat.)')
    #RHOS.plot(label='RHO (sat.)',lw=2,color='b')
    #plt.savefig('./'+well+'/'+well+'s/Density_sat_'+suffix+'.pdf',bbox_inches='tight')
    #
    #plt.figure('Rho (dry)')
    #RHOD.plot(label='RHO (dry)',lw=2,color='b')
    #plt.savefig('./'+well+'/'+well+'d/Density_dry_'+suffix+'.pdf',bbox_inches='tight')

    plt.figure('K (sat.)')
    KS.plot(label='$K$ (sat.)',lw=2,color='b')
    KS.plot(hist=True,label='$K$ (sat.) (MC)',color='g')
    plt.legend()
    plt.xlabel('$K$ (GPa)')
    plt.savefig('./'+well+'/'+well+'s/Bulk_sat_'+suffix+'.pdf',bbox_inches='tight')

    plt.figure('Bulk (dry)')
    KD.plot(label='$K$ (dry)',lw=2,color='b')
    KD.plot(hist=True,label='$K$ (dry) (MC)',color='g')
    plt.legend()
    plt.xlabel('$K$ (GPa)')
    plt.savefig('./'+well+'/'+well+'d/Bulk_dry_'+suffix+'.pdf',bbox_inches='tight')

    #Poisson's ratio figures
    plt.figure('Poisson (sat.)')
    PS.plot(label='$v$ (sat.)',lw=2,color='b')
    PS.plot(hist=True,label='$v$ (sat.) (MC)',color='g')
    plt.legend()
    plt.xlabel('$v$')
    plt.savefig('./'+well+'/'+well+'s/Poisson_sat_'+suffix+'.pdf',bbox_inches='tight')


    plt.figure('Poisson (dry.)')
    PD.plot(label='$v$ (dry.)',lw=2,color='b')
    PD.plot(hist=True,label='$v$ (dry.) (MC)',color='g')
    plt.legend()
    plt.xlabel('$v$')
    plt.savefig('./'+well+'/'+well+'d/Poisson_dry_'+suffix+'.pdf',bbox_inches='tight')

    #Young's modulus figures
    plt.figure('E (sat.)')
    ES.plot(label='$E$ (sat.)',lw=2,color='b')
    ES.plot(hist=True,label='$E$ (sat.) (MC)',color='g')
    plt.legend()
    plt.xlabel('$E$ (GPa)')
    plt.savefig('./'+well+'/'+well+'s/Young_sat_'+suffix+'.pdf',bbox_inches='tight')


    plt.figure('E (dry.)')
    ED.plot(label='$E$ (dry.)',lw=2,color='b')
    ED.plot(hist=True,label='$E$ (dry.) (MC)',color='g')
    plt.legend()
    plt.xlabel('$E$ (GPa)')
    plt.savefig('./'+well+'/'+well+'d/Young_dry_'+suffix+'.pdf',bbox_inches='tight')



    ## normalize by standard deviation (not necessary, but makes it easier
    ## to compare with plot on Interactive Wavelet page, at
    ## "http://paos.colorado.edu/research/wavelets/plot/"
    #
    #sst = templateh
    #ssq = queryh
    #dt = 0.01
    #timet = idxt
    #timeq = idxq
    ##xlim = ([1870, 2000])  # plotting range
    #pad = 1  # pad the time series with zeroes (recommended)
    #dj = 0.25  # this will do 4 sub-octaves per octave
    #s0 = 2 * dt  # this says start at a scale of 6 months
    #j1 = 7 / dj  # this says do 7 powers-of-two with dj sub-octaves each
    #lag1 = 0.72  # lag-1 autocorrelation for red noise background
    #mother = 'MORLET'#'PAUL'#'MORLET'#'PAUL'
    #
    ## Wavelet transform:
    #wavet, periodt, scalet, coit = wavelet(sst, dt, pad, dj, s0, j1, mother)
    #waveq, periodq, scaleq, coiq = wavelet(ssq, dt, pad, dj, s0, j1, mother)
    #powert = (np.abs(wavet)) ** 2  # compute wavelet power spectrum
    #phaset=(np.angle(wavet))
    #powerq= (np.abs(waveq)) ** 2  # compute wavelet power spectrum
    #phaseq=(np.angle(waveq))
    #
    ##plt.figure('Wavelet')
    ##plt.imshow((powert),extent=[np.min(timet),np.max(timet),0,1/(2*0.01)],aspect='auto',cmap='gist_stern')
    ##plt.colorbar()
    ##plt.show()
    ##plt.figure('Wavelet2')
    ##plt.imshow(phaset,extent=[np.min(timet),np.max(timet),0,1/(2*0.01)],aspect='auto',cmap='gist_stern')
    ##plt.colorbar()
    #
    ## start with a rectangular Figure
    #plt.figure('Summary Wavelet', figsize=(8, 8))
    #
    #axplot = plt.axes(rect_plot)
    #axx = plt.axes(rect_x)
    #axy = plt.axes(rect_y)
    #
    ## Plot the matrix
    ##axplot.pcolor(acc.T,cmap=cm.gray)
    #axplot.plot(idxqo,idxto)
    #axplot.plot(idxqoh,idxtoh,'m')
    #axplot.axvspan(np.min(tS),np.max(tS),alpha=0.5,color='red') #plot range of picked times
    #axplot.axhspan(np.min(tD),np.max(tD),alpha=0.5,color='blue') #plot range of picked times
    #axplot.grid(color='0.5')
    #
    ##Indices in P correspond to locations where S-times saturated are larger than S-times dry
    ##P=np.where(idxqo>=idxto) #Find the locations where time of S sat. is larger than time of S dry 
    ##axplot.scatter(idxqo[P],idxto[P],c='r',s=10)
    ##Define and plot the 1:1 line of match
    #x1=np.linspace(0,np.max([timeSd,timeSs]),10)
    #y1=x1
    #axplot.plot(x1,y1,'g') #plot the 1:1 line of match
    #
    #axplot.set_xlim((np.min(idxqo), np.max(idxqoh)))
    #axplot.set_ylim((np.min(idxto), np.max(idxtoh)))
    #axplot.tick_params(axis='both', which='major', labelsize=18)
    #
    ## Plot time serie horizontal
    #axx.imshow(powerq,extent=[np.min(timeq),np.max(timeq),0,1/(2*0.01)],aspect='auto',cmap='gist_stern')
    ##axx.plot(idxq,queryh,'-', color='m')
    #axx.axvspan(np.min(tS),np.max(tS),alpha=0.5,color='red') #plot range of picked times
    #axx.axvline(np.mean(tS),ymin=-10,ymax=10,linewidth=2, color='r') #plot mean of the picked times 
    #axx.axis('tight')
    #axx.grid(color='0.5')
    #axx.tick_params(axis='both', which='major', labelsize=18)
    #xloc = plt.MaxNLocator(4)
    #x2Formatter = FormatStrFormatter('%1.2f')
    #axx.yaxis.set_major_locator(xloc)
    #axx.yaxis.set_major_formatter(x2Formatter)
    #
    ## Plot time serie vertical
    #axy.imshow(np.fliplr(np.flipud(powert.T)),extent=[0,1/(2*0.01),np.min(timet),np.max(timet)],aspect='auto',cmap='gist_stern')#plot(template,idxt,'-',color='b')
    ##axy.plot(templateh,idxt,'-',color='m')
    #axy.axhspan(np.min(tD),np.max(tD),alpha=0.5,color='blue') #plot range of picked times
    #axy.axhline(np.mean(tD),xmin=-10,xmax=10,linewidth=2, color='b') #plot mean of the picked times 
    #axy.axis('tight')
    #axy.grid(color='0.5')
    #axy.invert_xaxis()
    #yloc = plt.MaxNLocator(4)
    #xFormatter = FormatStrFormatter('%1.2f')
    #axy.xaxis.set_major_locator(yloc)
    #axy.xaxis.set_major_formatter(xFormatter)
    #axy.tick_params(axis='both', which='major', labelsize=18)
    #
    ##Limits
    #axx.set_xlim(axplot.get_xlim())
    #axy.set_ylim(axplot.get_ylim())
    #
    #
    ##save the figure
    ##plt.savefig('./'+well+'/'+well+'d/Summary_Match_Pwaves_pick_'+suffix+'.pdf',bbox_inches='tight')
    ##plt.savefig('./'+well+'/'+well+'s/Summary_Match_Pwaves_pick_'+suffix+'.pdf',bbox_inches='tight')
    ##plt.show(block=True)#plt.show(block=True)
    #
    #
    ## start with a rectangular Figure
    #plt.figure('Summary Wavelet 2', figsize=(8, 8))
    #
    #axplot = plt.axes(rect_plot)
    #axx = plt.axes(rect_x)
    #axy = plt.axes(rect_y)
    #
    ## Plot the matrix
    ##axplot.pcolor(acc.T,cmap=cm.gray)
    #axplot.plot(idxqo,idxto)
    #axplot.plot(idxqoh,idxtoh,'m')
    #axplot.axvspan(np.min(tS),np.max(tS),alpha=0.5,color='red') #plot range of picked times
    #axplot.axhspan(np.min(tD),np.max(tD),alpha=0.5,color='blue') #plot range of picked times
    #axplot.grid(color='0.5')
    #
    ##Indices in P correspond to locations where S-times saturated are larger than S-times dry
    ##P=np.where(idxqo>=idxto) #Find the locations where time of S sat. is larger than time of S dry 
    ##axplot.scatter(idxqo[P],idxto[P],c='r',s=10)
    ##Define and plot the 1:1 line of match
    #x1=np.linspace(0,np.max([timeSd,timeSs]),10)
    #y1=x1
    #axplot.plot(x1,y1,'g') #plot the 1:1 line of match
    #
    #axplot.set_xlim((np.min(idxqo), np.max(idxqoh)))
    #axplot.set_ylim((np.min(idxto), np.max(idxtoh)))
    #axplot.tick_params(axis='both', which='major', labelsize=18)
    #
    ## Plot time serie horizontal
    #axx.imshow(phaseq,extent=[np.min(timeq),np.max(timeq),0,1/(2*0.01)],aspect='auto',cmap='gist_stern')
    ##axx.plot(idxq,queryh,'-', color='m')
    #axx.axvspan(np.min(tS),np.max(tS),alpha=0.5,color='red') #plot range of picked times
    #axx.axvline(np.mean(tS),ymin=-10,ymax=10,linewidth=2, color='r') #plot mean of the picked times 
    #axx.axis('tight')
    #axx.grid(color='0.5')
    #axx.tick_params(axis='both', which='major', labelsize=18)
    #xloc = plt.MaxNLocator(4)
    #x2Formatter = FormatStrFormatter('%1.2f')
    #axx.yaxis.set_major_locator(xloc)
    #axx.yaxis.set_major_formatter(x2Formatter)
    #
    ## Plot time serie vertical
    #axy.imshow(np.fliplr(np.flipud(phaset.T)),extent=[0,1/(2*0.01),np.min(timet),np.max(timet)],aspect='auto',cmap='gist_stern')#plot(template,idxt,'-',color='b')
    ##axy.plot(templateh,idxt,'-',color='m')
    #axy.axhspan(np.min(tD),np.max(tD),alpha=0.5,color='blue') #plot range of picked times
    #axy.axhline(np.mean(tD),xmin=-10,xmax=10,linewidth=2, color='b') #plot mean of the picked times 
    #axy.axis('tight')
    #axy.grid(color='0.5')
    #axy.invert_xaxis()
    #yloc = plt.MaxNLocator(4)
    #xFormatter = FormatStrFormatter('%1.2f')
    #axy.xaxis.set_major_locator(yloc)
    #axy.xaxis.set_major_formatter(xFormatter)
    #axy.tick_params(axis='both', which='major', labelsize=18)
    #
    ##Limits
    #axx.set_xlim(axplot.get_xlim())
    #axy.set_ylim(axplot.get_ylim())
    #
    #
    #
    #Font size for plots
    font= {'size' : 32}
    plt.rc('font',**font)
    plt.figure('Waveform Plot')
    plt.plot(idxq,query,color='b',label='P-waveform',lw=3)
    plt.plot(idxq,queryh,'--',color='m',label='Envelope',lw=3)
    plt.plot(idxq,-queryh,'--',color='m',lw=3)
    plt.axis([np.min(idxq),np.max(idxq),-1.05*np.max(queryh),1.05*np.max(queryh)])
    plt.grid(color='0.5')
    plt.legend()
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Normalized Amplitude (A.U.)')
    manager= plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.show()
    plt.savefig('./'+well+'/Waveform_P_'+well+'.pdf',bbox_inches='tight')
