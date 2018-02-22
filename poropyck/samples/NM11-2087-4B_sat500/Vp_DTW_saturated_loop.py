# -*- coding: utf-8 -*-
"""
Created on Sat Apr 08 14:27:30 2017

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

# Set up our R namespaces
R = rpy2.robjects.r
DTW = importr('dtw')

##Choose the folder with waveforms
#well=input("Type name of sample (e.g. 'NM11_2087-4B_sat3000'): \n")

#Read the waveforms 
#REORDER THE FILENAMES SO THE FIRST INDEX IS THE HIGHEST PRESSURE
files = sorted(glob.glob('*.csv'),key=lambda f: list(filter(str.isdigit, f)),reverse=True)#

Pc=[]
phases=[]
states=[]
cycles=[]    
for i in range(0,len(files)):
    name,Pconf,phase,state,cycle=rp.extract_parameters(files[i])
    Pc=np.append(Pc,Pconf)
    phases=np.append(phases,phase)
    states=np.append(states,state)
    cycles=np.append(cycles,cycle)

#CHOOSE THE PHASE TO PICK
indices=np.where((phases=='PP')& (cycles=='d2')& (states=='sat500'))[0]
indices_sorted=[x for _,x in sorted(zip(Pc[indices],indices))] #Sort according to pressures   

#PICK THE S-TIME FOR THE HIGHEST PRESSURE EXPERIMENT TO BACK-TRACK IT TO LOWER ONES
timeSd,Sd=rp.load_csv(files[indices[0]],0)
print('Pick the P-time for this pressure')
#Choose the two extremes to compare
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('click on points (Pc = '+np.str(Pc[indices[0]])+')')
coords=[]
ax.plot(timeSd,Sd,'-',picker=5),plt.xlabel('Time (microsecs)')
plt.grid('on')
fig.canvas.mpl_connect('pick_event', onpick)
plt.show(block=True)

tH=rp.P_values(coords) #PICK THE TIMES FOR THE HIGHEST PRESSURE
mtH,dtH=rp.histogram(tH,'Times for highest pressure Pc = '+np.str(Pc[indices[0]])+' psi')
M=np.array([mtH, dtH])
np.savetxt(''+phases[indices[0]]+'_time_std_'+states[indices[0]]+'_'+np.str(np.int(Pc[indices[0]]))+'.out',M) #save the picked times
plt.show(block=True)


#LOOP THROUGH ALL THE FILES OF THE FOLDER
for k in range(0,len(indices)-1): 

    #READ THE HIGHEST PRESSURE WAVEFORM FIRST
    timeSd,Sd=rp.load_csv(files[indices[k]],0)
    
    print('Choose the beginning and end that you want to compare')
    #Choose the two extremes to compare
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('click on points (Pc = '+np.str(Pc[indices[k]])+')')
    coords=[]
    ax.plot(timeSd,Sd,'-',picker=5),plt.xlabel('Time (microsecs)')
    plt.grid('on')
    fig.canvas.mpl_connect('pick_event', onpick)
    plt.show(block=True)
    #Pick the times
    t_id=np.min(coords[0][0])
    t_fd=np.min(coords[1][0])
    
    #Indices to crop the dry sample time vector
    i_d=np.min(np.where(timeSd>=t_id))
    f_d=np.min(np.where(timeSd>=t_fd))
    
    #Read the next lower pressure waveform
    timeSs,Ss=rp.load_csv(files[indices[k+1]],0)
    
    print('Choose the beginning and end that you want to compare')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('click on points (Pc = '+np.str(Pc[indices[k+1]])+')')
    coords=[]
    ax.plot(timeSs,Ss,'-',picker=5),plt.xlabel('Time (microsecs)')
    plt.grid('on')
    fig.canvas.mpl_connect('pick_event', onpick)
    plt.show(block=True)
    
    t_is=np.min(coords[0][0])
    t_fs=np.min(coords[1][0])
    
    #Indices to crop for the saturated sample time vector
    i_s=np.min(np.where(timeSs>=t_is))
    f_s=np.min(np.where(timeSs>=t_fs))
    
    #DYNAMIC TIME WARPING ALGORITHM BEGINS HERE!!!!
    idxt=timeSd[i_d:f_d]
    idxq=timeSs[i_s:f_s]
    
    #2.1) CHOOSE WHAT INFORMATION TO COMPARE (EITHER WAVEFORMS OR ENVELOPES)
    template=Sd[i_d:f_d]/np.max(np.abs(Sd[i_d:f_d]))
    query=Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s]))
    templateh=np.abs(hilbert(Sd[i_d:f_d]/np.max(np.abs(Sd[i_d:f_d]))))
    queryh=np.abs(hilbert(Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s]))))
    
    
    # Calculate the alignment vector and corresponding distance (DTW)
    alignment = R.dtw(query, template, keep=True)
    alignmenth = R.dtw(queryh, templateh, keep=True)
    dist = alignment.rx('distance')[0][0] #The lower the distance of alignment, the better the match
    disth = alignmenth.rx('distance')[0][0]
    
    print('Distance of the DTW algorithm (waveform): ',dist)
    print('Distance of the DTW algorithm (envelope): ',disth)
    
    #Indices of alignment
    I1=alignment.rx('index1')[0]
    I2=alignment.rx('index2')[0]
    I1h=alignmenth.rx('index1')[0]
    I2h=alignmenth.rx('index2')[0]
    
    print('Close the figures to continue running the code...')
    #Plots
    #Plot #1
    plt.figure('Points of match')
    plt.plot(idxt,template,label='Pd',c='b')
    plt.plot(idxq,query,label='Ps',c='g')
    plt.axis('tight')
    plt.legend()
    plt.xlabel('time (microsecs)')
    #Plot the matching points
    for i in np.arange(0,len(I1)-1,20):
        plt.plot([idxt[np.int(I2[i]-1)],idxq[np.int(I1[i]-1)]],[template[np.int(I2[i]-1)],query[np.int(I1[i]-1)]],'r-',lw=0.5)
    
    plt.show(block=True)
    
    #Plot #2 (mostly unnecessary): Compare a function over the other
    #plt.figure('One over the other')
    qo=[]
    qoh=[]
    to=[]
    toh=[]
    idxto=[]
    idxtoh=[]
    idxqo=[]
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
    axplot.plot(idxqo,idxto,'k',picker=5)
    axplot.plot(idxqoh,idxtoh,'m',picker=5)
    axplot.axis([t_id,t_fd,t_id,t_fd]) #Give same time scale as template
    axplot.grid('on')
    
    
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
    axx.plot(idxq,query,'-', color='b')
    axx.plot(idxq,queryh,'-', color='m') #The envelope
    axx.grid('on')
    axx.axis([t_id,t_fd,1.1*np.min(query),1.1*np.max(queryh)])
    axx.tick_params(axis='both', which='major', labelsize=18)
    
    
    # Plot time serie vertical
    axy.plot(template,idxt,'-',color='r')
    axy.plot(templateh,idxt,'-',color='m')
    axy.grid('on')
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
    coords=[]
    
    #Load the last picks and plot the picking band
    mean,std_error=np.loadtxt(''+phases[indices[k]]+'_time_std_'+states[indices[k]]+'_'+np.str(np.int(Pc[indices[k]]))+'.out')    
    
    axplot.axhline(mean+2*std_error, color='red', lw=2, alpha=0.5)
    axplot.axhline(mean-2*std_error, color='red', lw=2, alpha=0.5)  
    axplot.axis([t_id,t_fd,t_id,t_fd]) #Give same time scale as template    
    #plt.axis('equal'),plt.grid('on')
    
    fig.canvas.mpl_connect('pick_event', onpick)
    
    plt.show(block=True)
    
    #Extract the times from the DTW picking
    tD,tS=rp.S_values(coords)
    np.savetxt(''+phases[indices[k]]+'_time_Picks_'+states[indices[k]]+'_'+np.str(np.int(Pc[indices[k]]))+'.out',tD) #save the picked times
    np.savetxt(''+phases[indices[k+1]]+'_time_Picks_'+states[indices[k]]+'_'+np.str(np.int(Pc[indices[k+1]]))+'.out',tS) #save the picked times
    #Save the mean and standard error of the picking for the next iteration  
    M=np.array([np.mean(tS),np.std(tS)])
    np.savetxt(''+phases[indices[k+1]]+'_time_std_'+states[indices[k+1]]+'_'+np.str(np.int(Pc[indices[k+1]]))+'.out',M) #save the picked times    
    #Save the histograms
    rp.histogram(tD,'P-arrival ('+np.str(np.int(Pc[indices[k]]))+')')
    plt.xlabel('Time (microsecs)')
    plt.savefig(''+phases[indices[k]]+'_times_hist_'+states[indices[k]]+'_'+np.str(np.int(Pc[indices[k]]))+'.pdf',bbox_inches='tight')
    plt.show(block=True)
    rp.histogram(tS,'P-arrival ('+np.str(np.int(Pc[indices[k+1]]))+')')
    plt.xlabel('Time (microsecs)')
    plt.savefig(''+phases[indices[k+1]]+'_times_hist_'+states[indices[k+1]]+'_'+np.str(np.int(Pc[indices[k+1]]))+'.pdf',bbox_inches='tight')
    plt.show(block=True)

# Set default font size for plots:
font = {'size'   : 24}
plt.rc('font',**font)
Pp=500 #psi

    
fig = plt.figure('tp vs. Pressure')
ax = fig.add_subplot(111)
dP=50*np.ones(len(indices))  
Vs=[]
dVs=[]
for i in range(0,len(indices)):
    mean,std_error=np.loadtxt(''+phases[indices[i]]+'_time_std_'+states[indices[0]]+'_'+np.str(np.int(Pc[indices[i]]))+'.out')
    Vs=np.append(Vs,mean)    
    dVs=np.append(dVs,std_error)
    ax.scatter(Pc[indices[i]]-Pp,mean)
ax.errorbar(Pc[indices]-Pp, Vs, yerr=[2*dVs, 2*dVs], xerr=[2*dP, 2*dP], fmt='o')
ax.plot(Pc[indices]-Pp, Vs,'b')
plt.xlabel('Differential pressure (psi)'),plt.ylabel('P-time ($\mu s$)')
plt.grid('on')
plt.savefig('tp_'+phases[indices[i]]+'_picks_'+np.str(Pp)+'.pdf')
plt.show(block=True)

well='NM11_2087_4B'
#Load the caliper lengths
Lengths=np.loadtxt(''+well+'_Lengths.txt',delimiter=',')
L,StdL=rp.Length_sample(Lengths)


VP=[]
dVP=[]
for k in range(0,len(indices)):
    #Load the picked times file
    #Load the last picks and plot the picking band
    tDm,dtD=np.loadtxt(''+phases[indices[k]]+'_time_std_'+states[indices[k]]+'_'+np.str(np.int(Pc[indices[k]]))+'.out')    
    #CALCULATE THE VELOCITIES FROM THE TIMES PICKED FOR THE DRY SAMPLE
    Vd,dVd,Vd_down,Vdmc,Vd_up=rp.Velocity_S(L,tDm-3.5,StdL,dtD) #Velocity S is more general than Velocity_P
    np.savetxt(''+phases[indices[k]]+'_VP_std_'+states[indices[k]]+'_'+np.str(np.int(Pc[indices[k]]))+'.out', [Vd,dVd,Vd_down,Vdmc.mean,Vd_up], fmt='%1.2f',delimiter=',',header='Vp (dry), StdVp (dry), V_0025, V_mean, V_0975')
    VP=np.append(VP,Vd)
    dVP=np.append(dVP,dVd)

fig = plt.figure('Vp vs. Pressure')
ax = fig.add_subplot(111)
ax.errorbar(Pc[indices]-Pp, VP, yerr=[2*dVP, 2*dVP], xerr=[2*dP, 2*dP], fmt='o')
ax.plot(Pc[indices]-Pp, VP,'b')
plt.xlabel('Differential pressure (psi)'),plt.ylabel('$V_P$ (m/s)')
plt.grid('on')
plt.savefig('Vp_'+phases[indices[i]]+'_picks_'+np.str(Pp)+'.pdf')
plt.show(block=True)

#    #Save the velocities and intervals
#    #np.savetxt('./'+well+'/'+well+'d/Velocities_P_dry_'+suffix+'.out', [Vd,dVd,Vd_down,Vdmc.mean,Vd_up], fmt='%1.2f',delimiter=',',header='Vp (dry), StdVp (dry), V_0025, V_mean, V_0975')
#    #Plot the velocities distributions
#    plt.figure('Vp (dry)')
#    Veld=mc.N(Vd,dVd)
#    Veld.plot(label='Vp (dry)',lw=2,color='b')
#    Vdmc.plot(hist=True,label='Vp (dry) (MC)',color='g')
#    plt.legend()
#    plt.xlabel('Velocity (m/s)')
#    #Save the figure
#    manager= plt.get_current_fig_manager()
#    manager.window.showMaximized()
#    plt.show()
#    #plt.savefig('./'+well+'/'+well+'d/Vp_hist_dry_'+suffix+'.pdf',bbox_inches='tight')
#    plt.show(block=True)#plt.show(block=True)




##EDITS
##Plot #3: A Summary plot with the alignment function linking the two waveforms
## definitions for the axes
#left, width = 0.12, 0.60
#bottom, height = 0.08, 0.60
#bottom_h =  0.16 + width
#left_h = left + 0.27
#rect_plot = [left_h, bottom, width, height]
#rect_x = [left_h, bottom_h, width, 0.2]
#rect_y = [left, bottom, 0.2, height]
#
## start with a rectangular Figure
#fig=plt.figure('Summary', figsize=(10, 10))
#
##axplot = fig.add_subplot(111)
#axplot = plt.axes(rect_plot)
#axx = plt.axes(rect_x)
#axy = plt.axes(rect_y)
#
## Plot the matrix
#axplot.plot(idxqo,idxto,'k',picker=5)
#axplot.plot(idxqoh,idxtoh,'m',picker=5)
#axplot.axis([t_id,t_fd,t_id,t_fd]) #Give same time scale as template
#axplot.grid('on')
#
#
##Indices in P correspond to locations where S-times saturated are larger than S-times dry
#P=np.where(idxqo>=idxto) #Find the locations where time of S sat. is larger than time of S dry 
#Ph=np.where(idxqoh>=idxtoh)
##axplot.scatter(idxqo[P],idxto[P],c='r',s=10)
##axplot.scatter(idxqoh[Ph],idxtoh[Ph],c='r',s=10)
##Define and plot the 1:1 line of match
#x1=np.linspace(0,np.max([timeSd,timeSs]),10)
#y1=x1
#axplot.plot(x1,y1,'g') #plot the 1:1 line of match
#
#axplot.tick_params(axis='both', which='major', labelsize=18)
#
## Plot time serie horizontal
#axx.plot(idxq,query,'-', color='b')
#axx.plot(idxq,queryh,'-', color='m') #The envelope
#axx.grid('on')
#axx.axis([t_id,t_fd,1.1*np.min(query),1.1*np.max(queryh)])
#axx.tick_params(axis='both', which='major', labelsize=18)
#
#
## Plot time serie vertical
#axy.plot(template,idxt,'-',color='r')
#axy.plot(templateh,idxt,'-',color='m')
#axy.grid('on')
#axy.axis([1.1*np.min(template),1.1*np.max(templateh),t_id,t_fd])
#axy.invert_xaxis()
#axy.tick_params(axis='both', which='major', labelsize=18)
#
##Limits
#axx.set_xlim(axplot.get_xlim())
#axy.set_ylim(axplot.get_ylim())
#
#
##A plot to pick the S-wave arrival from the alignment function
##Pick the S-arrival time.
#
##Pick the times from the graph
#axplot.set_title('click on points')
#
##Global variable for storing the picked coordinates 
#coords=[]