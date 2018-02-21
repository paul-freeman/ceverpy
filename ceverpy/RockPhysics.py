# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:47:25 2017

@author: edur409
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import mcerp3 as mc
import mcerp3.umath as umath 
from scipy.stats import norm as normal
#import visa
import glob

def Voigt(M,fm,phi,sat,module):
    #Modules in Pa, fractions and porosities in percentages
    #sat is either 0 (dry) or 1 (wet)
    #Constants (MKS)
    if module=='Bulk':
        Mwater=2.2 #GPa
        Mair=0.000131 #GPa
    elif module=='Shear':
        Mwater=0
        Mair=0
    rhow=1000 
    phi=phi*0.01
    #Rescaling of fraction of rock
    if np.sum(fm)!=100:
        fm=fm*100/(np.sum(fm))
    fm=fm*0.01
    if sat==0:
        if Mair!=0:
            Meff=[phi*Mair]
        else:
            Meff=[]
        for i in range(0,len(M)):
            Meff=np.append(Meff,(1-phi)*fm[i]*M[i])
        Mveff=np.sum(Meff) #Voigt (upper bound)
    elif sat==1:
        if Mwater!=0:
            Meff=[phi*Mwater]
        else:
            Meff=[]
        for i in range(0,len(M)):
            Meff=np.append(Meff,(1-phi)*fm[i]*M[i])   
        Mveff=np.sum(Meff) #Voigt (upper bound)  
    else:
        print("type 0: for dry\n type 1: for saturated")
        return
    return Mveff
    
def Reuss(M,fm,phi,sat,module):
    #Modules in Pa, fractions and porosities in percentages
    #sat is either 0 (dry) or 1 (wet)
    #Constants (MKS)
    if module=='Bulk':
        Mwater=2.2 #GPa
        Mair=1.31e-4 #GPa
    elif module=='Shear':
        Mwater=0
        Mair=0
    rhow=1000 
    phi=phi*0.01
    #Rescaling of fraction of rock
    if np.sum(fm)!=100:
        fm=fm*100/(np.sum(fm))
    fm=fm*0.01
    if sat==0:
        if Mair!=0:
            Meff=[phi*(Mair**-1)]
        elif phi==1.0:
            Meff=np.inf
        else:
            Meff=[]
        for i in range(0,len(M)):
            Meff=np.append(Meff,(1-phi)*fm[i]*(M[i]**-1))
        Mreff=np.sum(Meff)**-1 #Reuss (lower bound)
    elif sat==1:
        if Mair!=0:
            Meff=[phi*(Mwater**-1)]
        elif phi==1.0:
            Meff=np.inf
        else:
            Meff=[]
        for i in range(0,len(M)):
            Meff=np.append(Meff,(1-phi)*fm[i]*(M[i]**-1))   
        Mreff=np.sum(Meff)**-1 #Reuss (lower bound)   
    else:
        print("type 0: for dry\n type 1: for saturated")
        return
    return Mreff
    
def V_p(K,Mu,Den):
    Vp=umath.sqrt((K*1e9+(4./3)*Mu*1e9)/(Den*1e3))  
    return Vp

def V_s(Mu,Den):
    Vs=umath.sqrt((Mu*1e9)/(Den*1e3))  
    return Vs
    
def Rock_den(den,fm,phi,sat):
    #sat= percentage of saturation 
    rhow=1.0
    rhoa=0.001225
    phi=phi*0.01
    sat=sat*0.01
        #Rescaling of fraction of rock
    if np.sum(fm)!=100:
        fm=fm*100/(np.sum(fm))
    fm=fm*0.01
    Den=np.sum(den*fm*(1.0-phi))+(1.0-sat)*phi*rhoa+sat*phi*rhow
    return Den,fm
    
def Wood(M,fm,phi,sat,module):
    #F in percentages
    #sat is either 0 (dry) or 1 (wet)
    #Rescaling of fraction of rock
    if np.sum(fm)!=100:
        fm=fm*100/(np.sum(fm))
    fm=fm*0.01  
    #Conditional on the module
    #Bounds
    if sat==0: #Dry case
        Mveff=Voigt(M,fm,phi,sat,module)
        Mreff=Reuss(M,fm,phi,sat,module)
    elif sat==1: #satuated case
        Mveff=Voigt(M,fm,phi,sat,module)
        Mreff=Reuss(M,fm,phi,sat,module)
    else:
        print("type 0: for dry\n type 1: for saturated")
        return
    #Den=Rock_den(den,fm,phi,sat)
    return Mveff,Mreff

def get_values():
    ax = plt.gca() # get axis handle
    line=ax.lines[0] # get the first line, there might be more
    x=line.get_xdata()
    y=line.get_ydata()
    plt.close()
    return x,y

def Voigt_Reuss_Hill(K,Mu,den,fm,sat):
    "Voigt-Reuss-Hill bounds given K, Mu, density, fractions, and saturation\
    saturation is either full (1) or none (0)."
    Vvoigt=[]
    Vreuss=[]
    Khill=[]
    Vhill=[]
    Dens=[]
    KR=[]
    KV=[]
    MR=[]
    MV=[]
    #Range of porosities
    phi_vals=np.arange(0,105,5)
    for i in range(0,len(phi_vals)):
        phi=phi_vals[i]
        if sat==0:
            dens,frac=Rock_den(den,fm,phi,0)
        elif sat==1:
            dens,frac=Rock_den(den,fm,phi,100)
        else:
            print('Set "sat" either to 0 or 1, please')
        Kv,Kr=Wood(K,frac,phi,sat,'Bulk')
        Mv,Mr=Wood(Mu,frac,phi,sat,'Shear')
        Khill=np.append(Khill,(Kv+Kr)/2)
        Vhill=np.append(Vhill,V_p((Kv+Kr)/2,0,dens))
        Vv=V_p(Kv,0,dens) #Mv
        Vr=V_p(Kr,0,dens) #Mr
        Vvoigt=np.append(Vvoigt,Vv)
        Vreuss=np.append(Vreuss,Vr)
        Dens=np.append(Dens,dens)
        KR=np.append(KR,Kr)
        KV=np.append(KV,Kv)
        MR=np.append(MR,Mr)
        MV=np.append(MV,Mv)
    #Output
    #plt.figure()
    plt.plot(phi_vals,Vreuss,'r')
    plt.plot(phi_vals,Vvoigt,'b')
    plt.plot(phi_vals,Vhill,'g')
    return
    
def Voigt_Reuss_Hill_mc(K2,Mu2,den2,fm2,sat):
    "Voigt-Reuss-Hill bounds given K, Mu, density, fractions, and saturation\
    saturation is either full (1) or none (0)."
    phi_vals=np.arange(2,82,5)
    Xr=[]
    Xr_int=[]
    Xr_mean=[]
    Yr=[]
    Xv=[]
    Xv_int=[]
    Xv_mean=[]
    Yv=[]
    Hill_mean=[]
    Hill_int=[]
    VHill=[]
    VHill_int=[]
    for i in range(0,len(phi_vals)):
        phi=mc.N(phi_vals[i],0.05*phi_vals[i])
        if sat==1:
            Dend2,fm=Rock_den(den2,fm2,phi,100)
        elif sat==0:
            Dend2,fm=Rock_den(den2,fm2,phi,0)
        else:
            print('Set "sat" either to 0 or 1, please')
        Kvd2,Krd2=Wood(K2,fm,phi,sat,'Bulk')  
        Mvd2,Mrd2=Wood(Mu2,fm,phi,sat,'Shear')
        #Kvs2,Krs2=Wood(K2,fm2,phi,1,'Bulk')  
        #Mvs2,Mrs2=Wood(Mu2,fm2,phi,1,'Shear')
        #Extract the values
        Vpdr2=V_p(Krd2,0,Dend2) #Mrd2
        Vpdr2.plot()
        x2,y2=get_values()
        Xr_mean=np.append(Xr_mean,Vpdr2.mean)
        Xr=np.append(Xr,x2)
        Xr_int=np.append(Xr_int,Vpdr2.percentile((0.025,0.975)))
        Yr=np.append(Yr,y2)
        Vpdv2=V_p(Kvd2,0,Dend2)#Mvd2
        Vpdv2.plot()
        x3,y3=get_values()
        Xv_mean=np.append(Xv_mean,Vpdv2.mean)
        Xv=np.append(Xv,x3)
        Xv_int=np.append(Xv_int,Vpdv2.percentile((0.025,0.975)))
        Yv=np.append(Yv,y3)
        #Hill Average
        Hill=(Kvd2+Krd2)/2
        Hill_mean=np.append(Hill_mean,Hill.mean)
        Hill_int=np.append(Hill_int,Hill.percentile((0.025,0.975)))
        VHill=np.append(VHill,V_p(Hill,0,Dend2).mean)
        VHill_int=np.append(VHill_int,V_p(Hill,0,Dend2).percentile((0.025,0.975)))
    #Reorder the arrays to extract the values
    XVv=np.reshape(Xv,[len(phi_vals),len(Xv)/len(phi_vals)])
    YVv=np.reshape(Yv,[len(phi_vals),len(Yv)/len(phi_vals)]) 
    XVv_int=np.reshape(Xv_int,[len(phi_vals),len(Xv_int)/len(phi_vals)])
    XVr_int=np.reshape(Xr_int,[len(phi_vals),len(Xr_int)/len(phi_vals)])
    HV_int=np.reshape(VHill_int,[len(phi_vals),len(Xr_int)/len(phi_vals)])
    #Output
    plt.figure('Voigt_Reuss limits')
    plt.plot(phi_vals,Xv_mean,'b')
    plt.plot(phi_vals,Xr_mean,'r')
    plt.plot(phi_vals,VHill,'g')
    plt.errorbar(phi_vals,Xv_mean, yerr=[np.abs(np.array(XVv_int)[:,0]-Xv_mean), np.abs(np.array(XVv_int)[:,1]-Xv_mean)], fmt='o',c='b')
    plt.errorbar(phi_vals,Xr_mean, yerr=[np.abs(np.array(XVr_int)[:,0]-Xr_mean), np.abs(np.array(XVr_int)[:,1]-Xr_mean)], fmt='o',c='r')
    plt.errorbar(phi_vals,VHill, yerr=[np.abs(np.array(HV_int)[:,0]-VHill), np.abs(np.array(HV_int)[:,1]-VHill)], fmt='o',c='g')
    plt.fill_between(phi_vals,np.array(XVv_int)[:,0],np.array(XVv_int)[:,1],alpha=0.3,color='b')
    plt.fill_between(phi_vals,np.array(XVr_int)[:,0],np.array(XVr_int)[:,1],alpha=0.3,color='r')
    plt.fill_between(phi_vals,np.array(HV_int)[:,0],np.array(HV_int)[:,1],alpha=0.3,color='g')
    plt.ylabel('Vp (m/s)'),plt.xlabel('Porosity')
    return

def Tek_read_waveform(channel,trigger,folder,name,plot):
    "Subroutine to read waveforms from the Tektronix DPO 3014 Oscilloscope \
    path: directory to save the files \
    name: sample name \
    plot: 1 = plot, 0 = no plot  \
    Example: \
    Tek_read_waveform(CH1,CH4,path,Al_big_cyl,plot) \
    Tek_read_waveform(CH1,CH4, C:\\Users\\edur409\\Desktop\\Waveforms, AL_bid_cyl.csv,1)"
    
    #Connect to Oscilloscope
    rm = visa.ResourceManager()
    rm.list_resources()
    inst = rm.open_resource('TCPIP::130.216.58.218::4000::SOCKET')
    #TCPIP[board]::host address::port::SOCKET
    print(inst.query('*IDN?'))
    
    #Read the trigger
    inst.write('DATA:ENCDG RIBINARY') #The encoding is RIBINARY (Big endian 16-bit signed)
    inst.write('DATA:SOU '+trigger+'')# Channel to read
    inst.write('DATA:WIDTH 2') #To read 2-bite data (16-bit numbers)
    inst.write('CURV?')
    #Information to convert the numbers to volts and time
    data = inst.read_raw()
    xincr = float(inst.ask('WFMPRE:XINCR?'))
    xzero=float(inst.ask('WFMPRE:XZERO?'))+xincr #Position of the trigger
    ymult = float(inst.ask('WFMPRE:YMULT?'))        
    yzero = float(inst.ask('WFMPRE:YZERO?'))
    yoff = float(inst.ask('WFMPRE:YOFF?'))        
    waveform =np.array(inst.query_binary_values('CURV?', datatype='h', is_big_endian=True))
    Volts = (waveform - yoff) * ymult  + yzero
    time = 1e6*(np.linspace(0, xincr * len(Volts), len(Volts))+xzero)
    Mt=np.array([time, Volts]).T #Array to save        
        
    #Set the variables to communicate
    inst.write('DATA:ENCDG RIBINARY') #The encoding is RIBINARY (Big endian 16-bit signed)
    inst.write('DATA:SOU '+channel+'')# Channel to read
    inst.write('DATA:WIDTH 2') #To read 2-bite data (16-bit numbers)
    inst.write('CURV?')
    
    #Information to convert the numbers to volts and time
    data = inst.read_raw()
    ymult = float(inst.ask('WFMPRE:YMULT?'))
    yzero = float(inst.ask('WFMPRE:YZERO?'))
    yoff = float(inst.ask('WFMPRE:YOFF?'))
    xincr = float(inst.ask('WFMPRE:XINCR?'))
    
    #Read the waveform ('h' works with a WIDTH of 2)
    waveform =np.array(inst.query_binary_values('CURV?', datatype='h', is_big_endian=True))
    
    #Transform to Volts
    Volts = (waveform - yoff) * ymult  + yzero
    time = 1e6*(np.linspace(0, xincr * len(Volts), len(Volts))+xzero)
    #PLot to check
    if plot==1:
        plt.plot(time,Volts) #What's on the screen
        plt.xlabel('time (microsecs)'),plt.ylabel('Volts')
    else:
        print('Set plot to 1 if you want to check the reading was OK')
    
    #Save the waveform
    M=np.array([time, Volts]).T #Array to save
    Rock=name #Name of your sample
    path=''+folder+'\\'+Rock+'.csv'
    #patht=''+folder+'\\'+Rock+'t.csv'
    np.savetxt(path, M, fmt='%1.8f',delimiter=',',header='time (microsecs), Volts')   
    #np.savetxt(patht, Mt, fmt='%1.8f',delimiter=',',header='time (microsecs), Volts')   
    return time, Volts

def plot_csv(name,plot):
    "plot_csv('tek0000ALL.csv')"
    A=np.loadtxt(name,delimiter=',',skiprows=21) #skiprows=1
    time=A[:,0]*1e6 #microsecs
    Signal=A[:,1]
    Trigger=A[:,2]
    #plot
    if plot==1:
        plt.figure(name)
        plt.plot(time,Signal,label='Signal')
        plt.plot(time,Trigger,label='Trigger')
        plt.xlabel('time (microsecs)')
        plt.legend()
    elif plot==0:
        print('no plot')
    else:
        print('Set 1 to plot the waveforms')
    return time,Trigger,Signal
    
def Length_sample(Lengths):
    Mean=np.around(np.mean(Lengths),decimals=3)
    Std=np.around(np.std(Lengths),decimals=3)
    print('('+np.str(Mean)+' +/- '+np.str(Std)+')')
    return Mean, Std

def Velocity_P(L,P_time,StdL,dtP):
    "Calulate the P wave velocity given length of sample and time of travel \
    \V,dV=Velocity_P(L,P_time,StdL,dtP), time in microseconds, length in cm"
    deci=0 #number of decimal places
    V=L*10000/P_time
    dV=V*np.sqrt((StdL/L)**2+(dtP/P_time)**2)
    print('P-Velocity: '+np.str(np.round(V,decimals=deci))+'+/-'+np.str(np.round(dV,decimals=deci))+' m/s')
    #Monte Carlo estimations
    Len=mc.N(L,StdL)    
    tP=mc.N(P_time,dtP)    
    Vmc=Len*1e4/tP
    V_down,V_up=Vmc.percentile((0.025,0.975)) #95% Confidence Interval
    print('P-Velocity (MC, 95% Credible Interval) : '+np.str(np.round(V_down,decimals=0))+' < '+np.str(np.round(Vmc.mean,decimals=0))+' < '+np.str(np.round(V_up,decimals=0))+' m/s')
    return V,dV,V_down,Vmc,V_up

def Dry_density(am,cm,da,dc):
    "Density and its error from Archimede's principle\
    a: mass of dry sample in air; c: mass of dry sample in water"
    rhow=1.#01
    drhow=0#.01
    den=am*rhow/(am-cm)
    den=np.round(den,decimals=2)
    delta_den=den*np.sqrt(((da/am)-(da/(am-cm)))**2+(dc/(am-cm))**2+(drhow/rhow)**2)
    delta_den=np.round(delta_den,decimals=3)
    print('Dry density (Archimedes): '+np.str(den)+'+/-'+np.str(delta_den)+'')    
    return den,delta_den

def Dry_density2(am,bm,dm,da,db,dd):
    "Density and its error from Archimede's principle\
    a: mass of dry sample in air; c: mass of dry sample in water"
    rhow=1.#01
    drhow=0#.01
    den=am*rhow/(bm-dm)
    den=np.round(den,decimals=2)
    delta_den=den*np.sqrt((da/am)**2+(db/(bm-dm))**2+(dd/(bm-dm))**2+(drhow/rhow)**2)
    delta_den=np.round(delta_den,decimals=3)
    print('Dry density (Archimedes): '+np.str(den)+'+/-'+np.str(delta_den)+'')    
    return den,delta_den
    
def Wet_density(bm,am,cm,db,da,dc):
    "Density and its error from Archimede's principle\
    b: mass of wet sample in air; a: mass of dry sample in air; \
    c: mass of dry sample in water"
    rhow=1.#01
    drhow=0#.01
    den=bm*rhow/(am-cm)
    den=np.round(den,decimals=2)
    delta_den=den*np.sqrt((db/bm)**2+(da/(am-cm))**2+(dc/(am-cm))**2+(drhow/rhow)**2)
    delta_den=np.round(delta_den,decimals=3)
    print('Wet density (Archimedes): '+np.str(den)+'+/-'+np.str(delta_den)+'')    
    return den,delta_den

def Wet_density2(bm,dm,db,dd):
    "Density and its error from Archimede's principle\
    b: mass of wet sample in air; a: mass of dry sample in air; \
    c: mass of dry sample in water"
    rhow=1.#01
    drhow=0#.01
    den=bm*rhow/(bm-dm)
    den=np.round(den,decimals=2)
    delta_den=den*np.sqrt(((db/bm)-(db/(bm-dm)))**2+(dd/(bm-dm))**2+(drhow/rhow)**2)
    delta_den=np.round(delta_den,decimals=3)
    print('Wet density (Archimedes): '+np.str(den)+'+/-'+np.str(delta_den)+'')    
    return den,delta_den
    
def Min_density(bm,am,cm,db,da,dc):
    "Density and its error from Archimede's principle\
    b: mass of wet sample in air; a: mass of dry sample in air; \
    c: mass of dry sample in water"
    rhow=1.#01
    drhow=0#.01
    den=am*rhow/((am-cm)-(bm-am))
    den=np.round(den,decimals=2)
    delta_den=den*np.sqrt(((da/am)-(2*da/(2*am-cm-bm)))**2+(db/(2*am-cm-bm))**2+(dc/(2*am-cm-bm))**2+(drhow/rhow)**2)
    delta_den=np.round(delta_den,decimals=3)
    print('Mineral density (Archimedes): '+np.str(den)+'+/-'+np.str(delta_den)+'')    
    return den,delta_den

def Min_density2(am,dm,da,dd):
    "Density and its error from Archimede's principle\
    b: mass of wet sample in air; a: mass of dry sample in air; \
    c: mass of dry sample in water"
    rhow=1.#01
    drhow=0#.01
    den=am*rhow/(am-dm)
    den=np.round(den,decimals=2)
    delta_den=den*np.sqrt(((da/am)-(da/(am-dm)))**2+(dd/(am-dm))**2+(drhow/rhow)**2)
    delta_den=delta_den#np.round(delta_den,decimals=3)
    print('Mineral density (Archimedes): '+np.str(den)+'+/-'+np.str(delta_den)+'')    
    return den,delta_den

def Porosity(bm,am,cm,db,da,dc):
    por=100*(bm-am)/(am-cm)
    por=np.round(por,decimals=2)
    delta_por=por*np.sqrt(((da/(bm-am))+(da/(am-cm)))**2+(db/(bm-am))**2+(dc/(am-cm))**2)
    delta_por=np.round(delta_por,decimals=2)
    print('Porosity (Archimedes): '+np.str(por)+'+/-'+np.str(delta_por)+'')
    return por, delta_por

def Porosity2(bm,am,dm,db,da,dd):
    por=100*(bm-am)/(bm-dm)
    por=np.round(por,decimals=2)
    delta_por=por*np.sqrt(((db/(bm-am))+(db/(bm-dm)))**2+(da/(bm-am))**2+(dd/(bm-dm))**2)
    delta_por=np.round(delta_por,decimals=2)
    print('Porosity (Archimedes): '+np.str(por)+'+/-'+np.str(delta_por)+'')
    return por, delta_por

def histogram(v,title): 
    x=np.linspace(np.min(v)-2*np.std(v),np.max(v)+2*np.std(v),50)
    y=normal.pdf(x,np.mean(v),np.std(v)) #range, mean, std
    plt.figure(title)
    plt.hist(v,normed=True,color='b')
    plt.plot(x,y,'r',lw=2)
    return np.mean(v),np.std(v)
    
def S_values(coords):
    'Input the time coordinates from the DTW picking'
    tS=[]
    tD=[]
    for i in range(0,len(coords)):
        tS=np.append(tS,coords[i][0])
        tD=np.append(tD,coords[i][1])
    #histogram(tD,'S-arrival (dry)')
    #histogram(tS,'S-arrival (sat.)')
    return tD,tS
    
def P_values(coords):
    'Input the time coordinates from the DTW picking'
    time=[]
    for i in range(0,len(coords)):
        time=np.append(time,coords[i][0])
    histogram(time,'P-arrival')
    return time
    
def Velocity_S(L,S_time,StdL,dtS):
    "Calulate the S wave velocity given length of sample and time of travel \
    \V,dV=Velocity_S(L,S_time,StdL,dtP), time in microseconds, length in cm"
    deci=0 #number of decimal places
    V=L*10000/S_time
    dV=V*np.sqrt((StdL/L)**2+(dtS/S_time)**2)
    print('S-Velocity: '+np.str(np.round(V,decimals=deci))+'+/-'+np.str(np.round(dV,decimals=deci))+' m/s')
    #Monte Carlo estimations
    Len=mc.N(L,StdL)    
    tS=mc.N(S_time,dtS)    
    Vmc=Len*1e4/tS
    V_down,V_up=Vmc.percentile((0.025,0.975)) #95% Confidence Interval
    print('S-Velocity (MC, 95% Credible Interval) : '+np.str(np.round(V_down,decimals=0))+' < '+np.str(np.round(Vmc.mean,decimals=0))+' < '+np.str(np.round(V_up,decimals=0))+' m/s')
    return V,dV,V_down,Vmc,V_up
    
def Vp_Vs(Vp,dVp,Vs,dVs):
    r=np.round(Vp/Vs,decimals=2)
    dr=np.round(r*np.sqrt((dVp/Vp)**2+(dVs/Vs)**2),decimals=2)
    return r,dr

def load_csv(filename,plot):
    "load a csv file removing the first line header"
    A=np.loadtxt(filename,delimiter=',',skiprows=1)
    time=A[:,0]
    Volts=A[:,1]
    if plot==1:
        plt.plot(time,Volts) #What's on the screen
        plt.xlabel('time (microsecs)'),plt.ylabel('Volts')
    else:
        print('No plot')
    return time,Volts

def extract_parameters(filename):
    "Extract Pressure, Cycle, Saturation state, and phase measured from file name\
    name,Pc,phase,state,cycle=extract_parameters(filename)\
    filename = string with the name of the file to read"
    c=0
    p=[]
    for i in range(0,len(filename)):
        if filename[i]=='_':
            c=c+1
            p=np.append(p,i)
    name=filename[0:np.int(p[0])]
    Pc=np.float(filename[np.int(p[0])+1:np.int(p[1])])
    phase=filename[np.int(p[1])+1:np.int(p[2])]
    state=filename[np.int(p[2])+1:np.int(p[3])] #state=filename[np.int(p[2])+1:np.int(p[2])+4]
    end=filename.find('.') #End of the string
    cycle=filename[np.int(p[3])+1:end]        
    return name,Pc,phase,state,cycle
    
def plot_waveforms(filenames,Pc,phases,cycles,phase,Pconf,Lengths):
    indices=np.where((Pc==Pconf)&(phases==phase))[0]
    "plot_waveforms('RR5_02','P',3000,Lengths)"
    plt.figure('Waveforms '+np.str(Pconf)+' psi')
    for i in range(0,len(indices)):    
        cycle=cycles[indices[i]]
        name=filenames[indices[i]]
        #print(name)
        timeP, P=load_csv(name,0)
        plt.subplot(3,1,1),plt.plot(timeP,P,label=cycle)
        plt.xlabel('time (microsecs)')
        plt.legend()
        tP=np.loadtxt('P_time_Picks_'+name+'.out')
        x=np.linspace(np.min(tP)-2*np.std(tP),np.max(tP)+2*np.std(tP),50)
        y=normal.pdf(x,np.mean(tP),np.std(tP)) #range, mean, std
        #subplot(2,1,2),plt.hist(tP,normed=True,color='b')
        plt.subplot(3,1,2),plt.plot(x,y,lw=2,label=cycle)
        plt.legend()
        plt.xlabel('time (microsecs)')
        #
        L,StdL=Length_sample(Lengths)
        VelP,DVelP=Velocity_P(L,np.mean(tP),StdL,np.std(tP))
        x=np.linspace(VelP-10*DVelP,VelP+10*DVelP,100)
        y=normal.pdf(x,VelP,DVelP)
        plt.subplot(3,1,3),plt.plot(x,y,lw=2,label=cycle)
        plt.xlabel(''+phase+' Velocity (m/s)')
        plt.legend()
    return

def plot_onthego(Figure,pos,wc,c,s):
    "Plot the waveforms collected on the spot. Inputs are strings \
    wc: wave component (PP,S1,S2)\
    c: cycle measured (u1,d1,u2,d2)\
    s: state of saturation (dry, sat500, etc..)\
    pos: position of plot in subplot figure (e.g. 111)"
    
    # unlike matlab, in python you have to call the functions/modules you
    # are going to use:
    import numpy as np # the numerical capabilities of python
    import matplotlib.pyplot as plt # for plotting 
    import glob
    from . import RockPhysics as rp
    #from obspy.core import Trace, Stats, Stream
    from scipy import signal
    
    # Set default font size for plots:
    font = {'size'   : 24}
    plt.rc('font',**font)

    # put all file names to read in a list. All files for one rock and one fluid (all pressures) should be in ONE folder as csv files
    files = sorted(glob.glob('*.csv'),key=lambda f: list(filter(str.isdigit, f)))#sort by pressure number
    
    #Extract the relevant information from the filename
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
    
    indices=np.where((phases==wc)& (cycles==c)& (states==s))[0]    
    indices_sorted=[x for _,x in sorted(zip(Pc[indices],indices))] #Sort according to pressures   
    
    # define values:
    labels = []
    offset=0
    fig= plt.figure(Figure)
    ax = fig.add_subplot(pos)
    # read and plot data in each file:
    for i in range(0,len(indices_sorted)):
        file=files[indices_sorted[i]]
        # read in time/amp:
        t,V=rp.load_csv(file,0)#readcsvdata(file)
        
        # rather arbitrary filter with params to tune:
        #b, a = signal.butter(2, 0.005, 'low')
        #V = signal.filtfilt(b, a, V)

        # normalises the traces:
        V= 1*V/np.max(V)

        # plot in \mu s:
        plt.plot(t,V+offset,'k')

        # offset traces vertically:
        offset = offset +1

        # create a custom y-label, based on pressure in file name:
        labels.append(Pc[indices_sorted[i]])  # 32:35 for 142,

    # this is very manual to set the pressure values on the y-axis. It
    # depends crucially on the offset value and the files/pressures
    # recorded:
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Pressure (psi)')
    plt.title(phases[indices_sorted[i]]+' from rock '+name)
    plt.grid()
    plt.axis('tight')
    #plt.savefig(name+'_'+phases[indices[i]]+'.pdf')
    plt.xlim(np.min(t),np.max(t))
    plt.show()
    return    
    
#def Poisson(Vp,dVp,Vs,dVs): #Check this subroutine
#    Pr=np.round((Vp**2-2*Vs**2)/2*(Vp**2-Vs**2),decimals=2)
#    dPr1=Pr*((2*Vp*dVp/(Vp**2-2*Vs**2))-((4*Vp**2*dVp)/2*(Vp**2-Vs**2))) #Partial Vp
#    dPr2=Pr*((4*Vs*dVs/(2*(Vp**2-Vs**2)))-4*Vs*dVs/(Vp**2-2*Vs**2)) #Partial Vs
#    dPr=np.round(np.sqrt(dPr1**2+dPr2**2),decimals=2)
#    return Pr,dPr
