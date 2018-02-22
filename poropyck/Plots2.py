# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 20:43:16 2017

@author: edur409
"""
def fourier_trans(dt,S):
    f=np.fft.fftfreq(len(S),d=dt)
    As=np.abs(np.fft.fft(S))
    Ap=np.angle(np.fft.fft(S))#np.unwrap(np.angle(np.fft.fft(S)))
    return f*1e-6,As,Ap


def fourier_phase(dt,Sd,Ss):
    fd,Asd,Psd=fourier_trans(dt,Sd)
    fs,Ass,Pss=fourier_trans(dt,Ss)
    #Give the functions the same spectra
    PS=200*np.ones((len(Sd),1))
    SSPS=PS*np.exp(1j*Psd)
    SDPS=PS*np.exp(1j*Pss)
    #Transforming back
    Sss=np.real(np.fft.ifft(np.fft.ifftshift(SSPS)))
    Sds=np.real(np.fft.ifft(np.fft.ifftshift(SDPS)))
    return Sds,Sss

    
#Import the Rock Physics subroutines
from . import RockPhysics as rp #Loads all the necessary files
import numpy as np
import matplotlib.pyplot as plt
import mcerp3 as mc

#Font size for plots
font= {'size' : 24}
plt.rc('font',**font)

A=np.loadtxt('Samples_Depth.txt',dtype={'names': ('sample', 'depth'),'formats': ('S20', 'f4')})
plt.close('all')
#value = input("Choose folder: Press 1 (dry rock), Press 2 (saturated rock) \n")
#if value==1:
#    print("Picking for dry rock (" + np.str(value) + ")")
#    folder=''+well+'d'
#elif value==2:
#    print("Picking for saturated rock (" + np.str(value) + ")")
#    folder=''+well+'s'
#else:
#    print("Please, follow the instructions. Try again!")
#    exit()

Vpd=[]
dVpd=[]
Vps=[]
dVps=[]
Vpsc=[] #Vp sat x-corr
dVpsc=[]
Vsd=[]
dVsd=[]
Vss=[]
dVss=[]
Vssc=[] #Vs sat x-corr
dVssc=[]
#Densities
rhod=[]
drhod=[]
rhos=[]
drhos=[]
rhom=[]
drhom=[]
rhomp=[]
drhomp=[]
#Porosities
phi=[]
dphi=[]
phip=[]
dphip=[]


for i in range(0,len(A)):
    sample=A[i][0].decode('UTF-8')
    #Read porosities    
    phi=np.append(phi,np.loadtxt('./'+sample+'/Porosity2.out')[0])
    dphi=np.append(dphi,np.loadtxt('./'+sample+'/Porosity2.out')[1])  
    #Read porosities    
    phip=np.append(phip,np.loadtxt('./'+sample+'/Porosity_Pycnometer.out')[0])
    dphip=np.append(dphip,np.loadtxt('./'+sample+'/Porosity_Pycnometer.out')[1])  
    #Read dry Vp velocities
    Vpd=np.append(Vpd,np.loadtxt('./'+sample+'/'+sample+'d''/Velocities_P_dry_both.out')[0])
    dVpd=np.append(dVpd,np.loadtxt('./'+sample+'/'+sample+'d''/Velocities_P_dry_both.out')[1])
    #Read sat Vp velocities
    Vps=np.append(Vps,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_P_sat_both.out')[0])
    dVps=np.append(dVps,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_P_sat_both.out')[1])
    #Read sat Vp velocities x-corr
    Vpsc=np.append(Vpsc,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_P_sat_Xcorr.out')[0])
    dVpsc=np.append(dVpsc,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_P_sat_Xcorr.out')[1])
    #Read the S Velocities dry
    Vsd=np.append(Vsd,np.loadtxt('./'+sample+'/'+sample+'d''/Velocities_S_dry_both.out')[0])
    dVsd=np.append(dVsd,np.loadtxt('./'+sample+'/'+sample+'d''/Velocities_S_dry_both.out')[1])    
    #read S-velocity saturated
    Vss=np.append(Vss,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_S_sat_both.out')[0])
    dVss=np.append(dVss,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_S_sat_both.out')[1])    
    #read S-velocity saturated x-corr
    Vssc=np.append(Vssc,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_S_sat_Xcorr.out')[0])
    dVssc=np.append(dVssc,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_S_sat_Xcorr.out')[1])        
    #Density#
    rhod=np.append(rhod,np.loadtxt('./'+sample+'/Dry_density.out')[0])
    drhod=np.append(drhod,np.loadtxt('./'+sample+'/Dry_density.out')[1])   
    #Density saturated
    rhos=np.append(rhos,np.loadtxt('./'+sample+'/Sat_density.out')[0])
    drhos=np.append(drhos,np.loadtxt('./'+sample+'/Sat_density.out')[1])
    #Density minerals
    rhom=np.append(rhom,np.loadtxt('./'+sample+'/Min_density2.out')[0])
    drhom=np.append(drhom,np.loadtxt('./'+sample+'/Min_density2.out')[1])
    #Density of pycnometer
    rhomp=np.append(rhomp,np.loadtxt('./'+sample+'/Density_Pycnometer.out')[0])
    drhomp=np.append(drhomp,np.loadtxt('./'+sample+'/Density_Pycnometer.out')[1])    
    
#plot S-waveforms (dry and saturated)
for i in range(0,len(A)):
    #Waveforms
    timeSd,T,Sd=rp.plot_csv('./'+sample+'/'+sample+'d/tek0002ALL.csv',0)
    timeSs,T,Ss=rp.plot_csv('./'+sample+'/'+sample+'s/tek0002ALL.csv',0)
    plt.figure('Waveforms')
    plt.plot(timeSd,Sd,label='Sd',c='b')
    plt.plot(timeSs,Ss,label='Ss',c='g')
    plt.axis('tight')
    plt.title('S-waveforms sample '+sample+'')
    plt.legend()
    plt.xlabel('time ($\mu$s)')
    plt.savefig('./'+sample+'/S-waveforms.pdf',bbox_inches='tight')
    plt.close()
    #Spectra of waveforms
    dt=1e-8
    fd,Asd,Psd=fourier_trans(dt,Sd)
    fs,Ass,Pss=fourier_trans(dt,Ss)
    plt.figure('Amplitude Spectra')
    plt.plot(fd,Asd,label='Amp.Spectrum Sd')
    plt.plot(fs,Ass,label='Amp. Spectrum Ss')
    plt.legend()
    plt.title('S-waveforms sample '+sample+'')
    plt.axis([0,1.0,0,np.max([Asd,Ass])])    
    plt.savefig('./'+sample+'/AmpSpectra-S-waveforms.pdf',bbox_inches='tight')
    plt.close()  
    #Phase Spectra
    plt.figure('Phase Spectra')
    plt.plot(fd,Psd,label='Ph.Spectrum Sd')
    plt.plot(fs,Pss,label='Ph. Spectrum Ss')
    plt.axis([0,1.0,np.min([Psd,Pss]),np.max([Psd,Pss])])
    plt.legend()
    plt.savefig('./'+sample+'/PhaseSpectra-S-waveforms.pdf',bbox_inches='tight')
    plt.close()
#    #Same Amplitude Spectrum
#    Sds,Sss=fourier_phase(dt,Sd,Ss)
#    plt.figure('Same Amplitude Spectrum')
#    plt.plot(timeSd,Sds,label='Sd same amp.')
#    plt.plot(timeSs,Sss,label='Ss same amp.')
#    #plt.axis([0,1.0,np.min([Psd,Pss]),np.max([Psd,Pss])])
#    plt.legend()
#    plt.savefig('./'+sample+'/Waveforms-S-same_amp.pdf',bbox_inches='tight')
#    plt.close()

font=9  
size=40    
#Plot Vp vs. phi (labels with wells)
fig = plt.figure('Vp (dry) vs. phi')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(phi[i],Vpd[i],s=size,zorder=1)
    ax.text(phi[i],Vpd[i], A[i][0].decode('UTF-8'), fontsize=font) 
ax.errorbar(phi, Vpd, yerr=[2*dVpd, 2*dVpd], xerr=[2*dphi, 2*dphi], fmt='o')
plt.scatter(phi,Vpd,s=4*size,c=rhom,zorder=2)
plt.xlabel('Porosity (%)'),plt.ylabel('P Velocity dry (m/s)')
plt.axis('tight')
plt.colorbar()  
plt.grid('on')  
plt.axis([0,np.max(phi)+np.max(dphi)+1,3500,5500])    
    
#Plot Vp vs. phi (labels with wells)
fig = plt.figure('Vp (sat) vs. phi')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(phi[i],Vps[i],s=size)
    ax.text(phi[i],Vps[i], A[i][0].decode('UTF-8'), fontsize=font) 
ax.errorbar(phi, Vps, yerr=[2*dVps, 2*dVps], xerr=[2*dphi, 2*dphi], fmt='o')
plt.scatter(phi,Vps,s=4*size,c=rhom,zorder=2)
plt.xlabel('Porosity (%)'),plt.ylabel('P Velocity sat. (m/s)')        
plt.axis('tight')
plt.colorbar()
plt.grid('on')  
plt.axis([0,np.max(phi)+np.max(dphi)+1,3500,5500])
    
fig = plt.figure('Vs (dry) vs. phi')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(phi[i],Vsd[i],s=size)
    ax.text(phi[i],Vsd[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(phi, Vsd, yerr=[2*dVsd, 2*dVsd], xerr=[2*dphi, 2*dphi], fmt='o')
plt.scatter(phi,Vsd,s=4*size,c=rhom,zorder=2)
plt.xlabel('Porosity (%)'),plt.ylabel('S Velocity dry (m/s)')    
plt.axis('tight')
plt.colorbar()
plt.grid('on')  
plt.axis([0,np.max(phi)+np.max(dphi)+1,2000,3500])

fig = plt.figure('Vs (sat) vs. phi')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(phi[i],Vss[i],s=size)
    ax.text(phi[i],Vss[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(phi, Vss, yerr=[2*dVss, 2*dVss], xerr=[2*dphi, 2*dphi], fmt='o')
plt.scatter(phi,Vss,s=4*size,c=rhom,zorder=2)
plt.xlabel('Porosity (%)'),plt.ylabel('S Velocity sat (m/s)')        
plt.axis('tight')
plt.colorbar()
plt.grid('on')  
plt.axis([0,np.max(phi)+np.max(dphi)+1,2000,3500])

fig = plt.figure('Vp (sat) vs. rhos')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(rhos[i],Vps[i],s=size)
    ax.text(rhos[i],Vps[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(rhos, Vps, yerr=[2*dVps, 2*dVps], xerr=[2*drhos, 2*drhos], fmt='o')
plt.xlabel('Sat. Density (g/cc)'),plt.ylabel('P Velocity sat (m/s)')      
plt.scatter(rhos,Vps,s=4*size,c=rhom,zorder=2)
plt.axis('tight')
plt.colorbar()
plt.grid('on')  
plt.axis([np.min(rhod)-np.min(drhod)-0.2,np.max(rhos)+np.max(drhos)+0.2,3500,5500])

fig = plt.figure('Vs (sat) vs. rhos')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(rhos[i],Vss[i],s=size)
    ax.text(rhos[i],Vss[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(rhos, Vss, yerr=[2*dVss, 2*dVss], xerr=[2*drhos, 2*drhos], fmt='o')
plt.xlabel('Sat. Density (g/cc)'),plt.ylabel('S Velocity sat (m/s)') 
plt.scatter(rhos,Vss,s=4*size,c=rhom,zorder=2)
plt.axis('tight')
plt.colorbar()
plt.grid('on')  
plt.axis([np.min(rhod)-np.min(drhod)-0.2,np.max(rhos)+np.max(drhos)+0.2,2000,3500])

fig = plt.figure('Vp (dry) vs. rhod')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(rhod[i],Vpd[i],s=size)
    ax.text(rhod[i],Vpd[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(rhod, Vpd, yerr=[2*dVpd, 2*dVpd], xerr=[2*drhod, 2*drhod], fmt='o')
plt.xlabel('Dry Density (g/cc)'),plt.ylabel('P Velocity dry (m/s)')      
plt.scatter(rhod,Vpd,s=4*size,c=rhom,zorder=2)
plt.axis('tight')
plt.colorbar()
plt.grid('on')  
plt.axis([np.min(rhod)-np.min(drhod)-0.2,np.max(rhos)+np.max(drhos)+0.2,3500,5500])

fig = plt.figure('Vs (dry) vs. rhod')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(rhod[i],Vsd[i],s=size)
    ax.text(rhod[i],Vsd[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(rhod, Vsd, yerr=[2*dVsd, 2*dVsd], xerr=[2*drhod, 2*drhod], fmt='o')
plt.xlabel('Dry Density (g/cc)'),plt.ylabel('S Velocity dry (m/s)') 
plt.scatter(rhod,Vsd,s=4*size,c=rhom,zorder=2)
plt.axis('tight')
plt.colorbar()
plt.grid('on')  
plt.axis([np.min(rhod)-np.min(drhod)-0.2,np.max(rhos)+np.max(drhos)+0.2,2000,3500])

Rd,dRd=rp.Vp_Vs(Vpd,dVpd,Vsd,dVsd)    
fig = plt.figure('Vp/Vs vs. phi')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(phi[i],Rd[i],s=size)
    ax.text(phi[i],Rd[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(phi, Rd, yerr=[2*dRd, 2*dRd], xerr=[2*dphi, 2*dphi], fmt='o')
plt.xlabel('Porosity (%)'),plt.ylabel('Vp/Vs dry')    
plt.grid('on')  
plt.axis('tight')

Rs,dRs=rp.Vp_Vs(Vps,dVps,Vss,dVss)    
fig = plt.figure('Vp/Vs vs. phi (sat)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(phi[i],Rs[i],s=size)
    ax.text(phi[i],Rs[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(phi, Rs, yerr=[2*dRs, 2*dRs], xerr=[2*dphi, 2*dphi], fmt='o')
plt.xlabel('Porosity (%)'),plt.ylabel('Vp/Vs sat')    
plt.grid('on')  
plt.axis('tight')


fig = plt.figure('Phi vs. rhom')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(rhom[i],phi[i],s=size)
    ax.text(rhom[i],phi[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(rhom, phi, yerr=[2*dphi, 2*dphi], xerr=[2*drhom, 2*drhom], fmt='o')
plt.scatter(rhom,phi,s=4*size,c=Rs,zorder=2)
plt.xlabel('Mineral density (g/cc)'),plt.ylabel('Porosity (%)')    
plt.axis('tight')

#PRd,dPRd=rp.Poisson(Vpd,dVpd,Vsd,dVsd) 
#fig = plt.figure('PR vs. phi')
#ax = fig.add_subplot(111)   
#for i in range(0,len(A)-1):
#    ax.scatter(phi[i],PRd[i],s=size)
#    ax.text(phi[i],PRd[i], A[i][0].decode('UTF-8'), fontsize=font)
#ax.errorbar(phi, PRd, yerr=[2*dPRd, 2*dPRd], xerr=[2*dphi, 2*dphi], fmt='o')
#plt.xlabel('Porosity (%)'),plt.ylabel('Vp/Vs dry')    
#plt.axis('tight')

#Plot Vp vs. phi (labels with wells)
fig = plt.figure('Vp (dry) vs. Vs (dry)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(Vpd[i],Vsd[i],s=size)
    ax.text(Vpd[i],Vsd[i], A[i][0].decode('UTF-8'), fontsize=font) 
ax.errorbar(Vpd, Vsd, yerr=[2*dVsd, 2*dVsd], xerr=[2*dVpd, 2*dVpd], fmt='o')
plt.scatter(Vpd,Vsd,s=4*size,c=phi,zorder=2)
plt.xlabel('Vp dry (m/s)'),plt.ylabel('Vs dry. (m/s)')        
plt.axis('tight')
plt.colorbar()
plt.grid('on')  
plt.axis([3500,5500,2000,3500])

fig = plt.figure('Vp (sat) vs. Vs (sat)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(Vps[i],Vss[i],s=size)
    ax.text(Vps[i],Vss[i], A[i][0].decode('UTF-8'), fontsize=font) 
ax.errorbar(Vps, Vss, yerr=[2*dVss, 2*dVss], xerr=[2*dVps, 2*dVps], fmt='o')
plt.scatter(Vps,Vss,s=4*size,c=phi,zorder=2)
plt.xlabel('Vp sat. (m/s)'),plt.ylabel('Vs sat. (m/s)')        
plt.axis('tight')
plt.colorbar()
plt.grid('on')  
plt.axis([3500,5500,2000,3500])

fig = plt.figure('DTW vs. X-corr')
ax = fig.add_subplot(111) 
for i in range(0,len(A)): 
    ax.scatter(Vss,Vssc,s=size)
    ax.scatter(Vps,Vpsc,s=size,c='r')
ax.errorbar(Vss, Vssc, yerr=[2*dVssc, 2*dVssc], xerr=[2*dVss, 2*dVss], fmt='o',c='b')
ax.errorbar(Vps, Vpsc, yerr=[2*dVpsc, 2*dVpsc], xerr=[2*dVps, 2*dVps], fmt='o',c='r')
x=np.linspace(1984,5250,10)
ax.plot(x,x,'g'),plt.axis('equal'),plt.grid('on'),plt.axis('tight')
plt.xlabel('Sat. Velocity DTW (m/s)'),plt.ylabel('Sat. Velocity x-corr (m/s)')
plt.savefig('DTW_Vs_XCORR.pdf',bbox_inches='tight')


fig = plt.figure('DTW vs. X-corr (Poro)')
ax = fig.add_subplot(111) 
for i in range(0,len(A)): 
    ax.scatter(Vss,Vssc,s=size)
    ax.scatter(Vps,Vpsc,s=size)
ax.errorbar(Vss, Vssc, yerr=[2*dVssc, 2*dVssc], xerr=[2*dVss, 2*dVss], fmt='o',c='b')
plt.scatter(Vss,Vssc,s=4*size,c=phi,zorder=2)
ax.errorbar(Vps, Vpsc, yerr=[2*dVpsc, 2*dVpsc], xerr=[2*dVps, 2*dVps], fmt='o',c='r')
plt.scatter(Vps,Vpsc,s=4*size,c=phi,zorder=2)
x=np.linspace(1984,5250,10)
ax.plot(x,x,'g'),plt.axis('equal'),plt.grid('on'),plt.axis('tight')
plt.xlabel('Sat. Velocity DTW (m/s)'),plt.ylabel('Sat. Velocity x-corr (m/s)')
plt.colorbar()
plt.savefig('DTW_Vs_XCORR.pdf',bbox_inches='tight')

#THE PLOTS FOR THE DTW AND XCORR COMPARISON
fig = plt.figure('Vs (sat.) vs. Vs (dry) DTW')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(Vss[i],Vsd[i],s=size)
    #ax.text(Vss[i],Vsd[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(Vss, Vsd, yerr=[2*dVsd, 2*dVsd], xerr=[2*dVss, 2*dVss], fmt='o',c='b',zorder=0)
plt.scatter(Vss,Vsd,s=4*size,c=phi,zorder=1,cmap='gray')
x=np.linspace(1800,3200,10)
plt.plot(x,x,'g')
plt.axis([1800,3200,1800,3200])
plt.grid('on')#plt.axis('equal'),plt.grid('on')#,plt.axis('tight'),
plt.xlabel('Vs (sat.) DTW (m/s)'),plt.ylabel('Vs (dry) (m/s)')
plt.colorbar()
plt.savefig('Vssat_Vs_Vsdry_DTW.pdf',bbox_inches='tight')

fig = plt.figure('Vs (sat.) vs. Vs (dry) XCORR')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(Vssc[i],Vsd[i],s=size)
    #ax.text(Vssc[i],Vsd[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(Vssc, Vsd, yerr=[2*dVsd, 2*dVsd], xerr=[2*dVss, 2*dVss], fmt='o',c='b',zorder=0)
plt.scatter(Vssc,Vsd,s=4*size,c=phi,zorder=1,cmap='gray')
x=np.linspace(1800,3200,10)
plt.plot(x,x,'g')
plt.axis([1800,3200,1800,3200])
plt.grid('on')#plt.axis('equal'),plt.axis('tight'),plt.grid('on')
plt.xlabel('Vs (sat.) XCORR (m/s)'),plt.ylabel('Vs (dry) (m/s)')
plt.colorbar()
plt.savefig('Vssat_Vs_Vsdry_XCORR.pdf',bbox_inches='tight')

fig = plt.figure('Vp (sat.) vs. Vp (dry) DTW')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(Vps[i],Vpd[i],s=size)
    #ax.text(Vps[i],Vpd[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(Vps, Vpd, yerr=[2*dVpd, 2*dVpd], xerr=[2*dVps, 2*dVps], fmt='o',c='b',zorder=0)
plt.scatter(Vps,Vpd,s=4*size,c=phi,zorder=1,cmap='gray')
x=np.linspace(3200,5300,10)
plt.plot(x,x,'g')
plt.axis([3200,5300,3200,5300])
plt.grid('on')#plt.axis('equal'),plt.axis('tight'),plt.grid('on')
plt.xlabel('Vp (sat.) DTW (m/s)'),plt.ylabel('Vp (dry) (m/s)')
plt.colorbar()
plt.savefig('Vpsat_Vs_Vpdry_DTW.pdf',bbox_inches='tight')

fig = plt.figure('Vp (sat.) vs. Vp (dry) XCORR')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(Vpsc[i],Vpd[i],s=size)
    #ax.text(Vpsc[i],Vpd[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(Vpsc, Vpd, yerr=[2*dVpd, 2*dVpd], xerr=[2*dVps, 2*dVps], fmt='o',c='b',zorder=0)
plt.scatter(Vpsc,Vpd,s=4*size,c=phi,zorder=1,cmap='gray')
x=np.linspace(3200,5300,10)
plt.plot(x,x,'g')
plt.axis([3200,5300,3200,5300])
plt.grid('on')#plt.axis('equal'),plt.axis('tight'),plt.grid('on')
plt.xlabel('Vp (sat.) XCORR (m/s)'),plt.ylabel('Vp (dry) (m/s)')
plt.colorbar()
plt.savefig('Vpsat_Vs_Vpdry_XCORR.pdf',bbox_inches='tight')


#Rho pycnometry vs. Rho Archimedes
fig = plt.figure('rhomp vs. rhom')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(rhom[i],rhomp[i],s=size)
    #ax.text(rhom[i],rhomp[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(rhom, rhomp, yerr=[2*drhomp, 2*drhomp], xerr=[2*drhom, 2*drhom], fmt='o',zorder=0)
plt.scatter(rhom,rhomp,s=4*size,c=phi,zorder=1,cmap='gray')
plt.colorbar()
x=np.linspace(np.min([rhom,rhomp]),np.max([rhom,rhomp]),20)
y=np.linspace(np.min([rhom,rhomp]),np.max([rhom,rhomp]),20)
plt.plot(x,y)
plt.xlabel('Mineral density Archimedes ($g/cm^3$)'),plt.ylabel('Mineral density Pycnometer ($g/cm^3$)')    
plt.axis('tight'),plt.grid('on')
plt.savefig('RhomPycVsRhom2.pdf',bbox_inches='tight')

#Porosity pycnometry vs. Porosity Archimedes
fig = plt.figure('phip vs. phi')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(phi[i],phip[i],s=size)
    #ax.text(phi[i],phip[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(phi, phip, yerr=[2*dphip, 2*dphip], xerr=[2*dphi, 2*dphi], fmt='o',zorder=0)
plt.scatter(phi,phip,s=4*size,c=Vps,zorder=1,cmap='gray')
plt.colorbar()
x=np.linspace(np.min([phi,phip]),np.max([phi,phip]),20)
y=np.linspace(np.min([phi,phip]),np.max([phi,phip]),20)
plt.plot(x,y)
plt.xlabel('Porosity Archimedes ($\%$)'),plt.ylabel('Porosity Pycnometer ($\%$)')    
plt.axis('tight'),plt.grid('on')
plt.savefig('PoroPycVsPoro2.pdf',bbox_inches='tight')


fig = plt.figure('Vs (sat.) DTW vs. Vs (sat) XCORR')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(Vssc[i],Vss[i],s=size)
    ax.text(Vssc[i],Vss[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(Vssc, Vss, yerr=[2*dVss, 2*dVss], xerr=[2*dVssc, 2*dVssc], fmt='o',c='b')
plt.scatter(Vssc,Vss,s=4*size,c=phi,zorder=2)
x=np.linspace(2000,3600,10)
plt.plot(x,x,'g')
plt.axis('equal'),plt.axis('tight'),plt.grid('on')
plt.xlabel('Vs (sat.) XCORR m/s'),plt.ylabel('Vs (sat.) DTW m/s')
plt.colorbar()

fig = plt.figure('Vp (sat.) vs. Vp (sat) XCORR & DTW')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(Vpsc[i],Vpd[i],s=size)
    ax.scatter(Vps[i],Vpd[i],s=size,c='r')
    ax.text(Vpsc[i],Vpd[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(Vpsc, Vpd, yerr=[2*dVpd, 2*dVpd], xerr=[2*dVpsc, 2*dVpsc], fmt='o',c='b')
plt.errorbar(Vps, Vpd, yerr=[2*dVpd, 2*dVpd], xerr=[2*dVps, 2*dVps], fmt='o',c='r')
plt.scatter(Vpsc,Vpd,s=4*size,c=phi,zorder=2)
x=np.linspace(3000,5600,10)
plt.plot(x,x,'g')
plt.axis('equal'),plt.axis('tight'),plt.grid('on')
plt.xlabel('Vp (sat.) XCORR m/s'),plt.ylabel('Vp (sat.) DTW m/s')
plt.colorbar()

#LOOK AT THIS PLOT AFTER MODIFYING
fig = plt.figure('Vs (sat.) vs. Vs (dry) XCORR & DTW')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(Vssc[i],Vsd[i],s=size)
    ax.scatter(Vss[i],Vsd[i],s=size,c='r')
    ax.text(Vssc[i],Vsd[i], A[i][0].decode('UTF-8')+np.str(phi[i]), fontsize=font) 
    ax.text(Vss[i],Vsd[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(Vssc, Vsd, yerr=[2*dVsd, 2*dVsd], xerr=[2*dVssc, 2*dVssc], fmt='o',c='b')
plt.errorbar(Vss, Vsd, yerr=[2*dVsd, 2*dVsd], xerr=[2*dVss, 2*dVss], fmt='o',c='r')
plt.scatter(Vssc,Vsd,s=4*size,c=phi,zorder=2)
x=np.linspace(2000,3200,10)
plt.plot(x,x,'g')
plt.axis('equal'),plt.axis('tight'),plt.grid('on')
plt.xlabel('Vs (sat.) m/s'),plt.ylabel('Vs (dry) m/s')
plt.colorbar()

fig = plt.figure('Vp (sat.) vs. Vp (dry) XCORR & DTW')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(Vpsc[i],Vpd[i],s=size)
    ax.scatter(Vps[i],Vpd[i],s=size,c='r')
    ax.text(Vpsc[i],Vpd[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(Vpsc, Vpd, yerr=[2*dVpd, 2*dVpd], xerr=[2*dVpsc, 2*dVpsc], fmt='o',c='b')
plt.errorbar(Vps, Vpd, yerr=[2*dVpd, 2*dVpd], xerr=[2*dVps, 2*dVps], fmt='o',c='r')
plt.scatter(Vpsc,Vpd,s=4*size,c=phi,zorder=2)
x=np.linspace(3000,5600,10)
plt.plot(x,x,'g')
plt.axis('equal'),plt.axis('tight'),plt.grid('on')
plt.xlabel('Vp (sat.) m/s'),plt.ylabel('Vp (dry) m/s')
plt.colorbar()


fig = plt.figure('Vp (sat) vs. Vs (sat) DTW & XCORR')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(Vps[i],Vss[i],s=size,c='b')
    ax.scatter(Vpd[i],Vsd[i],s=size,c='g')
    ax.scatter(Vpsc[i],Vssc[i],s=size,c='r')
    ax.text(Vps[i],Vss[i], A[i][0].decode('UTF-8'), fontsize=font) 
    ax.text(Vpsc[i],Vssc[i], A[i][0].decode('UTF-8'), fontsize=font) 
ax.errorbar(Vps, Vss, yerr=[2*dVss, 2*dVss], xerr=[2*dVps, 2*dVps], fmt='o',color='b')
ax.errorbar(Vpsc, Vssc, yerr=[2*dVssc, 2*dVssc], xerr=[2*dVpsc, 2*dVpsc], fmt='o',color='r')
#plt.scatter(Vps,Vss,s=4*size,c=phi,zorder=2)
plt.xlabel('Vp sat. (m/s)'),plt.ylabel('Vs sat. (m/s)')        
plt.axis('tight')
#plt.colorbar()
plt.grid('on')  
#plt.axis([3000,5500,2000,3500])

##Print statement
#print('%4.2f +/- %4.2f ' % (Vpd[0],dVpd[0]))
#s='%4.2f +/- %4.2f ' % (Vpd[0],dVpd[0]) #Print the same into a string
#
#filename='Table.txt'
#F=open(filename,'w')
#F.write('Well, V_p (dry), V_p (sat) XCORR, V_p (sat.) DTW, V_s (dry), V_s (sat.) XCORR, V_s (sat.) DTW')
#S=[]
#for i in range(0,len(A)):
#    name=A[i][0].decode('UTF-8')
#    s1='%4.0f +/- %3.0f ' % (Vpd[i],dVpd[i]) #Print the same into a string
#    s2='%4.0f +/- %3.0f ' % (Vpsc[i],dVpsc[i]) #Print the same into a string
#    s3='%4.0f +/- %3.0f ' % (Vps[i],dVps[i]) #Print the same into a string
#    s4='%4.0f +/- %3.0f ' % (Vsd[i],dVsd[i]) #Print the same into a string
#    s5='%4.0f +/- %3.0f ' % (Vssc[i],dVssc[i]) #Print the same into a string
#    s6='%4.0f +/- %3.0f ' % (Vss[i],dVss[i]) #Print the same into a string
#    F.write('%s, %s, %s, %s, %s, %s, %s '% (name, s1 ,  s2,  s3,  s4,  s5,  s6))
#F.close()

#Table of elastic constants



fig = plt.figure('Vs (sat.) vs. Vs (dry) XCORR & DTW no names')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(Vssc[i],Vsd[i],s=size)
    ax.scatter(Vss[i],Vsd[i],s=size,c='r')
    #ax.text(Vssc[i],Vsd[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(Vssc, Vsd, yerr=[2*dVsd, 2*dVsd], xerr=[2*dVssc, 2*dVssc], fmt='o',c='b')
plt.errorbar(Vss, Vsd, yerr=[2*dVsd, 2*dVsd], xerr=[2*dVss, 2*dVss], fmt='o',c='r')
#plt.scatter(Vssc,Vsd,s=4*size,c=phi,zorder=2)
x=np.linspace(2000,3600,10)
plt.plot(x,x,'g')
plt.axis('equal'),plt.axis('tight'),plt.grid('on')
plt.xlabel('Vs (sat.) m/s'),plt.ylabel('Vs (dry) m/s')
plt.colorbar()
