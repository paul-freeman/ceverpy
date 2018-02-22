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
Vsd=[]
dVsd=[]
phi=[]
dphi=[]
Vss=[]
dVss=[]
#Densities
rhod=[]
drhod=[]
rhos=[]
drhos=[]
rhom=[]
drhom=[]

for i in range(0,len(A)):
    sample=A[i][0].decode('UTF-8')
    #Read porosities    
    phi=np.append(phi,np.loadtxt('./'+sample+'/Porosity.out')[0])
    dphi=np.append(dphi,np.loadtxt('./'+sample+'/Porosity.out')[1])   
   #Read dry Vp velocities
    Vpd=np.append(Vpd,np.loadtxt('./'+sample+'/'+sample+'d''/Velocities_P.out')[0])
    dVpd=np.append(dVpd,np.loadtxt('./'+sample+'/'+sample+'d''/Velocities_P.out')[1])
    #Read sat Vp velocities
    Vps=np.append(Vps,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_P.out')[0])
    dVps=np.append(dVps,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_P.out')[1])
    #Read the S Velocities dry
    Vsd=np.append(Vsd,np.loadtxt('./'+sample+'/'+sample+'d''/Velocities_S_dry_wav.out')[0])
    dVsd=np.append(dVsd,np.loadtxt('./'+sample+'/'+sample+'d''/Velocities_S_dry_wav.out')[1])    
    #read S-velocity saturated
    Vss=np.append(Vss,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_S_sat_wav.out')[0])
    dVss=np.append(dVss,np.loadtxt('./'+sample+'/'+sample+'s''/Velocities_S_sat_wav.out')[1])    
    #Density#
    rhod=np.append(rhod,np.loadtxt('./'+sample+'/Dry_density.out')[0])
    drhod=np.append(drhod,np.loadtxt('./'+sample+'/Dry_density.out')[1])   
    #Density saturated
    rhos=np.append(rhos,np.loadtxt('./'+sample+'/Sat_density.out')[0])
    drhos=np.append(drhos,np.loadtxt('./'+sample+'/Sat_density.out')[1])
    #Density minerals
    rhom=np.append(rhom,np.loadtxt('./'+sample+'/Min_density.out')[0])
    drhom=np.append(drhom,np.loadtxt('./'+sample+'/Min_density.out')[1])
    
#plot S-waveforms (dry and saturated)
for i in range(0,len(A)):
    sample=A[i][0].decode('UTF-8')
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
1
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
