# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 09:12:59 2017

@author: leo
"""

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
rhob=[]
drhob=[]
#Porosities
phi=[]
dphi=[]
phip=[]
dphip=[]
#Elastic constants
Kd=[]
dKd=[]
Ks=[]
dKs=[]
Mud=[]
dMud=[]
Mus=[]
dMus=[]


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
    rhod=np.append(rhod,np.loadtxt('./'+sample+'/Dry_density2.out')[0])
    drhod=np.append(drhod,np.loadtxt('./'+sample+'/Dry_density2.out')[1])   
    #Density saturated
    rhos=np.append(rhos,np.loadtxt('./'+sample+'/Sat_density2.out')[0])
    drhos=np.append(drhos,np.loadtxt('./'+sample+'/Sat_density2.out')[1])
    #Density minerals
    rhom=np.append(rhom,np.loadtxt('./'+sample+'/Min_density2.out')[0])
    drhom=np.append(drhom,np.loadtxt('./'+sample+'/Min_density2.out')[1])
    #Density of pycnometer
    rhomp=np.append(rhomp,np.loadtxt('./'+sample+'/Density_Pycnometer.out')[0])
    drhomp=np.append(drhomp,np.loadtxt('./'+sample+'/Density_Pycnometer.out')[1]) 
    #Bulk Density of pycnometer
    rhob=np.append(rhob,np.loadtxt('./'+sample+'/Density_Pycnometer_bulk.out')[0])
    drhob=np.append(drhob,np.loadtxt('./'+sample+'/Density_Pycnometer_bulk.out')[1]) 
    #Bulk moduli
    Kd=np.append(Kd,np.loadtxt('./'+sample+'/'+sample+'d''/Bulk_dry_Pyc.out')[0])
    dKd=np.append(dKd,np.loadtxt('./'+sample+'/'+sample+'d''/Bulk_dry_Pyc.out')[1])
    Ks=np.append(Ks,np.loadtxt('./'+sample+'/'+sample+'s''/Bulk_sat_Pyc.out')[0])
    dKs=np.append(dKs,np.loadtxt('./'+sample+'/'+sample+'s''/Bulk_sat_Pyc.out')[1])
    #Shear moduli
    Mud=np.append(Mud,np.loadtxt('./'+sample+'/'+sample+'d''/Shear_dry_Pyc.out')[0])
    dMud=np.append(dMud,np.loadtxt('./'+sample+'/'+sample+'d''/Shear_dry_Pyc.out')[1])
    Mus=np.append(Mus,np.loadtxt('./'+sample+'/'+sample+'s''/Shear_sat_Pyc.out')[0])
    dMus=np.append(dMus,np.loadtxt('./'+sample+'/'+sample+'s''/Shear_sat_Pyc.out')[1])
    
#Plots
size=40
font=8
#
fig = plt.figure('Bulk moduli (porosity pyc.)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(phip[i],Ks[i],s=size)
    ax.scatter(phip[i],Kd[i],s=size,c='r')
    ax.text(phip[i],Kd[i], A[i][0].decode('UTF-8'), fontsize=font) 
    ax.text(phip[i],Ks[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(phip, Ks, yerr=[2*dKs, 2*dKs], xerr=[2*dphip, 2*dphip], fmt='o',c='b',label='Sat.')
plt.errorbar(phip, Kd, yerr=[2*dKd, 2*dKd], xerr=[2*dphip, 2*dphip], fmt='o',c='r',label='Dry')
plt.xlabel('Porosity Pyc. (%)'),plt.ylabel('Bulk modulus (GPa)')
plt.legend()

fig = plt.figure('Shear moduli (porosity pyc.)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(phip[i],Mus[i],s=size)
    ax.scatter(phip[i],Mud[i],s=size,c='r')
    ax.text(phip[i],Mud[i], A[i][0].decode('UTF-8'), fontsize=font)
    ax.text(phip[i],Mus[i], A[i][0].decode('UTF-8'), fontsize=font)
plt.errorbar(phip, Mus, yerr=[2*dKs, 2*dKs], xerr=[2*dphip, 2*dphip], fmt='o',c='b',label='Sat.')
plt.errorbar(phip, Mud, yerr=[2*dKd, 2*dKd], xerr=[2*dphip, 2*dphip], fmt='o',c='r',label='Dry')
plt.xlabel('Porosity Pyc. (%)'),plt.ylabel('Shear modulus (GPa)')
plt.legend()

plt.set_cmap('gray')
fig = plt.figure('Bulk moduli (rhom pyc.)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(rhom[i],Ks[i],s=size/2)
    #ax.scatter(rhom[i],Kd[i],s=size/2,c='r')
    #ax.text(rhom[i],Kd[i], A[i][0].decode('UTF-8'), fontsize=font) 
    #ax.text(rhom[i],Ks[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(rhom, Ks, yerr=[2*dKs, 2*dKs], xerr=[2*drhom, 2*drhom], fmt='o',c='b',label='Sat.',zorder=0,lw=2)
#plt.errorbar(rhom, Kd, yerr=[2*dKd, 2*dKd], xerr=[2*drhom, 2*drhom], fmt='o',c='r',label='Dry',zorder=0,lw=2)
#plt.scatter(rhom,Kd,s=3*size,c=phip,zorder=1)
plt.scatter(rhom,Ks,s=4*size,c=phip,zorder=1)
plt.xlabel('Mineral Density Pyc. ($g/cm^3$)'),plt.ylabel('Bulk modulus (GPa)')
plt.legend(loc=2)
plt.colorbar()

fig = plt.figure('Shear moduli (rhom pyc.)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(rhom[i],Mus[i],s=size/2)
    #ax.scatter(rhom[i],Mud[i],s=size/2,c='r')
    #ax.text(rhom[i],Mud[i], A[i][0].decode('UTF-8'), fontsize=font) 
    #ax.text(rhom[i],Mus[i], A[i][0].decode('UTF-8'), fontsize=font)
plt.errorbar(rhom, Mus, yerr=[2*dMus, 2*dMus], xerr=[2*drhom, 2*drhom], fmt='o',c='b',label='Sat.',zorder=0,lw=2)
#plt.errorbar(rhom, Mud, yerr=[2*dMud, 2*dMud], xerr=[2*drhom, 2*drhom], fmt='o',c='r',label='Dry',zorder=0,lw=2)
#plt.scatter(rhom,Mud,s=3*size,c=phip,zorder=1)
plt.scatter(rhom,Mus,s=4*size,c=phip,zorder=1)
plt.xlabel('Mineral Density Pyc. ($g/cm^3$)'),plt.ylabel('Shear modulus (GPa)')
plt.legend(loc=2)
plt.colorbar()


fig = plt.figure('Bulk moduli (densities)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(rhos[i],Ks[i],s=size)
    ax.scatter(rhod[i],Kd[i],s=size,c='r')
    ax.text(rhod[i],Kd[i], A[i][0].decode('UTF-8'), fontsize=font) 
    ax.text(rhos[i],Ks[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(rhos, Ks, yerr=[2*dKs, 2*dKs], xerr=[2*drhos, 2*drhos], fmt='o',c='b',label='Sat.')
plt.errorbar(rhod, Kd, yerr=[2*dKd, 2*dKd], xerr=[2*drhod, 2*drhod], fmt='o',c='r',label='Dry')
plt.xlabel('Density Arch. ($g/cm^3$)'),plt.ylabel('Bulk modulus (GPa)')
plt.legend(loc=2)

fig = plt.figure('Shear moduli (densities)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(rhos[i],Mus[i],s=size)
    ax.scatter(rhod[i],Mud[i],s=size,c='r')
    ax.text(rhod[i],Mud[i], A[i][0].decode('UTF-8'), fontsize=font) 
    ax.text(rhos[i],Mus[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(rhos, Mus, yerr=[2*dMus, 2*dMus], xerr=[2*drhos, 2*drhos], fmt='o',c='b',label='Sat.')
plt.errorbar(rhod, Mud, yerr=[2*dMud, 2*dMud], xerr=[2*drhod, 2*drhod], fmt='o',c='r',label='Dry')
plt.xlabel('Density Arch. ($g/cm^3$)'),plt.ylabel('Shear modulus (GPa)')
plt.legend(loc=2)

fig = plt.figure('Bulk moduli (bulk densities dry)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    #ax.scatter(rhob[i],Ks[i],s=size)
    ax.scatter(rhob[i],Kd[i],s=size,c='r')
    ax.text(rhob[i],Kd[i], A[i][0].decode('UTF-8'), fontsize=font) 
#plt.errorbar(rhos, Ks, yerr=[2*dKs, 2*dKs], xerr=[2*drhos, 2*drhos], fmt='o',c='b',label='Sat.')
plt.errorbar(rhob, Kd, yerr=[2*dKd, 2*dKd], xerr=[2*drhob, 2*drhob], fmt='o',c='r',label='Dry')
plt.xlabel('Bulk Density ($g/cm^3$)'),plt.ylabel('Bulk modulus (GPa)')
plt.legend(loc=2)

fig = plt.figure('Shear moduli (bulk densities dry)')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    #ax.scatter(rhob[i],Ks[i],s=size)
    ax.scatter(rhob[i],Mud[i],s=size,c='r')
    ax.text(rhob[i],Mud[i], A[i][0].decode('UTF-8'), fontsize=font) 
#plt.errorbar(rhos, Ks, yerr=[2*dKs, 2*dKs], xerr=[2*drhos, 2*drhos], fmt='o',c='b',label='Sat.')
plt.errorbar(rhob, Mud, yerr=[2*dMud, 2*dMud], xerr=[2*drhob, 2*drhob], fmt='o',c='r',label='Dry')
plt.xlabel('Bulk Density ($g/cm^3$)'),plt.ylabel('Shear modulus (GPa)')
plt.legend(loc=2)


Rd,dRd=rp.Vp_Vs(Vpd,dVpd,Vsd,dVsd)    
fig = plt.figure('Vp/Vs vs. phi')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(phi[i],Rd[i],s=size)
    ax.text(phi[i],Rd[i], A[i][0].decode('UTF-8'), fontsize=8)
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

fig = plt.figure('Vp (sat) vs. rhom')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(rhom[i],Vps[i],s=size)
    ax.text(rhom[i],Vps[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(rhom, Vps, yerr=[2*dVps, 2*dVps], xerr=[2*drhom, 2*drhom], fmt='o')
plt.xlabel('Sat.m Density (g/cc)'),plt.ylabel('P Velocity sat (m/s)')      
plt.scatter(rhom,Vps,s=4*size,c=rhom,zorder=2)
plt.axis('tight')
plt.colorbar()
plt.grid('on')  

fig = plt.figure('Vs (sat) vs. rhom')
ax = fig.add_subplot(111)   
for i in range(0,len(A)-1):
    ax.scatter(rhom[i],Vss[i],s=size)
    ax.text(rhom[i],Vss[i], A[i][0].decode('UTF-8'), fontsize=font)
ax.errorbar(rhom, Vss, yerr=[2*dVss, 2*dVss], xerr=[2*drhom, 2*drhom], fmt='o')
plt.xlabel('Sat.m Density (g/cc)'),plt.ylabel('S Velocity sat (m/s)')      
plt.scatter(rhom,Vss,s=4*size,c=rhom,zorder=2)
plt.axis('tight')
plt.colorbar()
plt.grid('on')  

#plt.axis([np.min(rhod)-np.min(drhod)-0.2,np.max(rhos)+np.max(drhos)+0.2,3500,5500])

plt.set_cmap('gray')
fig = plt.figure('Shear moduli (rhom pyc.) all')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(rhom[i],Mus[i],s=size/2)
    ax.scatter(rhom[i],Mud[i],s=size/2,c='r')
    #ax.text(rhom[i],Mud[i], A[i][0].decode('UTF-8'), fontsize=font) 
    #ax.text(rhom[i],Mus[i], A[i][0].decode('UTF-8'), fontsize=font)
plt.errorbar(rhom, Mus, yerr=[2*dMus, 2*dMus], xerr=[2*drhom, 2*drhom], fmt='o',c='b',label='Sat.',zorder=0,lw=2)
plt.errorbar(rhom, Mud, yerr=[2*dMud, 2*dMud], xerr=[2*drhom, 2*drhom], fmt='o',c='r',label='Dry',zorder=0,lw=2)
plt.scatter(rhom,Mud,s=3*size,c=phip,zorder=1)
plt.scatter(rhom,Mus,s=4*size,c=phip,zorder=1)
plt.xlabel('Mineral Density Pyc. ($g/cm^3$)'),plt.ylabel('Shear modulus (GPa)')
plt.legend(loc=2)
plt.colorbar()


fig = plt.figure('Bulk moduli (rhom pyc.) all')
ax = fig.add_subplot(111)   
for i in range(0,len(A)):
    ax.scatter(rhom[i],Ks[i],s=size/2)
    ax.scatter(rhom[i],Kd[i],s=size/2,c='r')
    #ax.text(rhom[i],Kd[i], A[i][0].decode('UTF-8'), fontsize=font) 
    #ax.text(rhom[i],Ks[i], A[i][0].decode('UTF-8'), fontsize=font) 
plt.errorbar(rhom, Ks, yerr=[2*dKs, 2*dKs], xerr=[2*drhom, 2*drhom], fmt='o',c='b',label='Sat.',zorder=0,lw=2)
plt.errorbar(rhom, Kd, yerr=[2*dKd, 2*dKd], xerr=[2*drhom, 2*drhom], fmt='o',c='r',label='Dry',zorder=0,lw=2)
plt.scatter(rhom,Kd,s=3*size,c=phip,zorder=1)
plt.scatter(rhom,Ks,s=4*size,c=phip,zorder=1)
plt.xlabel('Mineral Density Pyc. ($g/cm^3$)'),plt.ylabel('Bulk modulus (GPa)')
plt.legend(loc=2)
plt.colorbar()
