# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 19:51:59 2017

@author: edur409
"""

import numpy as np
import matplotlib.pyplot as plt
import mcerp3 as mc
from . import RockPhysics as rp

def Vol_porosity(sample):
    #Load the geometrical parmeters for volume 
    A=np.loadtxt('./Pycnometer_data/'+sample+'_10cycles/Vol_params.txt',skiprows=1)
    L=np.mean(A[:,0]) #Average length of sample
    dL=np.std(A[:,0]) #standard deviation of length of sample
    D=np.mean(A[:,1]) #Average length of sample
    dD=np.std(A[:,1]) #standard deviation of diameter of sample
    #Volume of sample and its error
    V=0.25*np.pi*L*D**2
    dV=V*np.sqrt(4*(dD/D)**2+(dL/L)**2)
    #Load the eskeletal volumes of the pycnometer
    B=np.loadtxt('./'+sample+'/Volume_Pycnometer.out')
    Vp=B[0] #Volume of pycnometer
    dVp=B[1] #error in volume of pycnometer
    #Calculate the porosity 
    phi=100*(V-Vp)/V
    dphi=phi*np.sqrt((dVp/(V-Vp))**2+(((1/(V-Vp))-1/V)*dV)**2)
    np.savetxt('./'+sample+'/Porosity_Pycnometer.out', [phi,dphi], fmt='%1.2f',delimiter=',',header='phi (%) , dphi (%)') 
    #Calculate Bulk Density from geometry and mass of sample
    M=np.loadtxt('./Pycnometer_data/'+sample+'_10cycles/Mass_'+sample+'.txt')    
    dM=0.002 
    rhob=M/V
    drhob=rhob*np.sqrt((dM/M)**2+(dV/V)**2)
    np.savetxt('./'+sample+'/Density_Pycnometer_bulk.out', [rhob,drhob], fmt='%1.2f',delimiter=',',header='rhob (g/cc) , drhob (g/cc)')     
    return V,Vp

#Read the name of the samples
F=np.loadtxt('Samples_Depth.txt',dtype={'names': ('sample', 'depth'),'formats': ('S20', 'f4')})

for i in range(0,len(F)):
    #Choose the sample
    sample=F[i][0].decode('UTF-8')#input("Type name of sample (e.g. 'NM11_2087_4A'): \n")
    
    #Density measurements
    a=np.loadtxt('./'+sample+'/'+sample+'_a.txt',delimiter=',')
    b=np.loadtxt('./'+sample+'/'+sample+'_b.txt',delimiter=',')
    c=np.loadtxt('./'+sample+'/'+sample+'_c.txt',delimiter=',')
    d=np.loadtxt('./'+sample+'/'+sample+'_d.txt',delimiter=',')
    
    a_mean=np.mean(a)
    b_mean=np.mean(b)
    c_mean=np.mean(c)
    d_mean=np.mean(d)
    da=1*np.std(a)
    db=1*np.std(b)
    dc=1*np.std(c)
    dd=1*np.std(d)
    
    #densities (Standard error)
    #Dry_den, delta_dry=rp.Dry_density(a_mean,c_mean,da,dc) #Dry density
    Dry_den, delta_dry=rp.Dry_density2(a_mean,b_mean,d_mean,da,db,dd)
    #Wet_den,delta_wet=rp.Wet_density(b_mean,a_mean,c_mean,db,da,dc) #Wet density
    Wet_den,delta_wet=rp.Wet_density2(b_mean,d_mean,db,dd)
    #Min_den,delta_min=rp.Min_density(b_mean,a_mean,c_mean,db,da,dc) #Mineral density
    Min_den,delta_min=rp.Min_density2(a_mean,d_mean,da,dd)
    #por,delta_por=rp.Porosity(b_mean,a_mean,c_mean,db,da,dc) #Porosity
    por,delta_por=rp.Porosity2(b_mean,a_mean,d_mean,db,da,dd)
    
    #Densities (MC)
    A=mc.N(a_mean,da)
    C=mc.N(c_mean,dc)
    B=mc.N(b_mean,db)
    D=mc.N(d_mean,dd)
    #Dry density (MC)
    RW=1
    Dry_denMC=A*RW/(B-D)
    DryDen_down,DryDen_up=Dry_denMC.percentile((0.025,0.975)) #95% Confidence Interval
    print(('Dry density (MC, 95% Credible Interval): '+np.str(np.round(DryDen_down,decimals=2))+' < '+np.str(np.round(Dry_denMC.mean,decimals=2))+' < '+np.str(np.round(DryDen_up,decimals=2))+''))    
    #Wet density (MC)    
    Wet_denMC=B*RW/(B-D)
    WetDen_down,WetDen_up=Wet_denMC.percentile((0.025,0.975)) #95% Confidence Interval
    print(('Wet density (MC, 95% Credible Interval): '+np.str(np.round(WetDen_down,decimals=2))+' < '+np.str(np.round(Wet_denMC.mean,decimals=2))+' < '+np.str(np.round(WetDen_up,decimals=2))+''))    
    #Mineral density (MC)    
    Min_denMC=A*RW/((B-D)-(B-A))
    MinDen_down,MinDen_up=Min_denMC.percentile((0.025,0.975)) #95% Confidence Interval
    print(('Min. density (MC, 95% Credible Interval): '+np.str(np.round(MinDen_down,decimals=2))+' < '+np.str(np.round(Min_denMC.mean,decimals=2))+' < '+np.str(np.round(MinDen_up,decimals=2))+''))    
    #Porosity (MC)
    Por_MC=100*(B-A)/(B-D)
    Por_down,Por_up=Por_MC.percentile((0.025,0.975)) #95% Confidence Interval
    print(('Porosity (MC, 95% Credible Interval): '+np.str(np.round(Por_down,decimals=2))+' < '+np.str(np.round(Por_MC.mean,decimals=2))+' < '+np.str(np.round(Por_up,decimals=2))+''))    
    
    #Save the densities in ASCII files
    #Save the velocities and intervals
    np.savetxt('./'+sample+'/Dry_density2.out', [Dry_den,delta_dry,DryDen_down,Dry_denMC.mean,DryDen_up], fmt='%1.2f',delimiter=',',header='Rho (dry), Std_Rho (dry), Rhod_0025, Rhod_mean, Rhod_0975')
    np.savetxt('./'+sample+'/Sat_density2.out', [Wet_den,delta_wet,WetDen_down,Wet_denMC.mean,WetDen_up], fmt='%1.2f',delimiter=',',header='Rho (sat.), Std_Rho (sat.), Rhos_0025, Rhos_mean, Rhos_0975')
    np.savetxt('./'+sample+'/Min_density2.out', [Min_den,delta_min,MinDen_down,Min_denMC.mean,MinDen_up], fmt='%1.2f',delimiter=',',header='Rho (min.), Std_Rho (min.), Rhom_0025, Rhom_mean, Rhom_0975')
    np.savetxt('./'+sample+'/Porosity2.out', [por,delta_por,Por_down,Por_MC.mean,Por_up], fmt='%1.2f',delimiter=',',header='Por (%), Std_Por (%), Por_0025, PorMC_mean, Por_0975')
    
    #Plots
    #Dry Density
    plt.figure('Dry Density')
    Dry_den=mc.N(Dry_den,delta_dry)
    Dry_den.plot(label='Density (dry)',lw=2,color='b')
    Dry_denMC.plot(hist=True,label='Density (dry) (MC)',color='g')
    plt.legend()
    plt.xlabel('Density (g/cm3)')
    plt.savefig('./'+sample+'/'+sample+'_Dry_Den2.pdf',bbox_inches='tight')
    #plt.show(block=True)
    
    #Saturated Density
    plt.figure('Saturated Density')
    Wet_den=mc.N(Wet_den,delta_wet)
    Wet_den.plot(label='Density (sat.)',lw=2,color='b')
    Wet_denMC.plot(hist=True,label='Density (sat.) (MC)',color='g')
    plt.legend()
    plt.xlabel('Density (g/cm3)')
    plt.savefig('./'+sample+'/'+sample+'_Saturated_Den2.pdf',bbox_inches='tight')
    #plt.show(block=True)
    
    #Mineral Density
    plt.figure('Mineral Density')
    Min_den=mc.N(Min_den,delta_min)
    Min_den.plot(label='Density (min.)',lw=2,color='b')
    Min_denMC.plot(hist=True,label='Density (min.) (MC)',color='g')
    plt.legend()
    plt.xlabel('Density (g/cm3)')
    plt.savefig('./'+sample+'/'+sample+'_Mineral_Den2.pdf',bbox_inches='tight')
    #plt.show(block=True)
    
    #Porosity
    plt.figure('Porosity')
    por=mc.N(por,delta_por)
    por.plot(label='Porosity (%)',lw=2,color='b')
    Por_MC.plot(hist=True,label='Porosity (%) (MC)',color='g')
    plt.legend()
    plt.xlabel('Porosity (%)')
    plt.savefig('./'+sample+'/'+sample+'_Porosity2.pdf',bbox_inches='tight')
    #plt.show(block=True)
    plt.close('all')
    
    #Porosity Pycnometer    
    #sample='NM11_2087_4C'
    V,Vp=Vol_porosity(sample)
    

    
    
    
