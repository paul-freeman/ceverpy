# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 19:51:59 2017

@author: edur409
"""

import numpy as np
import matplotlib.pyplot as plt
import mcerp3 as mc
from . import RockPhysics as rp

def densities_porosity():
    #Choose the well
    well=input("Type name of sample (e.g. 'NM11_2087_4A'): \n")

    #Density measurements
    a=np.loadtxt('./'+well+'/'+well+'_a.txt',delimiter=',')
    b=np.loadtxt('./'+well+'/'+well+'_b.txt',delimiter=',')
    c=np.loadtxt('./'+well+'/'+well+'_c.txt',delimiter=',')
    d=np.loadtxt('./'+well+'/'+well+'_d.txt',delimiter=',')

    a_mean=np.mean(a)
    b_mean=np.mean(b)
    c_mean=np.mean(c)
    d_mean=np.mean(d)
    da=1*np.std(a)
    db=1*np.std(b)
    dc=1*np.std(c)
    dd=1*np.std(d)

    #densities (Standard error)
    Dry_den, delta_dry=rp.Dry_density(a_mean,c_mean,da,dc) #Dry density
    Wet_den,delta_wet=rp.Wet_density(b_mean,a_mean,c_mean,db,da,dc) #Wet density
    Min_den,delta_min=rp.Min_density(b_mean,a_mean,c_mean,db,da,dc) #Mineral density
    por,delta_por=rp.Porosity(b_mean,a_mean,c_mean,db,da,dc) #Porosity

    #Densities (MC)
    A=mc.N(a_mean,da)
    C=mc.N(c_mean,dc)
    B=mc.N(b_mean,db)
    D=mc.N(d_mean,dd)
    #Dry density (MC)
    RW=1
    Dry_denMC=A*RW/(A-C)
    DryDen_down,DryDen_up=Dry_denMC.percentile((0.025,0.975)) #95% Confidence Interval
    print(('Dry density (MC, 95% Credible Interval): '+np.str(np.round(DryDen_down,decimals=2))+' < '+np.str(np.round(Dry_denMC.mean,decimals=2))+' < '+np.str(np.round(DryDen_up,decimals=2))+''))    
    #Wet density (MC)    
    Wet_denMC=B*RW/(A-C)
    WetDen_down,WetDen_up=Wet_denMC.percentile((0.025,0.975)) #95% Confidence Interval
    print(('Wet density (MC, 95% Credible Interval): '+np.str(np.round(WetDen_down,decimals=2))+' < '+np.str(np.round(Wet_denMC.mean,decimals=2))+' < '+np.str(np.round(WetDen_up,decimals=2))+''))    
    #Mineral density (MC)    
    Min_denMC=A*RW/((A-C)-(B-A))
    MinDen_down,MinDen_up=Min_denMC.percentile((0.025,0.975)) #95% Confidence Interval
    print(('Min. density (MC, 95% Credible Interval): '+np.str(np.round(MinDen_down,decimals=2))+' < '+np.str(np.round(Min_denMC.mean,decimals=2))+' < '+np.str(np.round(MinDen_up,decimals=2))+''))    
    #Porosity (MC)
    Por_MC=100*(B-A)/(A-C)
    Por_down,Por_up=Por_MC.percentile((0.025,0.975)) #95% Confidence Interval
    print(('Porosity (MC, 95% Credible Interval): '+np.str(np.round(Por_down,decimals=2))+' < '+np.str(np.round(Por_MC.mean,decimals=2))+' < '+np.str(np.round(Por_up,decimals=2))+''))    

    #Save the densities in ASCII files
    #Save the velocities and intervals
    np.savetxt('./'+well+'/Dry_density.out', [Dry_den,delta_dry,DryDen_down,Dry_denMC.mean,DryDen_up], fmt='%1.2f',delimiter=',',header='Rho (dry), Std_Rho (dry), Rhod_0025, Rhod_mean, Rhod_0975')
    np.savetxt('./'+well+'/Sat_density.out', [Wet_den,delta_wet,WetDen_down,Wet_denMC.mean,WetDen_up], fmt='%1.2f',delimiter=',',header='Rho (sat.), Std_Rho (sat.), Rhos_0025, Rhos_mean, Rhos_0975')
    np.savetxt('./'+well+'/Min_density.out', [Min_den,delta_min,MinDen_down,Min_denMC.mean,MinDen_up], fmt='%1.2f',delimiter=',',header='Rho (min.), Std_Rho (min.), Rhom_0025, Rhom_mean, Rhom_0975')
    np.savetxt('./'+well+'/Porosity.out', [por,delta_por,Por_down,Por_MC.mean,Por_up], fmt='%1.2f',delimiter=',',header='Por (%), Std_Por (%), Por_0025, PorMC_mean, Por_0975')

    #Plots
    #Dry Density
    plt.figure('Dry Density')
    Dry_den=mc.N(Dry_den,delta_dry)
    Dry_den.plot(label='Density (dry)',lw=2,color='b')
    Dry_denMC.plot(hist=True,label='Density (dry) (MC)',color='g')
    plt.legend()
    plt.xlabel('Density (g/cm3)')
    plt.savefig('./'+well+'/'+well+'_Dry_Den.pdf',bbox_inches='tight')
    plt.show(block=True)

    #Saturated Density
    plt.figure('Saturated Density')
    Wet_den=mc.N(Wet_den,delta_wet)
    Wet_den.plot(label='Density (sat.)',lw=2,color='b')
    Wet_denMC.plot(hist=True,label='Density (sat.) (MC)',color='g')
    plt.legend()
    plt.xlabel('Density (g/cm3)')
    plt.savefig('./'+well+'/'+well+'_Saturated_Den.pdf',bbox_inches='tight')
    plt.show(block=True)

    #Mineral Density
    plt.figure('Mineral Density')
    Min_den=mc.N(Min_den,delta_min)
    Min_den.plot(label='Density (min.)',lw=2,color='b')
    Min_denMC.plot(hist=True,label='Density (min.) (MC)',color='g')
    plt.legend()
    plt.xlabel('Density (g/cm3)')
    plt.savefig('./'+well+'/'+well+'_Mineral_Den.pdf',bbox_inches='tight')
    plt.show(block=True)

    #Porosity
    plt.figure('Porosity')
    por=mc.N(por,delta_por)
    por.plot(label='Porosity (%)',lw=2,color='b')
    Por_MC.plot(hist=True,label='Porosity (%) (MC)',color='g')
    plt.legend()
    plt.xlabel('Porosity (%)')
    plt.savefig('./'+well+'/'+well+'_Porosity.pdf',bbox_inches='tight')
    plt.show(block=True)
