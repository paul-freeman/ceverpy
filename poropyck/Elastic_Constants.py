# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 12:40:08 2017

@author: leo
"""
#Import the Rock Physics subroutines
from pkg_resources import resource_filename
from . import RockPhysics as rp #Loads all the necessary files
import numpy as np
import matplotlib.pyplot as plt
import mcerp3 as mc

def elastic_constants():
    #Font size for plots
    font= {'size' : 18}
    plt.rc('font',**font)

    filename = resource_filename('poropyck', 'Samples_Depth.txt')
    W=np.loadtxt(filename,dtype={'names': ('sample', 'depth'),'formats': ('S20', 'f4')})
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

    #Vpd=[]
    #dVpd=[]
    #Vps=[]
    #dVps=[]
    suffix='Arch' #Suffix for the files

    Mud=[]
    dMud=[]
    Mus=[]
    dMus=[]
    Kd=[]
    dKd=[]
    Ks=[]
    dKs=[]
    Ed=[]
    dEd=[]
    Es=[]
    dEs=[]
    Pd=[]
    dPd=[]
    Ps=[]
    dPs=[]


    for i in range(0,len(W)):
        well=W[i][0].decode('UTF-8') #Well
        #Read the density files
        A=np.loadtxt('./'+well+'/Dry_density2.out')
        B=np.loadtxt('./'+well+'/Sat_density2.out')
        rhod=A[0]
        drhod=A[1]
        rhos=B[0]
        drhos=B[1]
        #Create the distributions for the MC simulations
        RHOD=mc.N(rhod,drhod)
        RHOS=mc.N(rhos,drhos)
        print((''+well+' ($\rho_d$): ', RHOD.mean,RHOD.std))
        print((''+well+' ($\rho_s$): ', RHOS.mean,RHOS.std))    
        
        #MU    
        Vsd=np.loadtxt('./'+well+'/'+well+'d/Velocities_S_dry_both.out',skiprows=1)[0]
        dVsd=np.loadtxt('./'+well+'/'+well+'d/Velocities_S_dry_both.out',skiprows=1)[1]
        Vss=np.loadtxt('./'+well+'/'+well+'s/Velocities_S_sat_both.out',skiprows=1)[0]
        dVss=np.loadtxt('./'+well+'/'+well+'s/Velocities_S_sat_both.out',skiprows=1)[1]
        Velsd=mc.N(Vsd,dVsd)    
        Velss=mc.N(Vss,dVss)
        
        #Shear modulus    
        MUD=1e-6*RHOD*Velsd**2 #Values in GPa
        MUS=1e-6*RHOD*Velss**2 #Values in GPa    
        np.savetxt('./'+well+'/'+well+'d/Shear_dry_'+suffix+'.out',[MUD.mean,MUD.std],header='MUd (GPa) , dMUd (GPa)')    
        np.savetxt('./'+well+'/'+well+'s/Shear_sat_'+suffix+'.out',[MUS.mean,MUS.std],header='MUs (GPa) , dMUs (GPa)')        
        
        #P-wave velocities
        #Dry
        Vd=np.loadtxt('./'+well+'/'+well+'d/Velocities_P_dry_both.out',skiprows=1)[0]
        dVd=np.loadtxt('./'+well+'/'+well+'d/Velocities_P_dry_both.out',skiprows=1)[1]
        Veld=mc.N(Vd,dVd) #Velocity of dry rock
        #Saturated
        Vs=np.loadtxt('./'+well+'/'+well+'s/Velocities_P_sat_both.out',skiprows=1)[0]
        dVs=np.loadtxt('./'+well+'/'+well+'s/Velocities_P_sat_both.out',skiprows=1)[1]
        Vels=mc.N(Vs,dVs) #Velocity of dry rock
        
        print((''+well+' (Vpd): ', Veld.mean,Veld.std)) 
        print((''+well+' (Vps): ', Vels.mean,Vels.std))
        
        
        C=np.loadtxt('./'+well+'/'+well+'d/Shear_dry_'+suffix+'.out')
        D=np.loadtxt('./'+well+'/'+well+'s/Shear_sat_'+suffix+'.out')
        mud=C[0]
        dmud=C[1]
        mus=D[0]
        dmus=D[1]
        
        #Create distributions for Mud and Mus
        MUD=mc.N(mud,dmud)
        MUS=mc.N(mus,dmus)
        print((MUD.mean,MUS.std))
        
        #Bulk modulus
        KD=1e-6*RHOD*Veld**2-(4.0/3)*MUD
        KS=1e-6*RHOS*Vels**2-(4.0/3)*MUS
        print((KD.mean,KS.std))
        np.savetxt('./'+well+'/'+well+'d/Bulk_dry_'+suffix+'.out',[KD.mean,KD.std],header='Kd (GPa) , dKd (GPa)')    
        np.savetxt('./'+well+'/'+well+'s/Bulk_sat_'+suffix+'.out',[KS.mean,KS.std],header='Ks (GPa) , dKs (GPa)')
        
        #Poisson's Ratio  
        PD=(3*KD-2*MUD)/(2*(3*KD+MUD))
        PS=(3*KS-2*MUS)/(2*(3*KS+MUS))
        print((PD.mean,PS.std)) 
        np.savetxt('./'+well+'/'+well+'d/Poisson_dry_'+suffix+'.out',[PD.mean,PD.std],header='Pd  , dPd ')    
        np.savetxt('./'+well+'/'+well+'s/Poisson_sat_'+suffix+'.out',[PS.mean,PS.std],header='Ps  , dPs ')
        
        #Young Modulus
        ED=9*KD*MUD/(3*KD+MUD)
        ES=9*KS*MUS/(3*KS+MUS)
        print((ED.mean,ES.std))
        np.savetxt('./'+well+'/'+well+'d/Young_dry_'+suffix+'.out',[ED.mean,ED.std],header='Ed (GPa) , dEd (GPa)')    
        np.savetxt('./'+well+'/'+well+'s/Young_sat_'+suffix+'.out',[ES.mean,ES.std],header='Es (GPa) , dEs (GPa)')
        
        #Figures
        plt.figure('K (sat.)')
        KS.plot(label='$K$ (sat.)',lw=2,color='b')
        KS.plot(hist=True,label='$K$ (sat.) (MC)',color='g')
        plt.legend()
        plt.xlabel('$K$ (GPa)')
        plt.savefig('./'+well+'/'+well+'s/Bulk_sat_'+suffix+'.pdf',bbox_inches='tight')
        plt.close()    
        
        plt.figure('Bulk (dry)')
        KD.plot(label='$K$ (dry)',lw=2,color='b')
        KD.plot(hist=True,label='$K$ (dry) (MC)',color='g')
        plt.legend()
        plt.xlabel('$K$ (GPa)')
        plt.savefig('./'+well+'/'+well+'d/Bulk_dry_'+suffix+'.pdf',bbox_inches='tight')
        plt.close()    
        
        #Poisson's ratio figures
        plt.figure('Poisson (sat.)')
        PS.plot(label='$v$ (sat.)',lw=2,color='b')
        PS.plot(hist=True,label='$v$ (sat.) (MC)',color='g')
        plt.legend()
        plt.xlabel('$v$')
        plt.savefig('./'+well+'/'+well+'s/Poisson_sat_'+suffix+'.pdf',bbox_inches='tight')
        plt.close()    
        
        plt.figure('Poisson (dry.)')
        PD.plot(label='$v$ (dry.)',lw=2,color='b')
        PD.plot(hist=True,label='$v$ (dry.) (MC)',color='g')
        plt.legend()
        plt.xlabel('$v$')
        plt.savefig('./'+well+'/'+well+'d/Poisson_dry_'+suffix+'.pdf',bbox_inches='tight')
        plt.close()
        
        #Young's modulus figures
        plt.figure('E (sat.)')
        ES.plot(label='$E$ (sat.)',lw=2,color='b')
        ES.plot(hist=True,label='$E$ (sat.) (MC)',color='g')
        plt.legend()
        plt.xlabel('$E$ (GPa)')
        plt.savefig('./'+well+'/'+well+'s/Young_sat_'+suffix+'.pdf',bbox_inches='tight')
        plt.close()    
        
        plt.figure('E (dry.)')
        ED.plot(label='$E$ (dry.)',lw=2,color='b')
        ED.plot(hist=True,label='$E$ (dry.) (MC)',color='g')
        plt.legend()
        plt.xlabel('$E$ (GPa)')
        plt.savefig('./'+well+'/'+well+'d/Young_dry_'+suffix+'.pdf',bbox_inches='tight')
        plt.close()
        
        #Save values
        Mud=np.append(Mud,np.loadtxt('./'+well+'/'+well+'d''/Shear_dry_Arch.out')[0])
        dMud=np.append(dMud,np.loadtxt('./'+well+'/'+well+'d''/Shear_dry_Arch.out')[1])
        Mus=np.append(Mus,np.loadtxt('./'+well+'/'+well+'s''/Shear_sat_Arch.out')[0])
        dMus=np.append(dMus,np.loadtxt('./'+well+'/'+well+'s''/Shear_sat_Arch.out')[1])
        Kd=np.append(Kd,np.loadtxt('./'+well+'/'+well+'d''/Bulk_dry_Arch.out')[0])
        dKd=np.append(dKd,np.loadtxt('./'+well+'/'+well+'d''/Bulk_dry_Arch.out')[1])
        Ks=np.append(Ks,np.loadtxt('./'+well+'/'+well+'s''/Bulk_sat_Arch.out')[0])
        dKs=np.append(dKs,np.loadtxt('./'+well+'/'+well+'s''/Bulk_sat_Arch.out')[1])
        Ed=np.append(Ed,np.loadtxt('./'+well+'/'+well+'d''/Young_dry_Arch.out')[0])
        dEd=np.append(dEd,np.loadtxt('./'+well+'/'+well+'d''/Young_dry_Arch.out')[1])
        Es=np.append(Es,np.loadtxt('./'+well+'/'+well+'s''/Young_sat_Arch.out')[0])
        dEs=np.append(dEs,np.loadtxt('./'+well+'/'+well+'s''/Young_sat_Arch.out')[1])
        Pd=np.append(Pd,np.loadtxt('./'+well+'/'+well+'d''/Poisson_dry_Arch.out')[0])
        dPd=np.append(dPd,np.loadtxt('./'+well+'/'+well+'d''/Poisson_dry_Arch.out')[1])
        Ps=np.append(Ps,np.loadtxt('./'+well+'/'+well+'s''/Poisson_sat_Arch.out')[0])
        dPs=np.append(dPs,np.loadtxt('./'+well+'/'+well+'s''/Poisson_sat_Arch.out')[1])

    ##Print statement
    #print('%4.2f +/- %4.2f ' % (Vpd[0],dVpd[0]))
    #s='%4.2f +/- %4.2f ' % (Vpd[0],dVpd[0]) #Print the same into a string
    #

        
        
    filename='Table_Constants.txt'
    F=open(filename,'w')
    F.write('Well, $\mu$ (dry), $\mu$ (sat), $K$ (dry), $K$ (sat.), $E$ (dry), $E$ (sat.) , $\v$ (dry), $\v$ (sat.) \n')
    S=[]
    for i in range(0,len(W)):
        name=well
        s1='%3.2f +/- %3.2f ' % (Mud[i],dMud[i]) #Print the same into a string
        s2='%3.2f +/- %3.2f ' % (Mus[i],dMus[i]) #Print the same into a string
        s3='%3.2f +/- %3.2f ' % (Kd[i],dKd[i]) #Print the same into a string
        s4='%3.2f +/- %3.2f ' % (Ks[i],dKs[i]) #Print the same into a string
        s5='%3.2f +/- %3.2f ' % (Ed[i],dEd[i]) #Print the same into a string
        s6='%3.2f +/- %3.2f ' % (Es[i],dEs[i]) #Print the same into a string
        s7='%3.2f +/- %3.2f ' % (Pd[i],dPd[i]) #Print the same into a string
        s8='%3.2f +/- %3.2f ' % (Ps[i],dPs[i]) #Print the same into a string    
        F.write('%s, %s, %s, %s, %s, %s, %s, %s, %s \n'% (name, s1 ,  s2,  s3,  s4,  s5,  s6, s7, s8))
    F.close()

    size=40    
    fig = plt.figure('Poisson (dry) vs. Poisson (sat) ')
    ax = fig.add_subplot(111)   
    for i in range(0,len(W)):
        ax.scatter(Pd[i],Ps[i],s=size,c='b')
        #ax.text(Pd[i],Ps[i], well, fontsize=font) 
        #ax.text(Vpsc[i],Vssc[i], A[i][0], fontsize=font) 
    ax.errorbar(Pd, Ps, yerr=[2*dPs, 2*dPs], xerr=[2*dPd, 2*dPd], fmt='o',color='b')
    plt.xlabel('Poisson (dry)'),plt.ylabel('Poisson (sat.)')        
    plt.axis('tight')
    #plt.colorbar()
    plt.grid('on')  
