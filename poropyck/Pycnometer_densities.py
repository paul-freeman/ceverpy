# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 15:24:32 2017

@author: edur409
"""
def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def histogram(v,title): 
    x=np.linspace(np.min(v)-2*np.std(v),np.max(v)+2*np.std(v),50)
    y=normal.pdf(x,np.mean(v),np.std(v)) #range, mean, std
    plt.figure(title)
    plt.hist(v,normed=True,color='b')
    plt.plot(x,y,'r',lw=2)
    return


import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import mcerp3 as mc
import mcerp3.umath as umath 
#import urllib2
from urllib.request import urlopen
#from bs4 import BeautifulSoup
from scipy.stats import norm as normal

def pycnometer_densities():
    font = {'size'   : 18}
    plt.rc('font',**font)

    F=np.loadtxt('Samples_Depth.txt',dtype={'names': ('sample', 'depth'),'formats': ('S20', 'f4')})

    for i in range(0,len(F)):
        sample=F[i][0].decode('UTF-8')
        #filename='untitled-2'
        response = urlopen('file:./Pycnometer_data/'+sample+'_10cycles/Pycnometer_'+sample+'.html')
        #urllib.urlretrieve('http://130.216.47.208/3854/'+filename+'.html', ''+filename+'.html')
        html = response.read()
        A=html.splitlines() #Transform the html into lines
        #Number of purges
        n=np.int(A[9].split()[3])
        V=[]
        dV=[]
        Temp=[]
        
        c=0 #Count the number of cycles
        for i in range(0,len(A)-2): 
            if len(A[i])>0: #if the line is not empty
                if RepresentsInt(A[i].split()[0]): #read data as long as first character in line is a number
                    v=np.float(A[i].split()[1])
                    dv=np.float(A[i].split()[2])
                    T=np.float(A[i].split()[6])
                    V=np.append(V,v)
                    dV=np.append(dV,dv)
                    Temp=np.append(Temp,T)
                    c=c+1
            
        #Plots  
        histogram(V,'Volume')
        plt.xlabel('Volume ($cm^3$)'),plt.ylabel('Normalized frequency')
        plt.title('Cycles = '+np.str(c)+' , V = '+np.str(np.round(np.mean(V),decimals=4))+'$\pm$'+np.str(np.round(np.std(V),decimals=4))+' $cm^3$')
        Vmc=mc.N(np.mean(V),np.std(V)) 
        Vmin,Vmax=Vmc.percentile((0.025,0.975)) #95% Confidence Interval  
        plt.savefig('./'+sample+'/Volume_'+sample+'.pdf',bbox_inches='tight')
        
        histogram(Temp,'Temperature')
        plt.xlabel('Temperature ($^\circ C$)'),plt.ylabel('Normalized frequency')
        plt.title('Cycles = '+np.str(c)+' , T = '+np.str(np.round(np.mean(Temp),decimals=1))+'$\pm$'+np.str(np.round(np.std(Temp),decimals=1))+'$^\circ C$')
        plt.savefig('./'+sample+'/Temperature_'+sample+'.pdf',bbox_inches='tight')
        
        plt.figure('V vs. T')
        plt.scatter(Temp,V,s=80,c='b')
        plt.errorbar(Temp,V,yerr=dV,c='b',fmt='o')
        plt.xlabel('Temperature ($^\circ C$)'),plt.ylabel('Volume ($cm^3$)')
        plt.savefig('./'+sample+'/VolvsTemp_'+sample+'.pdf',bbox_inches='tight')
        
        M=np.loadtxt('./Pycnometer_data/'+sample+'_10cycles/Mass_'+sample+'.txt')#65.0216#55.4415#55.6390#52.2243#54.4906#55.9833#35.8836#19.7361#31.5223#33.8521#48.3861#27.62#48.3864#39.5264#51.0950#31.3210#51.0987
        #np.savetxt('Mass_'+sample+'.txt',M,delimiter=' ',fmt='%1.4e')
        RHO=M/V
        histogram(RHO,'Density')
        plt.xlabel('Density ($g/cm^3$)'),plt.ylabel('Normalized frequency')
        plt.title('Cycles = '+np.str(c)+' , rho = '+np.str(np.round(np.mean(RHO),decimals=4))+'$\pm$'+np.str(np.round(np.std(RHO),decimals=4))+' $g/cm^3$')
        plt.savefig('./'+sample+'/Density_Pycnometer_'+sample+'.pdf',bbox_inches='tight')
        
           
        dM=0.001 #error of the scale
        rho=M/np.mean(V)
        drho=rho*np.sqrt((dM/M)**2+(np.std(V)/np.mean(V)**2))
        x=np.linspace(rho-3*drho,rho+3*drho,50)
        y=normal.pdf(x,rho,drho) #range, mean, std
        plt.figure('Density_Evert')
        #plt.hist(v,normed=True,color='b')
        plt.plot(x,y,'r',lw=2)
        plt.xlabel('Density ($g/cm^3$)'),plt.ylabel('Normalized frequency')
        plt.title('$rho$ = '+np.str(np.round(rho,decimals=2))+'$\pm$'+np.str(np.round(drho,decimals=2))+'$g/cm^3$')
        plt.savefig('./'+sample+'/Density_Pycnometer2_'+sample+'.pdf',bbox_inches='tight')
        
         
            
        #SAVE THE VALUES IN TEXT FILES
        np.savetxt('./'+sample+'/Volume_Pycnometer.out', [np.mean(V),np.std(V),Vmin,Vmc.mean,Vmax], fmt='%1.2f',delimiter=',',header='V (cm^3), Std_V (cm^3), V_0025, VMC_mean, V_0975') 
        #if dRHO is below resolution  
        if np.std(RHO)<0.01:
            np.savetxt('./'+sample+'/Density_Pycnometer.out', [np.mean(RHO),0.01,rho,drho], fmt='%1.2f',delimiter=',',header='Rho_pyc (min.), Std_Rho_pyc (min.), Rho, dRho')
        else:     
            np.savetxt('./'+sample+'/Density_Pycnometer.out', [np.mean(RHO),np.std(RHO),rho,drho], fmt='%1.2f',delimiter=',',header='Rho_pyc (min.), Std_Rho_pyc (min.), Rho, dRho')
        plt.close('all')

    #def Vol_porosity(sample):
    #    #Load the geometrical parmeters for volume 
    #    A=np.loadtxt('./Pycnometer_data/'+sample+'_10cycles/Vol_params.txt',skiprows=1)
    #    L=np.mean(A[:,0]) #Average length of sample
    #    dL=np.std(A[:,0]) #standard deviation of length of sample
    #    D=np.mean(A[:,1]) #Average length of sample
    #    dD=np.std(A[:,1]) #standard deviation of diameter of sample
    #    #Volume of sample and its error
    #    V=0.25*np.pi*L*D**2
    #    dV=V*np.sqrt(4*(dD/D)**2+(dL/L)**2)
    #    #Load the eskeletal volumes of the pycnometer
    #    B=np.loadtxt('./'+sample+'/Volume_Pycnometer.out')
    #    Vp=B[0] #Volume of pycnometer
    #    dVp=B[1] #error in volume of pycnometer
    #    #Calculate the porosity 
    #    phi=100*(V-Vp)/V
    #    dphi=phi*np.sqrt((dVp/(V-Vp))**2+(((1/(V-Vp))-1/V)*dV)**2)
    #    np.savetxt('./'+sample+'/Porosity_Pycnometer.out', [phi,dphi], fmt='%1.2f',delimiter=',',header='phi (%) , dphi (%)') 
    #    return V,Vp
    #    
    #    
    #    
    #sample='NM11_2087_4C'
    #V,Vp=Vol_porosity(sample)
    #A=np.loadtxt('./Pycnometer_data/'+sample+'_10cycles/Vol_params.txt',skiprows=1)

        
    #file.close()    
    ##Number of cycles
    #n=np.int(A[8].split()[3])
    #V=[]
    #dV=[]
    #Temp=[]
    #for i in range(15,15+n):
    #    v=np.float(A[i].split()[1])
    #    dv=np.float(A[i].split()[2])
    #    T=np.float(A[i].split()[6])
    #    V=np.append(V,v)
    #    dV=np.append(dV,dv)
    #    Temp=np.append(Temp,T)
    #    
    #histogram(V,'Volume')
    #plt.xlabel('Volume ($cm^3$)'),plt.ylabel('Normalized frequency')
    #plt.title('Cycles = '+np.str(n)+' , V = '+np.str(np.round(np.mean(V),decimals=4))+'$\pm$'+np.str(np.round(np.std(V),decimals=4))+' $cm^3$')
    #Vmc=mc.N(np.mean(V),np.std(V)) 
    #Vmin,Vmax=Vmc.percentile((0.025,0.975))   
    #
    #histogram(Temp,'Temperature')
    #plt.xlabel('Temperature ($^\circ C$)'),plt.ylabel('Normalized frequency')
    #plt.title('Cycles = '+np.str(n)+' , T = '+np.str(np.round(np.mean(Temp),decimals=1))+'$\pm$'+np.str(np.round(np.std(Temp),decimals=1))+'$^\circ C$')
    #
    #plt.figure('V vs. T')
    #plt.scatter(Temp,V,s=80,c='b')
    #plt.errorbar(Temp,V,yerr=dV,c='b',fmt='o')
    #plt.xlabel('Temperature ($^\circ C$)'),plt.ylabel('Volume ($cm^3$)')
