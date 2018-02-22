""" A python script to read in a CSV file and plot data for lab XXX"""

# unlike matlab, in python you have to call the functions/modules you
# are going to use:
import csv # this is the module that can read, write and manipulate csv files
import numpy as np # the numerical capabilities of python
import matplotlib.pyplot as plt # for plotting 
import glob
import RockPhysics as rp
#from obspy.core import Trace, Stats, Stream
from scipy import signal

# Set default font size for plots:
font = {'size'   : 24}
plt.rc('font',**font)

def main():

    # put all file names to read in a list. All files for one rock and one fluid (all pressures) should be in ONE folder as csv files
    files = sorted(glob.glob('*.csv'),key=lambda f: int(filter(str.isdigit, f)))
    
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
    
    indices=np.where((phases=='PP')& (cycles=='d2')& (states=='sat500'))[0]    
    
    
    # define values:
    Pp=500 #psi
    labels = []
    offset=0
    fig, ax = plt.subplots(figsize=(12,8))

    # read and plot data in each file:
    for i in range(0,len(indices)):
        file=files[indices[i]]
        # read in time/amp:
        t,V=rp.load_csv(file,0)#readcsvdata(file)
        
        # rather arbitrary filter with params to tune:
        #b, a = signal.butter(2, 0.005, 'low')
        #V = signal.filtfilt(b, a, V)

        # normalises the traces:
        V= 1*V/np.max(V)
        
        #Load the picks
        mean,std_error=np.loadtxt(''+phases[indices[i]]+'_time_std_'+states[indices[0]]+'_'+np.str(np.int(Pc[indices[i]]))+'.out')        
        # plot in \mu s:
        plt.plot(t,V+offset,'k')
        p=np.min(np.where(t>=mean)) #y-position on the trace where the picking was made
        plt.scatter(mean,V[p]+offset,s=60,c='r',zorder=2)
        plt.errorbar(mean,V[p]+offset,xerr=2*std_error,fmt='o',c='r',lw=2,zorder=1)

        # offset traces vertically:
        offset = offset +1

        # create a custom y-label, based on pressure in file name:
        labels.append(Pc[indices[i]]-Pp)  # 32:35 for 142,
        
              
    # this is very manual to set the pressure values on the y-axis. It
    # depends crucially on the offset value and the files/pressures
    # recorded:
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Differential pressure (psi)')
    plt.title(phases[indices[i]]+' Normalised waveforms from rock '+name)
    plt.grid()
    plt.axis('tight')
    t1=11
    t2=17    
    plt.xlim(np.min(t1),np.max(t2))
    plt.savefig(name+'_'+phases[indices[i]]+'_picks.pdf')
    plt.show()

def readcsvdata(filename):
    """ function to read csv. Input filename, output floats time 
    and amplitude """
    csv_in = open(filename,'rU')
    myreader = csv.reader(csv_in,delimiter = ',')
    
    ###### skip the lines with header info: ##########################
    numberofheaderlines=0
    for i in range(numberofheaderlines):
        next(myreader) # skipping over the headers (change range to
                       # higher than 0 to skip lines)
    ##################################################################
    ####### reading in the csv values as strings, then converting to
    ####### floats: ##################################################
    x =[]
    y =[]
    for row in myreader:
        x.append(row[0])
        y.append(row[1])
    x=np.array([float(i) for i in x])
    y=np.array([float(i) for i in y])
    return x,y

# this will actually run the code:
if __name__ == '__main__':
    main()
    
