"""Poropyck common helper functions"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.signal import hilbert
from scipy.signal import blackman
import mcerp3 as mc
from . import rock_physics as rp
from .dtw import dtw # pylint: disable=no-name-in-module


COORDS = []


def crosscorr_lags(a_value, b_value):
    """crosscorr_lags"""
    c_value = np.correlate(a_value, b_value, mode='full')
    lags = np.linspace(-len(a_value), len(a_value), len(c_value))
    return lags, c_value


def onpick(event):
    """Subroutine to pick values from the active plot"""
    global COORDS
    thisline = event.artist
    xdata = thisline.get_xdata()
    ydata = thisline.get_ydata()
    ind = event.ind
    points = np.array([xdata[ind], ydata[ind]])
    print('onpick points:', points)
    COORDS.append((xdata[ind], ydata[ind]))
    return COORDS


# "Magic numbers"
X_MIN = 0
X_MAX = 40
X_DELTA = 0.01
DELTA = 0.0

A1_VALUES = [1.0, 1.5, 2.0, 3.0, 4.0, 2.0, 0.9]
F1_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
F2_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]

S_START = 500
S_END = 1500

T_START = 0
T_END = 1000
T_OFFSET = 10.0
T_OFFSET_2 = 2.0

def simple_harmonics(a2_values, cross_correlation, use_manager, print_time_lag, labels):
    """Define two functions to compared composed of simple harmonics"""
    global COORDS
    a1_values = A1_VALUES
    f1_values = F1_VALUES
    # a2_values = a2_values # these are passed into the function
    f2_values = F2_VALUES
    t_value = np.arange(X_MIN, X_MAX, X_DELTA)
    s1_values = np.zeros(len(t_value))
    s2_values = np.zeros(len(t_value))
    for i, values in enumerate(zip(a1_values, a2_values, f1_values, f2_values)):
        a1_value, a2_value, f1_value, f2_value = values
        s1_values += (
            a1_value * np.sin(2 * np.pi * f1_value * t_value)
            + a1_value * np.cos(2 * np.pi * f1_value * t_value)
            )
        s2_values += (
            a2_value * np.sin(2 * np.pi * f2_value * (t_value + DELTA))
            + a2_value * np.cos(2 * np.pi * f2_value * (t_value + DELTA))
            )
    s1_final = s1_values[S_START:S_END] * blackman(len(s1_values[S_START:S_END]))
    s2_final = s2_values[S_START:S_END] * blackman(len(s2_values[S_START:S_END]))

    well = input("Type name of sample (e.g. 'Dummy'): \n")

    time_sd = t_value[T_START:T_END] + T_OFFSET
    sd_values = s1_final

    print('Choose the beginning and end that you want to compare')
    #Choose the two extremes to compare
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('click on points')
    COORDS = []
    ax.plot(time_sd, sd_values, '-', picker=5)
    plt.xlabel(r'Time ($\mu$s)')
    plt.grid(color='0.5')
    _ = fig.canvas.mpl_connect('pick_event', onpick)
    #Work around for plt.show(block=True) in IPython using Spyder
    plt.show(block=True)#plt.show(block=True)
    #Pick the times
    t_id = np.min(COORDS[0][0])
    t_fd = np.min(COORDS[1][0])

    #Indices to crop the dry sample time vector
    i_d = np.min(np.where(time_sd >= t_id))
    f_d = np.min(np.where(time_sd >= t_fd))

    #1.2) READ THE SATURATED SAMPLE WAVEFORM SECOND
    #time_ss,T,Ss=rp.plot_csv('./'+well+'/'+well+'s/tek0001ALL.csv',0)
    time_ss = t_value[T_START:T_END] + T_OFFSET + T_OFFSET_2
    Ss = s2_final


    print('Choose the beginning and end that you want to compare')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('click on points')
    COORDS = []
    ax.plot(time_ss, Ss, '-', picker=5)
    plt.xlabel(r'Time ($\mu$s)')
    plt.grid(color='0.5')
    fig.canvas.mpl_connect('pick_event', onpick)
    plt.show(block=True)

    t_is = np.min(COORDS[0][0])
    t_fs = np.min(COORDS[1][0])

    #Indices to crop for the saturated sample time vector
    i_s = np.min(np.where(time_ss >= t_is))
    f_s = np.min(np.where(time_ss >= t_fs))

    #2) DO THE COMPARISON USING THE DTW ALGORITHM
    #THE DYNAMIC TIME WARPING ALGORITHM BEGINS HERE!!!!
    idxt = time_sd[i_d:f_d] #Indices of the template function (dry sample waveform)
    idxq = time_ss[i_s:f_s] #Indices of the query function (saturated sample waveform)

    #2.1) CHOOSE INFORMATION TO COMPARE (WAVEFORMS OR ENVELOPES)
    template = sd_values[i_d:f_d]/np.max(np.abs(sd_values[i_d:f_d]))
    query = Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s]))
    templatea = (np.angle(hilbert(sd_values[i_d:f_d]/np.max(np.abs(sd_values[i_d:f_d])))))/np.pi
    querya = (np.angle(hilbert(Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s])))))/np.pi
    templateh = np.abs(hilbert(sd_values[i_d:f_d]/np.max(np.abs(sd_values[i_d:f_d]))))
    queryh = np.abs(hilbert(Ss[i_s:f_s]/np.max(np.abs(Ss[i_s:f_s]))))
    suffix = 'both' #Suffix for the files

    # Calculate the alignment vector and corresponding distance (DTW)
    # Delete costs matrix because it can be quite large - and we don't use it.
    dist, indices1, indices2, costs = dtw(query, template)
    del costs
    disth, indices1h, indices2h, costsh = dtw(queryh, templateh)
    del costsh
    dista, indices1a, indices2a, costsa = dtw(querya, templatea)
    del costsa

    # The lower the distance of alignment, the better the match
    print('Distance of the DTW algorithm (waveform): {:.3f}'.format(dist))
    print('Distance of the DTW algorithm (phase): {:.3f}'.format(dista))
    print('Distance of the DTW algorithm (envelope): {:.3f}'.format(disth))

    #Indices of alignment
    I1 = indices1
    I2 = indices2
    I1h = indices1h
    I2h = indices2h
    I1a = indices1a
    I2a = indices2a

    print('Close the figures to continue running the code...')
    #3) PLOT THE OUTPUTS OF THE DTW ALGORITHM
    #PLOTTING THE OUTPUTS OF THE DTW ALGORITHM BEGINS HERE!!!
    #Plots
    #Plot #1: points of match between the waveforms
    plt.figure('Points of match')
    plt.plot(idxt, template, label='Pd', c='b')
    plt.plot(idxq, query, label='Ps', c='g')
    plt.axis('tight')
    plt.legend()
    plt.xlabel(r'time ($\mu$s)')
    #Plot the matching points
    for i in np.arange(0, len(I1)-1, 20):
        plt.plot(
            [idxt[np.int(I2[i]-1)], idxq[np.int(I1[i]-1)]],
            [template[np.int(I2[i]-1)], query[np.int(I1[i]-1)]],
            'r-',
            lw=0.5
            )
    #Save the figure
    if use_manager:
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
    plt.show()
    plt.savefig(
        './{0}/{0}d/DTW_Match_Pwaves_{1}.pdf'.format(well, suffix),
        bbox_inches='tight'
        )
    plt.savefig(
        './{0}/{0}s/DTW_Match_Pwaves_{1}.pdf'.format(well, suffix),
        bbox_inches='tight'
        )
    plt.show(block=True)#plt.show(block=True)

    #Plot #2 (mostly unnecessary): Compare a function over the other
    #plt.figure('One over the other')
    qo_ = []
    qoa = []
    qoh = []
    to_ = []
    toa = []
    toh = []
    idxto = []
    idxtoa = []
    idxtoh = []
    idxqo = []
    idxqoa = []
    idxqoh = []
    for i, _ in enumerate(I1):
        # time vector arranged by index
        idxto = np.append(idxto, idxt[np.int(I2[i])-1])
        # time vector arranged by index
        idxqo = np.append(idxqo, idxq[np.int(I1[i])-1])
        # Query function sampled by the index
        qo_ = np.append(qo_, query[np.int(I1[i])-1])
        # Template function sampled by the index
        to_ = np.append(to_, template[np.int(I2[i])-1])
    for i, _ in enumerate(I1h):
        # time vector arranged by index
        idxtoh = np.append(idxtoh, idxt[np.int(I2h[i])-1])
        # time vector arranged by index
        idxqoh = np.append(idxqoh, idxq[np.int(I1h[i])-1])
        # Query function sampled by the index
        qoh = np.append(qoh, queryh[np.int(I1h[i])-1])
        # Template function sampled by the index
        toh = np.append(toh, templateh[np.int(I2h[i])-1])
    for i, _ in enumerate(I1a):
        # time vector arranged by index
        idxtoa = np.append(idxtoa, idxt[np.int(I2a[i])-1])
        # time vector arranged by index
        idxqoa = np.append(idxqoa, idxq[np.int(I1a[i])-1])
        # Query function sampled by the index
        qoa = np.append(qoa, querya[np.int(I1a[i])-1])
        # Template function sampled by the index
        toa = np.append(toa, templatea[np.int(I2a[i])-1])
    #Plot #3: A Summary plot with the alignment function linking the two waveforms
    # definitions for the axes
    left, width = 0.12, 0.60
    bottom, height = 0.08, 0.60
    bottom_h = 0.16 + width
    left_h = left + 0.27
    rect_plot = [left_h, bottom, width, height]
    rect_x = [left_h, bottom_h, width, 0.2]
    rect_y = [left, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure('Summary', figsize=(10, 10))

    #axplot = fig.add_subplot(111)
    axplot = plt.axes(rect_plot)
    axx = plt.axes(rect_x)
    axy = plt.axes(rect_y)

    # Plot the matrix
    axplot.plot(idxqo, idxto, 'k', picker=5, lw=2)
    axplot.plot(idxqoa, idxtoa, 'm', picker=5, lw=2)
    axplot.plot(idxqoh, idxtoh, '--', color='m', picker=5, lw=2)
    # Give same time scale as template
    axplot.axis([t_id, t_fd, t_id, t_fd])
    axplot.grid(color='0.5')

    # Indices in P correspond to locations where S-times saturated are larger
    # than S-times dry
    # Find the locations where time of S sat. is larger than time of S dry
    P = np.where(idxqo >= idxto)
    Ph = np.where(idxqoh >= idxtoh)
    x1 = np.linspace(0, np.max([time_sd, time_ss]), 10)
    y1 = x1
    axplot.plot(x1, y1, 'g', lw=2) #plot the 1:1 line of match

    axplot.tick_params(axis='both', which='major', labelsize=18)

    # Plot time serie horizontal
    axx.plot(idxq,query, '-', color='b', lw=2)
    axx.plot(idxq,queryh, '--', color='m', lw=2) #The envelope
    axx.plot(idxq,querya, '-', color='m', lw=2)
    axx.grid(color='0.5')
    axx.axis([t_id, t_fd, 1.1 * np.min(query), 1.1 * np.max(queryh)])
    axx.tick_params(axis='both', which='major', labelsize=18)


    # Plot time serie vertical
    axy.plot(template, idxt, '-', color='r', lw=2)
    axy.plot(templateh, idxt, '--', color='m', lw=2)
    axy.plot(templatea, idxt, '-', color='m', lw=2)
    axy.grid(color='0.5')
    axy.axis([1.1 * np.min(template), 1.1 * np.max(templateh), t_id, t_fd])
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
    COORDS = []
    fig.canvas.mpl_connect('pick_event', onpick)

    plt.show(block=True)#plt.show(block=True)

    #ESTIMATION OF THE VELOCITIES FROM THE PICKED TIMES
    ##
    #Load the caliper lengths
    lengths = np.loadtxt('./{0}/{0}_Lengths.txt'.format(well), delimiter=',')
    L, StdL = rp.Length_sample(lengths)

    #Extract the times from the DTW picking
    tD, tS = rp.S_values(COORDS)
    np.savetxt('./{0}/{0}d/P_time_Picks_dry_{1}.out'.format(well, suffix), tD) #save the picked times
    np.savetxt('./{0}/{0}s/P_time_Picks_sat_{1}.out'.format(well, suffix), tS) #save the picked times
    #Save the histograms
    rp.histogram(tD,'P-arrival (dry)')
    plt.xlabel('Time ($\mu$s)')
    plt.savefig('./'+well+'/'+well+'d/P_times_hist_dry_'+suffix+'.pdf',bbox_inches='tight')
    plt.show(block=True)
    rp.histogram(tS,'P-arrival (sat.)')
    plt.xlabel('Time ($\mu$s)')
    #save the figure
    if use_manager:
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
    plt.show()
    plt.savefig('./'+well+'/'+well+'s/P_times_hist_sat_'+suffix+'.pdf',bbox_inches='tight')
    plt.show(block=True)

    #Time differece or lag from the DTW algorithm
    if print_time_lag:
        print('Time Lag (DTW): ', np.round(np.mean(tS)-np.mean(tD), decimals=2))

    #CALCULATE THE VELOCITIES FROM THE TIMES PICKED FOR THE DRY SAMPLE
    Vd,dVd,Vd_down,Vdmc,Vd_up=rp.Velocity_S(L,np.mean(tD),StdL,np.std(tD)) #Velocity S is more general than Velocity_P
    #Save the velocities and intervals
    np.savetxt('./'+well+'/'+well+'d/Velocities_P_dry_'+suffix+'.out', [Vd,dVd,Vd_down,Vdmc.mean,Vd_up], fmt='%1.2f',delimiter=',',header='Vp (dry), StdVp (dry), V_0025, V_mean, V_0975')
    #Plot the velocities distributions
    plt.figure('Vp (dry)')
    Veld=mc.N(Vd,dVd)
    Veld.plot(label='Vp (dry)',lw=2,color='b')
    Vdmc.plot(hist=True,label='Vp (dry) (MC)',color='g')
    plt.legend()
    plt.xlabel('Velocity (m/s)')
    #Save the figure
    if use_manager:
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
    plt.show()
    plt.savefig('./'+well+'/'+well+'d/Vp_hist_dry_'+suffix+'.pdf',bbox_inches='tight')
    plt.show(block=True)#plt.show(block=True)

    #CALCULATE THE VELOCITIES FROM THE TIMES PICKED FOR THE SATURATED SAMPLE
    Vp,dVp,Vp_down,Vpmc,Vp_up=rp.Velocity_S(L,np.mean(tS),StdL,np.std(tS)) #Velocity S is more general than Velocity_P
    #Save the velocities and intervals
    np.savetxt('./'+well+'/'+well+'s/Velocities_P_sat_'+suffix+'.out', [Vp,dVp,Vp_down,Vpmc.mean,Vp_up], fmt='%1.2f',delimiter=',',header='Vp (sat), StdVp (sat), V_0025, V_mean, V_0975')
    #Plot the velocities distributions
    plt.figure('Vp (sat.)')
    Vels=mc.N(Vp,dVp)
    Vels.plot(label='Vp (sat.)',lw=2,color='b')
    Vpmc.plot(hist=True,label='Vp (sat.) (MC)',color='g')
    plt.legend()
    plt.xlabel('Velocity (m/s)')
    #save the figure
    if use_manager:
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
    plt.show()
    plt.savefig('./'+well+'/'+well+'s/Vp_hist_sat_'+suffix+'.pdf',bbox_inches='tight')
    plt.show(block=True)

    #Plot #5: A Summary plot with the picked intervals
    # definitions for the axes
    left, width = 0.12, 0.60
    bottom, height = 0.08, 0.60
    bottom_h = 0.16 + width
    left_h = left + 0.27
    rect_plot = [left_h, bottom, width, height]
    rect_x = [left_h, bottom_h, width, 0.2]
    rect_y = [left, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure('Summary', figsize=(10, 10))

    #axplot = fig.add_subplot(111)
    axplot = plt.axes(rect_plot)
    axx = plt.axes(rect_x)
    axy = plt.axes(rect_y)

    # Plot the matrix
    axplot.plot(idxqo, idxto, 'k', lw=2)
    axplot.plot(idxqoa, idxtoa, 'm', lw=2)
    axplot.plot(idxqoh, idxtoh, '--', color='m', lw=2)
    # plot range of picked times
    axplot.axvspan(np.min(tS), np.max(tS), alpha=0.5, color='blue')
    # plot mean of the picked times
    axplot.axvline(np.mean(tS), ymin=t_id, ymax=t_fd, linewidth=2, color='b')
    # plot range of picked times
    axplot.axhspan(np.min(tD), np.max(tD), alpha=0.5, color='red')
    # plot mean of the picked times
    axplot.axhline(np.mean(tD), xmin=t_id, xmax=t_fd, linewidth=2, color='r')
    # Give same time scale as template
    axplot.axis([t_id, t_fd, t_id, t_fd])
    axplot.grid(color='0.5')


    #Indices in P correspond to locations where S-times saturated are larger than S-times dry
    P = np.where(idxqo >= idxto) #Find the locations where time of S sat. is larger than time of S dry
    Ph = np.where(idxqoh >= idxtoh)
    #axplot.scatter(idxqo[P],idxto[P],c='r',s=10)
    #axplot.scatter(idxqoh[Ph],idxtoh[Ph],c='r',s=10)
    #Define and plot the 1:1 line of match
    x1 = np.linspace(0, np.max([time_sd, time_ss]), 10)
    y1 = x1
    axplot.plot(x1, y1, 'g') #plot the 1:1 line of match
    axplot.tick_params(axis='both', which='major', labelsize=18)
    # Plot time serie horizontal
    axx.plot(idxq, query, '-', color='b', label='P (sat.)', lw=2)
    axx.plot(idxq, querya, '-', color='m', lw=2)
    # The envelope
    axx.plot(idxq, queryh, '--', color='m', lw=2)
    # plot range of picked times
    axx.axvspan(np.min(tS), np.max(tS), alpha=0.5, color='blue')
    # plot mean of the picked times
    axx.axvline(np.mean(tS), ymin=-10, ymax=10, linewidth=2, color='b')
    axx.grid(color='0.5')
    axx.axis([t_id, t_fd, 1.1 * np.min(query), 1.1 * np.max(queryh)])
    axx.tick_params(axis='both', which='major', labelsize=18)
    # Plot time series vertical
    axy.plot(template, idxt, '-', color='r', label='P (dry)', lw=2)
    axy.plot(templatea, idxt, '-', color='m', lw=2)
    axy.plot(templateh, idxt, '--', color='m', lw=2)
    # plot range of picked times
    axy.axhspan(np.min(tD), np.max(tD), alpha=0.5, color='red')
    # plot mean of the picked times
    axy.axhline(np.mean(tD), xmin=-10, xmax=10, linewidth=2, color='r')
    axy.grid(color='0.5')
    axy.axis([1.1 * np.min(template), 1.1 * np.max(templateh), t_id, t_fd])
    axy.invert_xaxis()
    axy.tick_params(axis='both', which='major', labelsize=18)
    # Limits
    axx.set_xlim(axplot.get_xlim())
    axy.set_ylim(axplot.get_ylim())
    # save the figure
    plt.savefig(
        './{0}/{0}d/Summary_Match_Pwaves_pick_{1}.pdf'.format(well, suffix),
        bbox_inches='tight')
    plt.savefig(
        './{0}/{0}s/Summary_Match_Pwaves_pick_{1}.pdf'.format(well, suffix),
        bbox_inches='tight')
    plt.show(block=True)
    # Show the Fourier Spectra of both signals
    plt.figure('Spectra')
    FS1 = np.abs(np.fft.rfft(s1_final))
    PS1 = np.unwrap(np.angle(np.fft.rfft(s1_final)))
    f1 = np.fft.rfftfreq(len(s1_final), d=0.01)
    FS2 = np.abs(np.fft.rfft(s2_final))
    PS2 = np.unwrap(np.angle(np.fft.rfft(s2_final)))
    f2 = np.fft.rfftfreq(len(s2_final), d=0.01)
    plt.plot(f1, FS1, 'r', lw=2)
    plt.plot(f2, FS2, 'b', lw=2)
    plt.axis([0, 1.2, 0, np.max([FS1, FS2])])
    if labels:
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Power')
    plt.figure('Phase  Spectra')
    plt.plot(f1, PS1, 'r', lw=2)
    plt.plot(f2, PS2, 'b', lw=2)
    plt.axis([0, 1.2, np.min([PS1, PS2]), np.max([PS1, PS2])])
    if labels:
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Unwrapped phase')
    if cross_correlation:
        _cross_correlation(idxt, idxq, query, template, well)


def _cross_correlation(idxt, idxq, query, template, well):
    """Cross-correlation"""
    dtt = np.mean(np.diff(idxt)) #sampling rate of template
    dtq = np.mean(np.diff(idxq)) #sampling rate of query
    if dtt <= dtq:
        dt = dtt
        query = np.interp(idxt, idxq, query)
        time = idxt
    elif dtt > dtq:
        dt = dtq
        template = np.interp(idxq, idxt, template)
        time = idxq
    WS = [5]
    TAU = []
    for j in range(0, len(WS)):
        # window size in microsecs
        window_size_ms = WS[j]
        # window size in number of samples
        window_size_n = np.int(np.round(window_size_ms / dt))
        p = np.arange(np.min(time), np.max(time), 5)
        pos = np.min(np.where(time >= np.max(p)))
        # slide the correlation window
        Corr = []
        Taus = []
        Th = [] #half time of the sliding window
        for i in range(0, len(time)-window_size_n):
            timew = time[i:i+window_size_n]
            C = np.correlate(template[i:i+window_size_n], query[i:i+window_size_n], mode='same') #dry is template, sat is query
            Corr = np.append(Corr, np.max(C))
            lags = np.linspace(-dt*len(C)/2, dt*len(C)/2, len(C))
            pmax = np.where(C == np.max(C))[0]
            Taus = np.append(Taus, lags[pmax])
            Th = np.append(Th, np.min(timew)+(np.max(timew)-np.min(timew))/2)

        print('Mode of Tau for window size '+np.str(window_size_ms)+': ',sp.stats.mode(Taus).mode[0])
        C=np.correlate(template,query,mode='same') #dry is template, sat is query
        pmax=np.where(C==np.max(C))[0]
        lag=np.linspace(-dt*len(C)/2,dt*len(C)/2,len(C))
        tau=lag[pmax] #delta time between waveforms
        TAU=np.append(TAU,tau)

    plt.figure('Windowed waveforms')
    plt.plot(time,template,label='P (dry)',color='r',lw=2)
    plt.plot(time,query,label='P (sat)',color='b',lw=2)
    plt.grid(color='0.5')
    plt.legend()
    plt.xlabel('time ($\mu s$)'),plt.ylabel('A.U.')
    plt.savefig('./'+well+'/'+well+'s/Vp_waveforms_xcorr.pdf',bbox_inches='tight')
    #Cross-correlation
    plt.figure('Cross-correlation')
    plt.plot(lag,C,lw=2)
    plt.axvline(tau,ymin=1.1*np.min(C),ymax=1.1*np.max(C),c='b',lw=2,label='$\Delta t$: '+np.str(np.round(tau[0],decimals=2))+' $\mu s$')
    plt.axvline(sp.stats.mode(Taus).mode[0],ymin=1.1*np.min(C),ymax=1.1*np.max(C),c='r',lw=2,label='$\Delta t_w$: '+np.str(np.round(sp.stats.mode(Taus).mode[0],decimals=2))+' $\mu s$')
    plt.xlabel('lag ($\mu s$)'),plt.ylabel('Correlation Coefficient')
    plt.legend()
    plt.grid(color='0.5')
    #save the figure
    plt.show()
    plt.savefig('./'+well+'/'+well+'s/Vp_sat_Xcorr.pdf',bbox_inches='tight')
    plt.show(block=True)
