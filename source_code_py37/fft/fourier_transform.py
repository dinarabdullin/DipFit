'''
Calculate FFT and inverse FFT
'''

import os
import sys
import wx
import math
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['axes.facecolor']= 'white'
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['lines.linewidth'] = 2
rcParams['xtick.major.size'] = 8
rcParams['xtick.major.width'] = 1.5
rcParams['ytick.major.size'] = 8
rcParams['ytick.major.width'] = 1.5
rcParams['font.size'] = 14
linestyles = ['k-', 'r-', 'b-', 'm-', 'c-']


def get_path(message):
    app = wx.App(None) 
    dialog = wx.FileDialog(None, message, wildcard='*.dat', style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
    if dialog.ShowModal() == wx.ID_OK:
        path = dialog.GetPath()
    return path


def load_data(filename, column_number):
    x = []
    y = []
    file = open(filename, 'r')
    for line in file:
        str = line.split()
        x.append( float(str[0]) )
        y.append( float(str[column_number]) )
    file.close()
    xv = np.array(x)
    yv = np.array(y)
    return [xv, yv]


def crop_signal(x, y, xrange):
    yCropped = []
    for i in range(x.size):
        if (x[i] >= xrange[0]) and (x[i] <= xrange[1]):
            yCropped.append(y[i])
    yvCropped = np.array(yCropped)
    return yvCropped
    
    
def plot_timetrace(X, Y, legend_text = []):
    xmin = np.amin([np.amin(x) for x in X])
    xmax = np.amax([np.amax(x) for x in X])
    ymin = np.amin([np.amin(y) for y in Y])
    ymax = np.amax([np.amax(y) for y in Y])
    fig = plt.figure(facecolor='w', edgecolor='w')
    axes = fig.gca()
    for i in range(len(X)):
        axes.plot(X[i], Y[i], linestyles[i])
    axes.set_xlim(xmin, xmax)
    axes.set_ylim(ymin-0.1, ymax+0.1)
    axes.set_xlabel(r'$\mathit{t}$ ($\mathit{\mu s}$)')
    axes.set_ylabel(r'$Signal$')
    if not (legend_text == []):
        axes.legend(legend_text)
    plt.tight_layout()
    plt.draw()
    plt.show(block=False)


def plot_spectrum(X, Y, legend_text = []):
    xmin = np.amin([np.amin(x) for x in X])
    xmax = np.amax([np.amax(x) for x in X])
    ymin = np.amin([np.amin(y) for y in Y])
    ymax = np.amax([np.amax(y) for y in Y])
    fig = plt.figure(facecolor='w', edgecolor='w')
    axes = fig.gca()
    for i in range(len(X)):
        axes.plot(X[i], Y[i], linestyles[i])
    axes.set_xlim(xmin, xmax)
    axes.set_ylim(ymin-0.1, ymax+0.1)
    axes.set_xlabel(r'$\mathit{f}$ ($\mathit{MHz}$)')
    axes.set_ylabel(r'$Amplitude$')
    if not (legend_text == []):
        axes.legend(legend_text)
    plt.tight_layout()
    plt.draw()
    plt.show(block=False)
    return [fig, axes]
    

def keep_figures_live():
	plt.show()


def save_spectrum(filename, x, y):
    file = open(filename, 'w')
    if isinstance(y[0], complex):
        for i in range(x.size):    
            file.write('{0:<12.6f}'.format(x[i]))
            file.write('{0:<12.6f}'.format(np.real(y[i])))
            file.write('{0:<12.6f}'.format(np.imag(y[i])))
            file.write('\n')  
    else:
        for i in range(x.size):    
            file.write('{0:<12.6f}'.format(x[i]))
            file.write('{0:<12.6f}'.format(y[i]))
            file.write('\n')  
    file.close()


def offset_correction(sig, bckg_start):
    bckg_level = np.mean(sig[-bckg_start:])
    bckg = bckg_level * np.ones(sig.size)
    sigShift = sig - bckg
    return [sigShift, bckg]


def appodization(sig, mode='None', parameters={}):
    sigApp = []
    if (mode == 'None'):
        window = np.ones(sig.size)
        sigApp = sig
    if (mode == 'Hamming'):
        hamming = np.hamming(2*sig.size-1)
        window = hamming[-sig.size:]
        sigApp = window * sig
    return [sigApp, window]


def zerofilling(t, sig, length = 0):
    tZF = []
    sigZF = []
    if (length == 0):
        tZF = t
        sigZF = sig
    else:
        tmax = t[0] + (t[1]-t[0]) * float((length+1)*t.size-1)
        tZF = np.linspace(t[0], tmax, (length+1)*t.size)
        sigZF = np.append(sig, np.zeros(length*sig.size))
    return [sigZF, tZF]


def scale_first_point(sig):
    sigScale = sig
    sigScale[0] = 0.5 * sigScale[0]
    return sigScale
    

def FFT(t, sigRe, sigIm, parameters, fRef = [], spcRef = []): 
    if (parameters['complex']):
        sig = sigRe + sigIm * 1j
    else:
        sig = sigRe
    
    # Shift the spectrum such that the last time points have a zero amplitude
    [sigShift, bckg] = offset_correction(sig, parameters['bckg_start'])
    plot_timetrace([t, t, t, t], [sigRe, sigIm, bckg, sigShift], ['Re(sig)', 'Im(sig)', 'bckg', 'Re(sig)-bckg'])
    
    # Apply the appodization
    sigApp, window = appodization(sigShift, mode=parameters['appod_func'], parameters={})
    plot_timetrace([t, t, t], [sigShift, window, sigApp], ['sig', 'window', 'window*sig']) 
    
    # Apply zero-padding
    [sigZF, tZF] = zerofilling(t, sigApp, parameters['zerofil_length'])
    plot_timetrace([t, tZF], [sigApp, sigZF], ['sig', 'sig+zeros'])

    # Scale first point by 1/2
    if parameters['scale_first_point']:
        sigScale = scale_first_point(sigZF)
    else:
        sigScale = sigZF
    
    # Apply FFT
    spc = np.fft.fft(sigScale)
    dt = t[1] - t[0]
    f = np.fft.fftfreq(spc.size, dt)
    
    # Re-organize the data
    fShift = np.fft.fftshift(f)
    spcShift = np.fft.fftshift(spc)
    
    # Calculate the Re, Im, and Abs of the spectrum
    spcRe = np.real(spcShift)
    spcIm = np.imag(spcShift)
    spcAbs = np.abs(spcShift)
    plot_spectrum([fShift, fRef], [spcRe/np.amax(spcRe), spcRef/np.amax(spcRef)], ['Re(spc)', 'DeerAnalysis'])
    
    # Set the data to be returned
    if (parameters['output_data'] == 'orig'):
        ans = spc
    elif (parameters['output_data'] == 'real'):
        ans = spcRe
    elif (parameters['output_data'] == 'imag'):
        ans = spcIm
    elif (parameters['output_data'] == 'abs'):
        ans = spcAbs
    
    return [fShift, ans, f, spc]
    

def inverse_FFT(f, spc, tRef = [], sigRef = []): 
    # Apply inverse FFT
    sig = np.fft.ifft(spc)
    fNuq = np.amax(abs(f))
    dt = 1 / (2 * fNuq)
    Nt = f.size
    t = np.linspace(0, dt*(Nt-1), Nt)
    
    # Restore the amplitude of first data point
    sig[0] = 2.0 * sig[0]
    
    # Shift the time trace
    sig = sig + (1.0-sig[0])*np.ones(sig.size)
    
    # Calculate the Re of the time trace
    sigRe = np.real(sig)
    sigIm = np.real(sig)
    plot_timetrace([t, tRef], [sigRe, sigRef], ['Re(signal)', 'reference'])

    return [t, sigRe, sigIm]


if __name__ == '__main__':
    # Load the time trace
    filename = get_path("Open the file with a background-corrected time trace...")
    #filename = "D:/Project/Measurements/Data/Brehm/PulseEPR/2018_08_13_Iron(III)_HS_LS/Analysis/ChK52_RIDME_11984G_2/ChK52_10000eq_HIm_200uM_THFd8_RIDME_11984G_10K_fit.dat"
    [t, sigRe] = load_data(filename, 1)
    sigIm = np.zeros(sigRe.size)

    # Load the spectrum
    filename = get_path("Open the file with a spectrum...")
    #filename = "D:/Project/Measurements/Data/Brehm/PulseEPR/2018_08_13_Iron(III)_HS_LS/Analysis/ChK52_RIDME_11984G_2/ChK52_10000eq_HIm_200uM_THFd8_RIDME_11984G_10K_spc.dat"
    [fRef, spcRef] = load_data(filename, 1)
    
    # Calculate the FFT of the time trace
    # parameters| complex           | True / False
    #             bckg_start        | Any integer number below the number of points
    #             appod_func        | None / Hamming
    #             zerofil_length    | Any integer number
    #             scale_first_point | True / False
    #             output_data       | orig / real / imag / abs
    parameters = {}
    parameters['complex'] = False 
    parameters['bckg_start'] = -50 
    parameters['appod_func'] = 'Hamming'
    parameters['zerofil_length'] = 3 
    parameters['scale_first_point'] = True
    parameters['output_data'] = 'real' 
    [f, spc, fRaw, spcRaw] = FFT(t, sigRe, sigIm, parameters, fRef, spcRef)
    
    # Save the FFT
    filename = os.path.dirname(os.path.abspath(filename)) + "\spc.dat"
    save_spectrum(filename, f, spc)
    
    ## The inverse FFT of the spectrum
    #[tInv, sigReInv, sigImInv] = inverse_FFT(fRaw, spcRaw, t, sigRe)

    keep_figures_live()