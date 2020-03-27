'''
fourier_transform.py

Calculates FFT and inverse FFT

Requirements: Python3, numpy, scipy, matplotlib
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


default_parameters = {
    'complex': True, 
    'background_start': -50,
    'appodization': True,
    'zerofilling': 1,
    'scale_first_point': True,   
    'output': 'real'
}


def get_path(message):
    app = wx.App(None) 
    dialog = wx.FileDialog(None, message, wildcard='*.dat', style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
    if dialog.ShowModal() == wx.ID_OK:
        path = dialog.GetPath()
    else:
        path=""
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


def read_fft_parameters():
    parameters = {}
    print('\nEnter FFT settings:\n')
    var = input("Complex input data: True or False (default: %s)\n" % default_parameters['complex'])
    if (var == ""):
        parameters['complex'] = default_parameters['complex']
    else:
        if (var == 'True') or (var == 'true'):
            parameters['complex'] = True
        elif (var == 'False') or (var == 'false'):
            parameters['complex'] = False
        else:
            raise ValueError('Illelgible value!')
            sys.exit(1)
    var = input("Background start counted from the last data point (default: %d)\n" % default_parameters['background_start'])
    if (var == ""):
        parameters['background_start'] = default_parameters['background_start']
    else:
        val = [int(i) for i in var.split(' ')]
        if len(val) == 1:
            if val[0] > 0:
                val[0] = -val[0]
            parameters['background_start'] = val[0]
        else:
            raise ValueError('Illelgible value!')
            sys.exit(1)
    var = input("Appodization: True or False (default: %s)\n" % default_parameters['appodization'])
    if (var == ""):
        parameters['appodization'] = default_parameters['appodization']
    else:
        if (var == 'True') or (var == 'true'):
            parameters['appodization'] = True
        elif (var == 'False') or (var == 'false'):
            parameters['appodization'] = False
        else:
            raise ValueError('Illelgible value!')
            sys.exit(1)
    var = input("Zero filling: [entered value] x [number of data points] (default: %d)\n" % default_parameters['zerofilling'])
    if (var == ""):
        parameters['zerofilling'] = default_parameters['zerofilling']
    else:
        val = [int(i) for i in var.split(' ')]
        if len(val) == 1:
            if val[0] >= 0:
                parameters['zerofilling'] = val[0]
            else:
                raise ValueError('Illelgible value!')
                sys.exit(1)
        else:
            raise ValueError('Illelgible value!')
            sys.exit(1)
    var = input("Scale firt data point: True or False (default: %s)\n" % default_parameters['scale_first_point'])
    if (var == ""):
        parameters['scale_first_point'] = default_parameters['scale_first_point']
    else:
        if (var == 'True') or (var == 'true'):
            parameters['scale_first_point'] = True
        elif (var == 'False') or (var == 'false'):
            parameters['scale_first_point'] = False
        else:
            raise ValueError('Illelgible value!')
            sys.exit(1)
    var = input("Output data format: original, real, imaginary, absolute (default: %s)\n" % default_parameters['output'])
    if (var == ""):
        parameters['output'] = default_parameters['output']
    else:
        if (var == 'Original') or (var == 'original'):
            parameters['output'] = 'original'
        elif (var == 'Real') or (var == 'real'):
            parameters['output'] = 'real'
        elif (var == 'imaginary') or (var == 'imaginary'):
            parameters['output'] = 'imaginary'
        elif (var == 'Absolute') or (var == 'absolute'):
            parameters['output'] = 'absolute' 
        else:
            raise ValueError('Illelgible value!')
            sys.exit(1)
    print(parameters)
    return parameters


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


def appodization(sig, activate):
    sigApp = []
    if (activate):
        hamming = np.hamming(2*sig.size-1)
        window = hamming[-sig.size:]
        sigApp = window * sig
    else:
        window = np.ones(sig.size)
        sigApp = sig
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
    # Specify the signal
    if (parameters['complex']):
        sig = sigRe + sigIm * 1j
    else:
        sig = sigRe
    # Shift the spectrum such that the last time points have a zero amplitude
    [sigShift, bckg] = offset_correction(sig, parameters['background_start'])
    #plot_timetrace([t, t, t], [sigRe, bckg, sigShift], ['real(signal)', 'bckg', 'real(signal) - bckg'])
    # Apply the appodization
    sigApp, window = appodization(sigShift, parameters['appodization'])
    #plot_timetrace([t, t, t], [sigShift, window, sigApp], ['signal', 'window', 'signal * window']) 
    # Apply zero-padding
    [sigZF, tZF] = zerofilling(t, sigApp, parameters['zerofilling'])
    #plot_timetrace([t, tZF], [sigApp, sigZF], ['signal', 'signal + zeros'])
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
    # Set the data to be returned
    if (parameters['output'] == 'original'):
        ans = spc
    elif (parameters['output'] == 'real'):
        ans = spcRe
    elif (parameters['output'] == 'imaginary'):
        ans = spcIm
    elif (parameters['output'] == 'absolute'):
        ans = spcAbs
    # Plot the result
    if (fRef == []):
        plot_spectrum([fShift], [ans/np.amax(ans)])
    else:    
        plot_spectrum([fShift, fRef], [ans/np.amax(ans), spcRef/np.amax(spcRef)], ['calculated', 'reference'])
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
    # Plot the result
    if (tRef == []):
        plot_timetrace([t], [sigRe])
    else: 
        plot_timetrace([t, tRef], [sigRe, sigRef], ['calculated', 'reference'])
    return [t, sigRe, sigIm]


if __name__ == '__main__':
    # Load the time trace
    filename_input_timetrace = get_path("Load a time trace")
    if (filename_input_timetrace == ""):
        raise ValueError('Could not load the time trace!')
        sys.exit(1)
    else:
        print('\nThe time trace is loaded from \n%s' % filename_input_timetrace)
        [t, sigRe] = load_data(filename_input_timetrace, 1)
        sigIm = np.zeros(sigRe.size)
    # Load the spectrum
    filename_input_spectrum = get_path("Load a spectrum (optional)")
    if (filename_input_spectrum == ""):
        print('\nNo spectrum was specified!')
        fRef = []
        spcRef = []
    else:
        print('\nThe spectrum is loaded from \n%s' % filename_input_spectrum)
        [fRef, spcRef] = load_data(filename_input_spectrum, 1) 
    # Enter the FFT parameters
    parameters = read_fft_parameters() 
    # Calculate the FFT of the dipolar time trace
    [f, spc, fRaw, spcRaw] = FFT(t, sigRe, sigIm, parameters, fRef, spcRef)
    filename_output_spectrum = os.path.splitext(filename_input_timetrace)[0] + "_fft.dat"
    save_spectrum(filename_output_spectrum, f, spc)
    ## The inverse FFT of the spectrum
    #[tInv, sigReInv, sigImInv] = inverse_FFT(fRaw, spcRaw, t, sigRe)
    keep_figures_live()