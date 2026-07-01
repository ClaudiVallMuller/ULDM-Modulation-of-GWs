# -*- coding: utf-8 -*-
"""
Set of useful functions for python. From plots to mathematical analysis.
Author: Claudi Vall Müller - University of Amsterdam
Last modified: 20.04.2026

To import the functions of this file in another one, use

import sys
sys.path.append(<directory where this .py file is>)
from python_functions_utils import * 

"""
   
# Libraries, constants and functions
import matplotlib.pyplot as plt
import numpy as np
import math

G = 6.6743e-11 # SI
c = 299792458 # SI
h = 6.626070e-34  # SI
hbar = h/(2*np.pi)  # SI
Mpl = 1.2e19 # GeV - Planck mass


Msun = 1.9884e30 # SI
pc = 3.0857e16  # number of meters in a parsec
eV = 1.602177e-19 # number of joules in an electronvolt
yr = 3.1557e7 # number of seconds in a year
pi = np.pi 

def plot_1d_function(xValues, yValues, x2Fun = False, x2FunInv = False, fontsizeV = 15, xLabel = '$x$', x2Label = '', 
                     yLabel = '$y$', titleLabel = "", textLabels = "", nBins = 10, save = False, saveName = "plot_1d.png", savedpi = 500,
                     xScale = 'linear', yScale = 'linear', 
                     colors = ['blue'], plotLabels = [''], xLines = [], LineStyles = ['-'], legendSize = 1, LineWidths = [1], 
                     fillRegions = [], fillColors = [], fillAlphas = [], 
                     fillRegionsX = [], fillColorsX = [], fillAlphasX = [],
                     plot2D = [], plot2Dcolormap = 'viridis', plot2DLabel = '', plot2Dalpha = 1, plot2DScale = 'linear', 
                     plot2DMin = None, plot2DMax = None, plot2DInv = False, contourValues = '', contourColors = 'black',
                    show = True, titleSize = 1, textSize = 1.3, xMin = None, xMax = None, yMin = None, yMax = None, xInv = False, yInv = False,
                    xSizePlot = 1, ySizePlot = 1):
    
    "Plot 1D function (discrete array)."
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib.ticker import MaxNLocator
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator
    from matplotlib.colors import TABLEAU_COLORS, same_color

    try: 
        import scienceplots
        plt.style.use('science')
    except:
#        %pip install SciencePlots
#   
        print("I recommend installing SciencePlots using -- pip install SciencePlots--")
    
    # Parameters of the plot
    if fontsizeV != 0:
        plt.rcParams.update({'font.size': fontsizeV})
        plt.rcParams['figure.figsize'] = [fontsizeV, fontsizeV]
        if len(plot2D) == 0:
            f, ax = plt.subplots(1, 1, figsize = (fontsizeV/2.5 * xSizePlot, fontsizeV/2.5 * ySizePlot))
        else:
            f, ax = plt.subplots(1, 1, figsize = (fontsizeV/2.5 * 1.217, fontsizeV/2.5))
        textSize = fontsizeV * 1.3
    else:
        f, ax = plt.subplots(1, 1)

    
    # Axes and labels
    ax.set_xlabel(xLabel, fontsize = textSize)
    ax.set_ylabel(yLabel, fontsize = textSize)
    ax.set_title(titleLabel, fontsize = textSize)
    ax.set_xscale(xScale)
    ax.set_yscale(yScale)

    if xMin != None or xMax != None:
        ax.set_xlim(xMin, xMax)
    if yMin != None or yMax != None:
        ax.set_ylim(yMin, yMax)

   
    # Axis inversion
    ax.yaxis.set_inverted(yInv)
    ax.xaxis.set_inverted(xInv)
    
    # Additional parameters
    # ax.set_aspect('auto', adjustable = 'box')
    # # Axes ticks
    # plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=nBins))
    # plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=nBins))
    # ax.xaxis.set_minor_locator(AutoMinorLocator())
    # ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.tick_params(axis="both", direction="in", which = "major", width = fontsizeV/10, length = fontsizeV/2, top = True, right = True)
    # ax.tick_params(axis="both", direction="in", which = "minor", width = fontsizeV/10, length = fontsizeV/5, top = True, right = True)

    # Plot data
    try:   # if a set of arrays is given
        LEN = len(yValues[0])
        multiplots = True
        numberPlots = len(yValues)
    except:  # if only one array
        multiplots = False

    if multiplots == True:
        if len(plotLabels) != numberPlots:
            plotLabels = ['' for n in range(numberPlots)]
        if len(colors) != numberPlots:
            colors = [list(TABLEAU_COLORS)[n][4:] for n in range(numberPlots)]
        if len(LineStyles) != numberPlots:
            LineStyles = ['-' for n in range(numberPlots)]
        if len(LineWidths) != numberPlots:
            LineWidths = [1 for n in range(numberPlots)]
        for n in range(numberPlots):
            if LineStyles[n] == "." or LineStyles[n] == "x":
                plot = ax.plot(xValues[n], yValues[n], label = plotLabels[n], c = colors[n], marker = LineStyles[n], markersize = 7 * LineWidths[n], lw = 0)        
            else:    
                plot = ax.plot(xValues[n], yValues[n], label = plotLabels[n], c = colors[n], ls = LineStyles[n], lw = LineWidths[n])        
    
    else:
        if LineStyles == "." or LineStyles == "x":
                plot = ax.plot(xValues, yValues, label = plotLabels, c = colors, marker = LineStyles, markersize = 7 * LineWidths, lw = 0)        
        else:    
                plot = ax.plot(xValues, yValues, label = plotLabels[0], c = colors[0], ls = LineStyles[0], lw = LineWidths[0])        
    

    
    # Plot 2D colormaps
    if contourValues == '':
        contourValues = [None for k in range(len(plot2D))]
    if len(plot2D) != 0:
        if np.size(plot2D[0][0][0]) != 1:    # Second 2D plot is countour
            for k in range(len(plot2D)):
                xGrid = plot2D[k][0]
                yGrid = plot2D[k][1]
                zGrid = plot2D[k][2]
                if plot2DMin[k] == None:
                    plot2DMin[k] = np.min(zGrid)
                if plot2DMin[k] == 0:
                    plot2DMin[k] = 1e-20
                if plot2DMax[k] == None:
                    plot2DMax[k] = np.max(zGrid)
                if plot2DInv[k] == True:
                    plot2Dcolormap[k] = plot2Dcolormap[k] + "_r"
                
                colorPlot = ax.pcolormesh(xGrid, yGrid, zGrid, cmap = plot2Dcolormap[k], alpha = plot2Dalpha[k], norm = plot2DScale[k],
                                         vmax = plot2DMax[k], vmin = plot2DMin[k])

                if contourValues[k] != None:
                    ax.contour(xGrid, yGrid, zGrid, contourValues[k], linewidths = 1, colors = contourColors[k])

                if k==0:
                    cbar = f.colorbar(colorPlot)
                    cbar.ax.set_ylabel(plot2DLabel, fontsize = textSize)
                   
        else:
            xGrid = plot2D[0]
            yGrid = plot2D[1]
            zGrid = plot2D[2]
            if plot2DMin == None:
                plot2DMin = np.min(zGrid)
            if plot2DMin == 0:
                plot2DMin = 1e-20
            if plot2DMax == None:
                plot2DMax = np.max(zGrid)
            if plot2DInv == True:
                plot2Dcolormap = plot2Dcolormap + "_r"
            colorPlot = ax.pcolormesh(xGrid, yGrid, zGrid, cmap = plot2Dcolormap, alpha = plot2Dalpha, norm = plot2DScale,
                                     vmax = plot2DMax, vmin = plot2DMin)
            cbar = f.colorbar(colorPlot)
            cbar.ax.set_ylabel(plot2DLabel, fontsize = textSize)
            ax.contour(xGrid, yGrid, zGrid, contourValues, linewidths = 1, colors = contourColors)



     # Second (upper) x axis
    if (x2Fun and x2FunInv) != False:
        secax = ax.secondary_xaxis('top', functions=(x2Fun, x2FunInv))
        secax.set_xlabel(x2Label, fontsize = textSize)

        # Make sure only the secondary axis has top ticks
        ax.tick_params(axis='x', which='both', top=False)
        secax.tick_params(axis='x', which='both', bottom=False)
            
    
    # Texts in the plot
    if multiplots == True:
        whiteX = xValues[0]
        whiteY = yValues[0]
    else:
        whiteX = xValues
        whiteY = yValues
    for k in range(len(textLabels)):
        plot = ax.plot(whiteX, whiteY, color = 'white', linewidth = 0, label = textLabels[k])    

    # Legend (only for 1D plots for now)
    if len(plot2D) == 0:
        ax.legend(fontsize = fontsizeV * 1.1 * legendSize)
    
    lineTop = ax.get_ylim()[0]
    lineBot = ax.get_ylim()[1]
    # Vertical lines
    for k in range(len(xLines)):
        ax.vlines(xLines[k], lineTop, lineBot, color = 'grey', linestyle = 'dashed')


    # Fill regions Y direction
    if len(fillColors) != len(fillRegions):
        fillColors = [list(TABLEAU_COLORS)[k][4:] for k in range(len(fillRegions))]
    if len(fillAlphas) != len(fillRegions):
        fillAlphas = [0.5 for k in range(len(fillRegions))]
    for k in range(len(fillRegions)):
        ax.fill_between(fillRegions[k][0], fillRegions[k][1], fillRegions[k][2], color = fillColors[k], alpha = fillAlphas[k])

    # Fill regions X direction
    if len(fillColorsX) != len(fillRegionsX):
        fillColorsX = [list(TABLEAU_COLORS)[k][4:] for k in range(len(fillRegionsX))]
    if len(fillAlphasX) != len(fillRegionsX):
        fillAlphasX = [0.5 for k in range(len(fillRegionsX))]
    for k in range(len(fillRegionsX)):
        ax.fill_betweenx(fillRegionsX[k][0], fillRegionsX[k][1], fillRegionsX[k][2], color = fillColorsX[k], alpha = fillAlphasX[k])
    
    # Output
    plt.tight_layout()
    if save == True:
        plt.savefig(saveName, dpi = savedpi)
    if show == False:
        plt.close() 
    else:
        plt.show()


def plot_2d_function(xGrid, yGrid, zGrid, colormap = "hot", colorAlpha = 1, minValue = 0, maxValue = 0, fontsizeV = 15,
                     xScale = 'linear', yScale = 'linear', zScale = 'linear',
                     xLabel = '$x$', yLabel = '$y$', zLabel = '$z$', titleLabel = '', textLabels = '', xLines = [], 
                     xMin = None, xMax = None, yMin = None, yMax = None,
                     x2Fun = False, x2FunInv = False, x2Label = '', xInv = False, yInv = False,
                     nBins = 10, save = False, saveName = "plot_2d.png", savedpi = 200, show = True):
    "Plot 2D function (discrete 2D array) in a colormap."
    from matplotlib.ticker import MaxNLocator, AutoMinorLocator

    try:
        import scienceplots
        plt.style.use('science')
    except:
#        %pip install SciencePlots   
        print("I'd recommend installing SciencePlots using -- pip install SciencePlots--")
    
    # Parameters of the plot
    plt.rcParams.update({'font.size': fontsizeV})
    f, ax = plt.subplots(1, 1, figsize = (fontsizeV/2, fontsizeV/2))
    textSize = fontsizeV * 1.3

    # Axes and labels
    ax.set_xlabel(xLabel, fontsize = textSize)
    ax.set_ylabel(yLabel, fontsize = textSize)
    ax.set_title(titleLabel, fontsize = textSize)

    ax.set_xscale(xScale)
    ax.set_yscale(yScale)
    
    if xMin != None or xMax != None:
        ax.set_xlim(xMin, xMax)
    if yMin != None or yMax != None:
        ax.set_ylim(yMin, yMax)

    # Second (upper) x axis
    if (x2Fun and x2FunInv) != False:
        secax = ax.secondary_xaxis('top', functions=(x2Fun, x2FunInv))
        secax.set_xlabel(x2Label, fontsize = textSize)

        # Make sure only the secondary axis has top ticks
        ax.tick_params(axis='x', which='both', top=False)
        secax.tick_params(axis='x', which='both', bottom=False)
    
    
    # Axis inversion
    ax.yaxis.set_inverted(yInv)
    ax.xaxis.set_inverted(xInv)
    
    ax.set_aspect('auto')

    # Axes ticks
    # plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=nBins))
    # plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=nBins))    
    # ax.xaxis.set_minor_locator(AutoMinorLocator())
    # ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.tick_params(axis="both", direction="in", which = "major", width = fontsizeV/10, length = fontsizeV/2, top = True, right = True)
    # ax.tick_params(axis="both", direction="in", which = "minor", width = fontsizeV/10, length = fontsizeV/5, top = True, right = True)


    # # Vertical lines
    # for k in range(len(xLines)):
    #     ax.vlines(xLines[k], ax.get_ylim()[0], ax.get_ylim()[1], color = 'gray', linestyle = 'dashed')

    
    # Plot data
    if minValue == 0:
        minValue = np.min(zGrid)
    if maxValue == 0:
        maxValue = np.max(zGrid)
        
    plot = ax.pcolormesh(xGrid, yGrid, zGrid, cmap = colormap, vmax = maxValue, vmin = minValue, alpha = colorAlpha, norm = zScale)
    cbar = f.colorbar(plot)
    cbar.ax.set_ylabel(zLabel, fontsize = textSize)


    # Texts in the plot
    whiteX = xGrid[0]
    whiteY = yGrid[0]
    for k in range(len(textLabels)):
        plot = ax.plot(whiteX, whiteY, color = 'white', linewidth = 0, label = textLabels[k])    

    
    # Legend
    ax.legend(fontsize = fontsizeV * 1.1)
    
    # Output
    plt.tight_layout()
    if save == True:
        plt.savefig(saveName, dpi = savedpi)

    if show == False:
        plt.close()
    else:
        plt.show()


# Define function for string formatting of scientific notation [https://stackoverflow.com/questions/18311909/how-do-i-annotate-with-power-of-ten-formatting]
def sci_notation(num, decimal_digits=0, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    
    if num == 0:
        return r"$0$"
    else:
        if exponent is None:
            exponent = int(np.floor(np.log10(abs(num))))
        coeff = round(num / float(10**exponent), decimal_digits)
        if precision is None:
            precision = decimal_digits
    
        if exponent == 0:
            return r"${0:.{1}f}$".format(coeff, precision)
    
        if coeff == 1:
            if exponent == 1:
                return r"$10$"
            else:
                return r"$10^{{{1:d}}}$".format(coeff, exponent, precision)
    
        else:
            if exponent == 1:
                return r"${0:.{2}f} \times 10$".format(coeff, exponent, precision)
            else:
                return r"${0:.{2}f} \times 10^{{{1:d}}}$".format(coeff, exponent, precision)
    
    
    
    
### Compute slope from an array
def slope_array(array, delta):
    '''Returs array with the value of the slope at each point of an input array (array). (delta) is the 
    spacing, assumed constant, of the independent variable between two points of (array).'''
    
    arrayMinus = np.insert(array, 0, array[0])
    arrayMinus = np.delete(arrayMinus, len(arrayMinus) - 1)
    arrayPlus = np.delete(array, 0)
    arrayPlus = np.insert(arrayPlus, len(array) - 1, array[-1])
    
    # Inner
    arraySlope = (arrayPlus - arrayMinus)/(2*delta)
    
    # Extrema
    arraySlope[0] = (arrayPlus[0] - array[0])/delta
    arraySlope[-1] = (array[-1] - arrayMinus[-1])/delta
    
    return arraySlope

### Detect indices of positions where the two consecutive elements in an array change sign 
def sign_change(array):
    '''Returns indices array of the positions where the values of an array change their sign.'''

#     arrayMinus = np.insert(array, 0, array[0])
#     arrayMinus = np.delete(arrayMinus, len(arrayMinus) - 1)
    arrayPlus = np.delete(array, 0)
    arrayPlus = np.insert(arrayPlus, len(array) - 1, array[-1])

    # Inner
    arrayProd = arrayPlus * array
    arrayProd[np.where(arrayProd == 0)] = 1e-5
    arraySign = arrayProd/np.abs(arrayProd)
    arrayIndices = np.where(arraySign == -1)
    return arrayIndices[0]




def FT_to_time_waveform(freqs, Ampl, Phase, t_coal=0.0):
    """Computes and returns the inverse Fourier transform, to go from the
    waveform in Fourier space to the waveform in time domain.

    Arguments:
    - freqs:  array with (equally spaced) frequency values (Hz)
    - Ampl:   amplitude (modulus) of the Fourier transform at freqs
    - Phase:  phase (argument) of the Fourier transform at freqs (radians)
    - t_coal: coalescence time shift in seconds (default: 0.0)

    Returns (time_values, ht_values):
    - time_values: array of time values (s)
    - ht_values:   waveform in time domain
    """

    # Base Fourier-domain waveform (for t_coal = 0)
    h_fd = Ampl * np.exp(1j * Phase)

    # Apply time shift t_coal in frequency domain:
    # h(t - t_coal)  <->  e^{-2π i f t_coal} h̃(f)
    h_fd *= np.exp(-2j * np.pi * freqs * t_coal)

    df   = freqs[1] - freqs[0]
    fmin = np.min(freqs)
    fcut = np.max(freqs)

    # Choose sampling and length
    f_Nyquist = max(2 * fcut, 4096.0)       # Hz
    N         = int(round(f_Nyquist / df))  # number of time samples
    dt        = 1.0 / (N * df)
    T_total   = N * dt

    # Full FFT array H[k], k = 0 … N-1
    H_full = np.zeros(N, dtype=complex)

    k_min = int(round(fmin / df))
    k_max = k_min + len(freqs)
    k_max = min(k_max, N // 2 + 1)
    n_use = k_max - k_min

    H_full[k_min:k_max] = h_fd[:n_use] / dt

    # Hermitian symmetry
    H_full[N - k_max + 1 : N - k_min + 1] = np.conj(h_fd[:n_use][::-1]) / dt

    # Inverse FFT -> time series
    h_td_complex = np.fft.ifft(H_full)
    ht_values    = np.real(h_td_complex)

    # Time array (still centered so that coalescence ~ 0 for t_coal = 0)
    time_values = np.arange(N) * dt - T_total + t_coal

    return time_values, ht_values



def FT_to_time_waveform_no_tshift(freqs, Ampl, Phase):
    """
    Computes and returns the inverse Fourier transform, to go from the
    waveform in Fourier space to the waveform in time domain.

    Arguments:
    - freqs:  array with (equally spaced) frequency values (Hz)
    - Ampl:   amplitude (modulus) of the Fourier transform at freqs
    - Phase:  phase (argument) of the Fourier transform at freqs (radians)

    Returns (time_values, ht_values):
    - time_values: array of time values (s), starting at t = 0
    - ht_values:   waveform in time domain

    NOTE:
    - No explicit time shift is imposed.
    - The effective coalescence time is entirely determined by Phase(freqs).
    """

    # Fourier-domain waveform (no imposed time shift)
    h_fd = Ampl * np.exp(1j * Phase)

    df   = freqs[1] - freqs[0]
    fmin = np.min(freqs)
    fcut = np.max(freqs)

    # Choose sampling and length (same logic as your original function)
    f_Nyquist = max(2 * fcut, 4096.0)       # Hz (target 2 * Nyquist)
    N         = int(round(f_Nyquist / df))  # number of time samples
    dt        = 1.0 / (N * df)
    T_total   = N * dt

    # Full FFT array H[k], k = 0 … N-1
    H_full = np.zeros(N, dtype=complex)

    k_min = int(round(fmin / df))
    k_max = k_min + len(freqs)
    k_max = min(k_max, N // 2 + 1)
    n_use = k_max - k_min

    # Positive frequencies (including DC / low frequency part)
    H_full[k_min:k_max] = h_fd[:n_use] / dt

    # Hermitian symmetry for real time series
    H_full[N - k_max + 1 : N - k_min + 1] = np.conj(h_fd[:n_use][::-1]) / dt

    # Inverse FFT -> time series
    h_td_complex = np.fft.ifft(H_full)
    ht_values    = np.real(h_td_complex)

    # Time array: no explicit centering or coalescence alignment
    time_values = np.arange(N) * dt  # t = 0, dt, ..., (N-1)*dt

    return time_values, ht_values


