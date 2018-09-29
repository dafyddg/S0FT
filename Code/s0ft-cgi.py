#!/usr/pkgsrc/20140707/bin/python2.7
# -*- coding: utf-8 -*-
from __future__ import division
print "Content-type: text/html\n\n"

##!/usr/pkgsrc/20140707/bin/python2.7
##!/usr/bin/python
# s0ft-cgi.py
# D. Gibbon
# 2018-08-08
# S0FT - Simple F0 Tracker

#==============================================================
# Hard-wired stuff

opsys = "linux"
# opsys = "solaris"

def inithtml():
	print "<html><head><title>S0FT - Simple F0 Tracker - D. Gibbon</title></head>"
	print "<body>"
	return

cgifields = ['metadata','filebase','fft','f0x','peaks','polymodel','figwidth','figheight','f0min','f0max','f0window','f0stepfactor','f0clipmin','f0clipmax','f0median','downsample','cutofflo','orderlo','cutoffhi','orderhi','centreclip','voicewin','voicediff','voicetype']

#==============================================================

import os, sys, re
import cgi, cgitb; cgitb.enable()
from cgi import escape

import numpy as np
from numpy.fft import rfft, irfft       # irfft not used
from numpy import argmax, sqrt, mean, absolute, linspace, log10, logical_and, average, diff, correlate

import matplotlib
matplotlib.use('Agg')   # WWW and saved graphics files
from matplotlib.mlab import find
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

if opsys == "linux":
	import mpld3	# Stream to web

from scipy.signal import medfilt
from scipy.signal import blackmanharris, fftconvolve
import scipy.io.wavfile as wave
from scipy.signal import butter, lfilter, freqz, filtfilt, hann, hamming
from scipy.interpolate import interp1d

#==============================================================

def polyregline(x,y,d):
        x = range(len(y))
        fit, res, _, _, _ = np.polyfit(x, y, d, full=True)
        yfit = np.polyval(fit,x)
        return yfit
# fit with np.polyfit; linear (degree 1) regression only

def polyregline2(x,y,d):
        m, b, _, _, _ = stats.linregress(x, y)
        poly = [ m * xx + b for xx in x]
        return poly

#==============================================================
# Butterworth filters

def butter_lowpass(cutoff, fs, order):
	nyq = 0.5 * fs
	normal_cutoff = cutoff / nyq
	b, a = butter(order, normal_cutoff, btype='low', analog=False)
	return b, a

def butter_lowpass_filter(data, cutoff, fs, order):
	b, a = butter_lowpass(cutoff, fs, order=order)
	y = lfilter(b, a, data)
	return y

def butter_highpass(cutoff, fs, order):
	nyq = 0.5 * fs
	normal_cutoff = cutoff / nyq
	b, a = butter(order, normal_cutoff, btype='high', analog=False)
	return b, a

def butter_highpass_filter(data, cutoff, fs, order):
	b, a = butter_highpass(cutoff, fs, order=order)
	y = lfilter(b, a, data)
	return y

def butterworth(signal,fs,cutoffhi,orderhi,cutofflo,orderlo):
	buttery = signal
	buttery = butter_highpass_filter(buttery, cutoffhi, fs, orderhi)
	buttery = butter_lowpass_filter(buttery, cutofflo, fs, orderlo)
	return buttery

#==============================================================
# F0 detection

def parabolic(f, x):
    xv = 1.0/2.0 * (f[x-1] - f[x+1]) / float(f[x-1] - 2 * f[x] + f[x+1]) + x
    yv = f[x] - 1.0/4.0 * (f[x-1] - f[x+1]) * (xv - x)
    return (xv, yv)

#=====================================

def freq_from_fft(sig, fs):
    windowed = sig * blackmanharris(len(sig))
    f = rfft(windowed)
    i = argmax(abs(f)) # Just use this for less-accurate, naive version
    true_i = parabolic(abs(f), i)[0]
    return fs * true_i / len(windowed)

#=====================================

def freq_from_crossings(sig, fs):
	indices = find((sig[1:] >= 0) & (sig[:-1] < 0))
	crossings = [i - sig[i] / (sig[i+1] - sig[i]) for i in indices]
	diffc = diff(crossings)

	# BUG: np.average does not like empty lists
	if len(diffc) > 0:
		return 1.0 * fs / average(diff(crossings))
	else:
		return 0

#=====================================
# Voice detection

def voicedetector(f0vector,voicewin,voicediff):
    print voicewin,voicediff
    f0new = []
    for i in range(len(f0vector[:-voicewin])):
        dispwin = f0xvector[i:i+voicewin]
        if np.mean(np.abs(np.diff(dispwin))) < voicediff:
            f0new += [f0xvector[i]]
        else:
            f0new += [0]
    return f0new

#=====================================

def freq_from_peaks(sig,fs):
	sig = diff(sig**2)
	return freq_from_crossings(sig,fs)/2.0	

def f0x(y,fs,f0winsamples,f0step,f0median,f0min,f0max,centreclip,f0clipmin,f0clipmax):
        yabs = np.abs(y)
        ymax = np.max(yabs)
        ymin = np.min(yabs)
        ydiff = ymax-ymin

        clip = ymin + int(round(ydiff * centreclip))
#       perc = np.percentile(y,clipfactor)
        y = np.array([ yy if abs(yy)>clip else 0 for yy in y ])

        hannwin = hann(f0winsamples)    # Zero ends
        hammingwin = hamming(f0winsamples,sym=True)    # Non-zero ends
        hammingzwin = hamming(f0winsamples,sym=False)   # Zero ends

	fft = cgiparams['fft']
	f0x = cgiparams['f0x']
	peaks = cgiparams['peaks']

        f0xvector = []
        for i in range(0,len(y[:])-f0winsamples,f0step):
                sigwin = y[i:i+f0winsamples]
##               sigwin = hannwin * sigwin
##               sigwin = hammingwin * sigwin
##               sigwin = hammingzwin * sigwin
#                 f1 = freq_from_fft(sigwin,fs)
                f2 = freq_from_crossings(sigwin,fs)
                f3 = freq_from_peaks(sigwin,fs)

#                if fft == "on" and f0x == "on" and peaks == "on":
#                        f = np.array((f1 * f2 * f3)**0.333)
#                elif fft == "on" and peaks == "on":
#                        f = np.array((f1 * f3)**0.5)
#                elif fft == "on" and f0x == "on":
#                        f = np.array((f1 * f2)**0.5)
                if f0x == "on" and peaks == "on":
                        f = np.array((f2 * f3)**0.5)
#                elif fft == "on":
#                	f = f1
		elif f0x == "on":
                        f = f2
		elif peaks == "on":
			f = f3
                try:
                        fint = int(round(f))
                except:
#                        print "Indeterminate (nan exception)"
                        fint = np.median(sigwin)
#                       fint =f0min
#                       fint = 0
                f0xvector += [fint]
# F0 median filter
        f0xvector = medfilt(f0xvector,f0median)
# F0 low and high clipping
        f0xvector = np.array([ yy if yy>f0clipmin else 0 for yy in f0xvector])
        f0xvector = np.array([ yy if yy<f0clipmax else 0 for yy in f0xvector])

        return f0xvector

#==============================================================
# RAPT

if opsys == "linux":

    def f0rapt(filebase,f0min,f0max,framerate,freqweight):

        wavfilename = "Data/"+filebase + ".wav"
        paramfilename = "params"

	if opsys == "linux":
            f0filename = "/var/www/Webtmp/"+filebase + ".f0"
            f0csvfilename = "/var/www/Webtmp/"+filebase + ".csv"
            f0logfilename = "/var/www/Webtmp/"+filebase + ".f0log"
        else:
            f0filename = "../../docs/Webtmp/"+filebase + ".f0"
            f0csvfilename = "../../docs/Webtmp/"+filebase + ".csv"
            f0logfilename = "../../docs/Webtmp/"+filebase + ".f0log"

        paramstring = "float min_f0 = "+str(f0min)+";\nfloat max_f0 = " + str(f0max)+";\nfloat frame_step = " + str(framerate) + ";\nfloat freq_weight = " + str(freqweight) + ";\n"
        handle = open(paramfilename, "w")
        handle.write(paramstring)
        handle.close()
        getf0filename = "/opt/esps/bin/get_f0"
        getpplainfilename = "/opt/esps/bin/pplain"
        getf0command = getf0filename+" -P "+paramfilename+" "+wavfilename + " " + f0filename
        getcsvcommand = getpplainfilename+" "+f0filename + " > " + f0csvfilename
        os.system(getf0command + " 2> " + f0logfilename)
        os.system(getcsvcommand + " 2>> " + f0logfilename)
        cleanf0filescommand = "rm -f *.f0 *.csv *log"
#       os.system(cleanf0filescommand)
        handle = open(f0csvfilename,"r")
        csvlines = handle.readlines()
        handle.close()
        csvtable = [ line.split(" ")[:-1] for line in csvlines if line != '' ]
        if len(csvtable) < 1:
                print "No f0 output detected."
                exit()
        csvarray = np.array([ map(float,x) for x in csvtable])
        f0list = csvarray[:,0]
        voicelist = csvarray[:,1]
        return f0list,voicelist

#==============================================================
# Fetch CGI parameters

def cgitransferlines(cgifields):
    fs = cgi.FieldStorage()
    fieldvalues = []
    for field in cgifields:
        if fs.has_key(field):
            fieldname = fs[field].value
        else:
            fieldname = '9876'
        fieldvalues = fieldvalues + [(field,fieldname)]
    return dict(fieldvalues)

#==============================================================
#==============================================================
# Fetch CGI parameters as a dictionary

cgiparams = cgitransferlines(cgifields)

paramlist = map(str,[ cgiparams[x] for x in cgifields[1:] ])
paramzipped = zip(cgifields[1:],paramlist)

paramjoined = [ ':'.join(x) for x in paramzipped ]
paramsep01 = 8
paramsep02 = 13
paramsepall = '; '
paramstringio = paramsepall.join(paramjoined[:paramsep01])
paramstringf0proc = paramsepall.join(paramjoined[paramsep01:paramsep02])
paramstringwaveproc = paramsepall.join(paramjoined[paramsep02:])

#==============================================================
#==============================================================
# Fetch waveform and sample rate

try:
    filebase = cgiparams['filebase']
    wavfilename = "Data/" + filebase + ".wav"
    fs,signal = wave.read(wavfilename)

    nykvistfrequency = 2000
    downsamplemax = int(round(1.0*fs/nykvistfrequency))

    try:
        downsample = int(cgiparams['downsample'])
    except:
        print "Warning: downsampling factor is "+str(downsample)+", must be a positive integer.<br>"

    if downsample < 1:
        print "Warning: downsampling factor is "+str(downsample)+", must be > 1<br>."

    if downsample > downsamplemax:
        print "Warning: Downsampling factor ",downsample,"infringes a relatively safe Nykvist condition of",nykvistfrequency,"Hz for sampling speech F0.<br>Infringement of the Nykvist condition (F<sub>sample</sub> > 2&#8901;F<sub>max</sub>) will distort the F0 track and narrow the spectrogram range, and may trigger a processing error.<br>"

#==============================================================
# Load signal

    signal = signal[::downsample] 
    fs = int(round(1.0*fs/downsample))

    signallength = len(signal)
    signalduration = float(signallength) / float(fs)

#==============================================================
# Decode CGI parameters

    voicetype = cgiparams['voicetype']

    if voicetype == 'low':
        f0min = 60
        f0max = 250
        cutoffhi = 95
        orderhi = 3
        cutofflo = 135
        orderlo = 5
        f0clipmin = 60
        f0clipmax = 250
        centreclip = 0.1

    elif voicetype == 'high':
        f0min = 120
        f0max = 350
        cutoffhi = 150
        orderhi = 3
        cutofflo = 250
        orderlo = 5
        f0clipmin = 120
        f0clipmax = 350
        centreclip = 0.15

    else:
        f0min = int(cgiparams['f0min'])
        f0max = int(cgiparams['f0max'])
        cutoffhi = int(cgiparams['cutoffhi'])
        orderhi = int(cgiparams['orderhi'])
        cutofflo = int(cgiparams['cutofflo'])
        orderlo = int(cgiparams['orderlo'])
        f0clipmin = int(cgiparams['f0clipmin'])
        f0clipmax = int(cgiparams['f0clipmax'])
        centreclip = float(cgiparams['centreclip'])

    f0window = float(cgiparams['f0window'])
    f0winsamples = int(round(fs * f0window))
    f0stepfactor = int(cgiparams['f0stepfactor'])
    f0step = int(round(f0winsamples/f0stepfactor))

    f0median = int(cgiparams['f0median'])

    polydegree = int(cgiparams['polymodel'])
    if polydegree > 40:
        print "Polydegree ",polydegree," reset to 40."
        polydegree = 40
    if polydegree < 0:
        print "Polydegree ",polydegree," reset to 0."
        polydegree = 0

    voicewin = int(cgiparams['voicewin'])
    voicediff = int(cgiparams['voicediff'])

#==============================================================
# High and low pass filtering

    buttery = butterworth(signal,fs,cutoffhi,orderhi,cutofflo,orderlo)

#==============================================================
# S0FT F0 estimation with zero-crossing and peak-picking

    f0xvector = f0x(buttery,fs,f0winsamples,f0step,f0median,f0min,f0max,centreclip,f0clipmin,f0clipmax)

#===============
# Calculation of polynomial model

    f0xmed = np.median(f0xvector)
    f0xmedfill = [ x if x<f0clipmax and x>f0clipmin else f0xmed  for x in f0xvector ]

#    f0xnew = voicedetector(f0xvector,voicewin,voicediff)
#    f0xvector = np.array(f0xnew)

    f0xpoly = np.array(polyregline(range(len(f0xmedfill)),f0xmedfill,polydegree))

#==============================================================
# RAPT

    if opsys == "linux":

        f0raptvector = f0rapt(filebase,f0min,f0max,0.01,0.02)[0]

# Calculation of polynomial model

        f0raptmed = np.median(f0raptvector)
        f0raptmedfill = [ x if x<f0clipmax and x>f0clipmin else f0raptmed for x in f0raptvector ]
#        f0raptpoly = np.array(polyregline(range(len(f0raptvector)),f0raptvector,polydegree))

#==============================================================

except:
    print "Parameter value error."
    exit()

#==============================================================
#==============================================================

inithtml()

# <h2 align=\"center\">S0FT - Simple F0 Tracker</h2>

metadata = cgiparams['metadata']

print '<font size="-1">'
print '<table align="center" valign="top"><tr valign="top"><td>'

print '<table align="center" valign="top">'
print '<tr valign="top"><td width="250"><b>',metadata,'</b></td></tr>'
print '</table>'

print '</td><td>'

#print '<table align="center" valign="top">'
#print '<tr><td>Main parameters:</td><td>',paramstringio,'</td></tr>'
#print '<tr><td>F0 processing:</td><td>',paramstringf0proc,'</td></tr>'
#print '<tr><td>Waveform processing:</td><td>',paramstringwaveproc,'</td></tr>'
# print '</table>'

# print '</td><td>'

print '<table align="center" valign="top">'
print '<tr><td>Sampling rate:</td><td>',fs,'Hz</td></tr>'
print '<tr><td>Signal length:</td><td>',signallength,'samples</td></tr>'
print '<tr><td>Signal duration:</td><td>',signalduration,'s</td></tr>'
print '</table>'

print '</td><td>'

print '<table align="center" valign="top">'
print '<tr><td>F0 sample window:</td><td>',f0window,'(s),',f0winsamples,'samples</td></tr>'
print '<tr><td>F0 window step:</td><td>factor',f0stepfactor,",", f0step,'samples</td></tr>'
print '<tr valign="top"><td>Voice selection:</td><td>"'+voicetype+'"</td></tr>'
print '</table>'

print '</td></tr>'

print '</font>'

#==============================================================================
#==============================================================================

figwidth = int(cgiparams['figwidth'])
figheight = int(cgiparams['figheight'])

if opsys == "linux":
    figrows = 7; figcols = 1
else:
    figrows = 5; figcols = 1

fig = plt.figure(1,figsize=(figwidth,figheight))

fontsize = 18

#==============================================================================
# Waveform

rownum = 0 ; colnum = 0
rowspan = 1; colspan = 1
ax1 = plt.subplot2grid((figrows,figcols), (rownum,colnum),rowspan=colspan,colspan=colspan)
ax1.set_title("S0FT: F0 estimator ("+filebase+")",fontsize=fontsize)
x = [ 1.0*x/fs for x in range(len(signal)) ]
y = signal
ax1.set_yticks([])
ax1.set_xlim(-0.1,np.ceil(x[-1])+0.01)
plt.plot(x,y,color='lightblue')
# plt.plot(range(len(buttery)),buttery-5000,linewidth=4,color='r')

#==============================================================================
# F0X

rownum += rowspan  ; colnum = 0
rowspan = 2; colspan = 1
ax2 = plt.subplot2grid((figrows,figcols), (rownum,colnum),rowspan=rowspan,colspan=colspan)
ax2.set_title("S0FT",fontsize=fontsize)
y = f0xvector
x = [ (1.0 * x * f0window) / f0stepfactor for x in range(len(y))]
ax2.set_xlim(-0.1,np.ceil(x[-1])+0.01)
ax2.set_ylim((f0min,f0max))
plt.grid(color='b', linestyle='--', linewidth=0.5)
plt.scatter(x,y,color="b",s=4)

for xx,yy,pp,ss in zip(x,y,f0xpoly,y[1:]+1):
    if yy > f0clipmin and yy<f0clipmax:
        plt.scatter(xx,pp,color='r',s=0.5)

#==============================================================================
# F0 RAPT

if opsys == "linux":
    rownum += rowspan ; colnum = 0
    rowspan = 2; colspan = 1
    ax3 = plt.subplot2grid((figrows,figcols), (rownum,colnum),rowspan=rowspan,colspan=colspan)
    ax3.set_title("RAPT", fontsize=fontsize)
    y = f0raptvector
    x = [ x * 0.01 for x in range(len(y)) ]
    ax3.set_xlim(-0.1,np.ceil(x[-1])+0.01)
    ax3.set_ylim((f0min,f0max))
    plt.grid(color='b', linestyle='--', linewidth=0.5)
    plt.scatter(x,y,s=3,color='g')

#==============================================================================
# Spectrogram

specmin = 0
specmax = 3000
nykvistspect = int(round(fs/2.0))
if nykvistspect < 3000:
    specmax = nykvistspect
specwinfactor = 0.005
specwin = int(float(65536*specwinfactor))

rownum += rowspan ; colnum = 0
rowspan = 2; colspan = 1
ax4 = plt.subplot2grid((figrows,figcols), (rownum,colnum),rowspan=rowspan,colspan=colspan)
ax4.set_title("Spectrogram", fontsize=fontsize)

plt.ylabel('Freq [Hz]',fontsize=fontsize)

NFFT = specwin
ax4.specgram(signal, NFFT=NFFT, Fs=fs)
plt.axis(ymin=specmin, ymax=specmax)

#==============================================================================

plt.tight_layout()

if opsys == "linux":
	figlocaladdress = '/var/www/Webtmp/'+filebase+'.png'
	figwebaddress = '/Webtmp/'+filebase+'.png'
	audioaddress = wavfilename
else:
	figlocaladdress = '../../docs/Webtmp/'+filebase+'.png'
	figwebaddress = '/gibbon/Webtmp/'+filebase+'.png'
	audioaddress = '/gibbon/S0FT/'+wavfilename

plt.savefig(figlocaladdress)

print '<br><table align="center" width="100%"><tr align="center"><td colspan="2">'

print '<audio controls style="width: 90%;">'
print '<source src="'+audioaddress+'" type="audio/wav" length="100%">'
print 'Your browser does not support the audio element.</audio>'
print '</td></tr>'

print '<tr align="left"><td width=5%></td><td width=100%>'

print '<img src='+figwebaddress+' width="85%">'

"""
# mpld3 has very poor resolution of the spectrogram

if opsys == "linux":
#	print '<a href=\"' + figwebaddress + '\">Downloadable figure</a>'
	htmlstr = mpld3.fig_to_html(fig)
	print htmlstr
else:
	print '<img src='+figwebaddress+' width="100%">'
"""
print '</td></tr></table>'

plt.close()

#==============================================================================

print "</body></html>"

#==================================================================
# EOF =============================================================
