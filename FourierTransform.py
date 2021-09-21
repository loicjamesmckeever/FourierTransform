# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 12:36:39 2021

@author: Loïc James McKeever
"""

import numpy as np
import scipy as sp
import sympy as smp
from skimage import color, io
from scipy.fft import fftfreq, fft, ifft, fft2, ifft2

from bokeh.layouts import column, row
from bokeh.models import CustomJS, Slider, HoverTool, Label, Button
from bokeh.plotting import figure, output_file, show, ColumnDataSource

## Discrete Fourier Transform
## Nyquist frequency and undersampling

#Define various parameters
T = 40 #seconds
N = 100 #measurements
t = np.linspace(0, T, N)
dt = np.diff(t)[0]

#Define a few specific frequencies
f1 = 20/(N*dt)
f2 = 10/(N*dt)
f3 = (10+5*N)/(N*dt)

#Create a few signals
x1 = np.sin(2*np.pi*f1*t) + 0.3 * np.sin(2*np.pi*f2*t) + 0.3*np.random.random(len(t))
x2 = np.sin(2*np.pi*f2*t) + 0.1*np.random.random(len(t))
x3 = np.sin(2*np.pi*f3*t) + 0.1*np.random.random(len(t))

#Plot in bokeh
source = ColumnDataSource(data=dict(t=t,x1=x1))
source2 = ColumnDataSource(data=dict(t=t, x2=x2, x3=x3))

#Plot x1 in it's own graph
plot1 = figure(title="Signal x1 over time", width=1200, height=600)
plot1.line('t','x1', source=source)
x_label = "Time in seconds"
y_label = "Signal amplitude"
plot1.xaxis.axis_label = x_label
plot1.yaxis.axis_label = y_label

#Fourier Tranform plot for x1
f = fftfreq(len(t), dt)[:N//2]
x1_FFT = np.abs(fft(x1)[:N//2])

source_FFT = ColumnDataSource(dict(f=f, x1_FFT=x1_FFT))

plot1_FFT = figure(title="Fourier transform of x1", width=1200, height=600)
plot1_FFT.line('f', 'x1_FFT', source=source_FFT)
x_label_FFT = 'Frequency in Hz'
y_label_FFT = 'Intensity'
plot1_FFT.xaxis.axis_label = x_label_FFT
plot1_FFT.yaxis.axis_label = y_label_FFT

#Plot x2 and x3 together
plot2 = figure(title="Signals x2 and x3 over time", width=1200, height=600)
plot2.line('t','x2', source=source2, line_color='red')
plot2.line('t', 'x3', source=source2, line_color='green')
plot1.xaxis.axis_label = x_label
plot1.yaxis.axis_label = y_label

#Fourier transform of x2 and x3 together
x2_FFT = np.abs(fft(x2))[:N//2]
x3_FFT = np.abs(fft(x3))[:N//2]

source_FFT2 = ColumnDataSource(dict(f=f, x2_FFT=x2_FFT, x3_FFT=x3_FFT))

plot2_FFT = figure(title="Fourier transform of x1", width=1200, height=600)
plot2_FFT.line('f', 'x2_FFT', source=source_FFT2, line_color='red')
plot2_FFT.line('f', 'x3_FFT', source=source_FFT2, line_color='green')
plot2_FFT.xaxis.axis_label = x_label_FFT
plot2_FFT.yaxis.axis_label = y_label_FFT

N_slider = Slider(start=50, end=1500, value=100, step=50, title="N")

callback = CustomJS(args=dict(source=source2, f2=f2, f3=f3, N_slider=N_slider), code = """
                    
                    const data = source.data;
                    var N = N_slider.value;
                    
                    var t = makeArr(0,40, N);
                    var x2 = [];
                    var x3 = [];
                    
                   
                    for (var i = 0; i < N; i++){
                            x2[i] = Math.sin(2*Math.PI*f2*t[i]) + 0.1*Math.random();
                            x3[i] = Math.sin(2*Math.PI*f3*t[i]) + 0.1*Math.random();
                    }
                    

                    function makeArr(startValue, stopValue, cardinality) {
                        var arr = [];
                        var step = (stopValue - startValue) / (cardinality - 1);
                        for (var i = 0; i < cardinality; i++) {
                                arr.push(startValue + (step * i));
                                }
                        return arr;
                    }
                    
                    data['x2'] = x2;
                    data['x3'] = x3;
                    data['t'] = t;
                    
                    source.change.emit();
                    
                    """)

N_slider.js_on_change('value',callback)

## makeArr function from user mhodges on Stackoverflow

#Define output file
output_file("FourierTranform.html", title="Fourier Transform")

#Define layout
layout = column(plot1, plot1_FFT, row(plot2, N_slider), plot2_FFT)

#Display file
show(layout)