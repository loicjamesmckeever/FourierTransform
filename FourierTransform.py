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
from bokeh.models import CustomJS, Slider, HoverTool, Label, Button, Div
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

plot2_FFT = figure(title="Fourier transform of x2 and x3", width=1200, height=600)
plot2_FFT.line('f', 'x2_FFT', source=source_FFT2, line_color='red')
plot2_FFT.line('f', 'x3_FFT', source=source_FFT2, line_color='green')
plot2_FFT.xaxis.axis_label = x_label_FFT
plot2_FFT.yaxis.axis_label = y_label_FFT

#Define slider object
N_slider = Slider(start=100, end=1500, value=100, step=10, title="N (number of measurements)")

#FFT for x2 and x3 with each step of N
x2_FFT_N = []
x3_FFT_N = []
f_N = []
for n in range(100, 1550, 10):
    t_N = np.linspace(0, T, n)
    dt_N = np.diff(t_N)[0]
    x2_N = np.sin(2*np.pi*f2*t_N) + 0.1*np.random.random(len(t_N))
    x3_N = np.sin(2*np.pi*f3*t_N) + 0.1*np.random.random(len(t_N))
    x2_FFT_N.append(np.abs(fft(x2_N))[:n//2])
    x3_FFT_N.append(np.abs(fft(x3_N))[:n//2])
    f_N.append(fftfreq(len(t_N), dt_N)[:n//2])

#Custom JavaScript for slider
callback = CustomJS(args=dict(source=source2, source_FFT=source_FFT2, f2=f2, f3=f3, N_slider=N_slider, x2_FFT_N=x2_FFT_N, x3_FFT_N=x3_FFT_N, f_N = f_N), code = """
                    
                    const data = source.data;
                    const data_FFT = source_FFT.data;
                    
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
                    
                    data_FFT['f'] = f_N[((N/10)-10)];
                    data_FFT['x2_FFT'] = x2_FFT_N[((N/10)-10)];
                    data_FFT['x3_FFT'] = x3_FFT_N[((N/10)-10)];
                    
                    source_FFT.change.emit();
                    """)

N_slider.js_on_change('value',callback)

## makeArr function from user mhodges on Stackoverflow

#Define output file
output_file("FourierTranform.html", title="Discrete Fourier Transform")

#HTML headers and paragrphs
header1 = Div(text = """
             <header> Fourier Transform </header>
             """
             , style = {'font-size': '300%', 'font-family':'Georgia, serif'})

header2 = Div(text = """
             <header> Introduction </header>
             """
             , style = {'font-size': '200%', 'font-family':'Georgia, serif'})

p1 = Div(text = """
         <div> 
         The Fourier transform is a method of decomposing a function that depends on space or time into
         a function depending on spatial of temporal frequency.  For example you can use a Fourier transform
         on an audio signal that varies in amplitude over time to determine the frequencies the audio signal
         is composed of. The transform is defined by the following equation: 
             </div>
         """
         , style = {'font-size': '200%', 'font-family':'Georgia, serif'})

equation1 = Div(text = """
          <math> 
          <mover><mi>f</mi><mo>&Hat;</mo></mover> (&xi;) = 
          <msubsup><mo>&Integral;</mo><mn>-&infin;</mn><mo>&infin;</mo></msubsup> f(x) <msup><mi>e</mi><mo>-2&pi;ix&xi;</mo></msup> <mi>dx,  (Eq. 1)</mi>
          _______________________________________________________________________________________________
          </math>
         """
         , style = {'font-size': '200%', 'text-align':'center'})

header3 = Div(text = """
             <header> Discrete Fourier Transform </header>
             """
             , style = {'font-size': '200%', 'font-family':'Georgia, serif'})

p2 = Div(text = """
         <div>The discrete Fourier transform is a version of the Fourier transform that takes a finite number 
         of equally spaced samples of a time dependant function and returns a finite set of equally spaced
         samples of the frequency dependant NEEDS EDITING.
         </div>
         """
         , style = {'font-size': '200%', 'font-family':'Georgia, serif'})

equation2 = Div(text = """
          <math> 
          <mover><mi>x</mi><mo>&Hat;</mo></mover> (<msub><mi>f</mi><mn>n</mn></msub>) = 
          <msubsup><mo>&sum;</mo><mn>k</mn><mo>N</mo></msubsup> <msub><mi>x</mi><mn>t</mn></msub> <msup><mi>e</mi><mo>-2&pi;i(fn)k&Delta;t</mo></msup> <mi>,  (Eq. 2)</mi>
          _______________________________________________________________________________________________
          </math>
         """
         , style = {'font-size': '200%', 'text-align':'center'})

p3 = Div(text = """
         <div>The discrete Fourier transform is a version of the Fourier transform that takes a finite number 
         of equally spaced samples of a time dependant function and returns a finite set of equally spaced
         samples of the frequency dependant NEEDS EDITING.
         </div>
         """
         , style = {'font-size': '200%', 'font-family':'Georgia, serif'})

#Define layout
layout = column(header1, header2, p1, equation1, header3, p2, equation2, p3, plot1, plot1_FFT, plot2, N_slider, plot2_FFT)

#Display file
show(layout)