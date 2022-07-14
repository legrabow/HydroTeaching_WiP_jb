#!/usr/bin/env python
# coding: utf-8

# <img src='../IMG/TUD_logo.png' align='right' width='15%'></img>
# 
# [Zurück zur Übersicht MWW01 - Grundwasserbewirtschaftung mit Computermodellen](../MWW01/T02/.ipynb)   
# [Zurück zur Übersicht T02 - Konzeptionelle Modelle](./MWW01_T02_00_index.ipynb)

# # Analytical solution for 1D unconfined flow with two defined head boundaries
# 
# ## Equations
# 
# $h(x)=\sqrt{H_l^2-\frac{H_l^2-H_r^2}{L}x+\frac{q_n}{K}x(L-x)}$
# 

# Developed by: Thomas.Reimann@tu-dresden.de / Sophie.Pfoertner@mailbox.tu-dresden.de / Anne.Pfoertner@mailbox.tu-dresden.de
# Last change: 2020 11 03
# Current state: Funktioniert und sieht schon gut aus
# ToDo:
# 
# - Gitter / Hilfslinien optional
# - Wasserscheide markieren / x-Wert plotten
# - Wasserdreieck skalieren (wenn y_scale variiert)
# - Textbox (recharge): position skalieren (wenn y_scale ändert)
# - Code optimieren (TR)
# - Text ergänzen (TR)

# In[1]:


# Initialize librarys
from scipy.special import erfc, erf
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from ipywidgets import *


# In[2]:


# Definition of the function
def head(hl, hr, L, R, K,y_scale):
    x = np.arange(0, L,L/1000)
    R=R/1000/365.25/86400
    h=(hl**2-(hl**2-hr**2)/L*x+(R/K*x*(L-x)))**0.5
    
    # PLOT FIGURE
    fig = plt.figure(figsize=(9,6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x,h)
    ax.set(xlabel='x', ylabel='head',title='Hydraulic head for 1D unconfined flow')
    ax.fill_between(x,0,h, facecolor='lightblue')
    
    # BOUNDARY CONDITIONS hl, hr
    ax.vlines(0, 0, hl, linewidth = 10, color='b')
    ax.vlines(L, 0, hr, linewidth = 10, color='b')
    
    # MAKE 'WATER'-TRIANGLE
    h_arrow = (hl**2-(hl**2-hr**2)/L*(L*0.96)+(R/K*(L*0.96)*(L-(L*0.96))))**0.5  #water level at arrow
    ax.arrow(L*0.96,(h_arrow+(h_arrow*0.002)), 0, -0.01, fc="k", ec="k", head_width=(L*0.015), head_length=(h_arrow*0.0015))
    ax.hlines(y= h_arrow-(h_arrow*0.0005), xmin=L*0.95, xmax=L*0.97, colors='blue')   
    ax.hlines(y= h_arrow-(h_arrow*0.001), xmin=L*0.955, xmax=L*0.965, colors='blue')

    #ARROWS FOR RECHARGE 
    if R != 0:
        head_length=(R*86400*365.25*1000*0.002)*y_scale/2
        h_rch1 = (hl**2-(hl**2-hr**2)/L*(L*0.25)+(R/K*(L*0.25)*(L-(L*0.25))))**0.5  #water level at arrow for Recharge Posotion 1
        ax.arrow(L*0.25,(h_rch1+head_length), 0, -0.01, fc="k", ec="k", head_width=(head_length*300/y_scale), head_length=head_length)
        h_rch2 = (hl**2-(hl**2-hr**2)/L*(L*0.50)+(R/K*(L*0.50)*(L-(L*0.50))))**0.5  #water level at arrow for Recharge Postition 2
        ax.arrow(L*0.50,(h_rch2+head_length), 0, -0.01, fc="k", ec="k", head_width=(head_length*300/y_scale), head_length=head_length)
        h_rch3 = (hl**2-(hl**2-hr**2)/L*(L*0.75)+(R/K*(L*0.75)*(L-(L*0.75))))**0.5  #water level at arrow for Recharge Position 3
        ax.arrow(L*0.75,(h_rch3+head_length), 0, -0.01, fc="k", ec="k", head_width=(head_length*300/y_scale), head_length=head_length)

    
    #Grundwasserscheide
    max_y = max(h)
    max_x = x[h.argmax()]
    R_min_ms=K*abs(hl**2-hr**2)/L**2
    if R>R_min_ms:
        plt.vlines(max_x,0,max_y, color="r")

    plt.ylim(hl*(1-y_scale/100),hr*(1+y_scale/100))
    plt.xlim(-50,L+50)
    plt.text(L, (hr*1.016), 'R: {} m/s '.format(R), horizontalalignment='right', bbox=dict(boxstyle="square", facecolor='grey'))
    ax.grid()
    plt.show()
    #print('R: ',R, ' m/s')


# In[3]:


# Computation

interact(head,
         y_scale = widgets.BoundedFloatText(value=3, min=1, max=100, step=1, description='y_scale:', disabled=False),
         hl=widgets.BoundedFloatText(value=150, min=0, max=1000, step=1, description='h left:', disabled=False),
         hr=widgets.BoundedFloatText(value=152, min=0, max=1000, step=1, description='h right:', disabled=False),
         L= widgets.BoundedFloatText(value=2500,min=0, max=20000,step=100, description='L:' , disabled=False),
         R=(-500,500,10),
         K=widgets.FloatLogSlider(value=0.0001,base=10,min=-6, max=-2, step=0.1,readout=True,readout_format='.2e'))


# <hr>
# &copy; 2021 | Thomas Reimann
# <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img style="float: right" alt="Creative Commons Lizenzvertrag" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a>
