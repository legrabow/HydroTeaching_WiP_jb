#!/usr/bin/env python
# coding: utf-8

# <img src='../IMG/TUD_logo.png' align='right' width='15%'></img>
# 
# [Zurück zur Übersicht MWW01 - Grundwasserbewirtschaftung mit Computermodellen](../MWW01_00_index.ipynb)   
# [Zurück zur Übersicht T02 - Konzeptionelle Modelle](./MWW01_T02_00_index.ipynb)

# Developed by: Thomas.Reimann@tu-dresden.de / Sophie.Pfoertner@mailbox.tu-dresden.de / Anne.Pfoertner@mailbox.tu-dresden.de / Leonard.Grabow@mailbox.tu-dresden.de
# <br>Last change: 2022 05 16
# <br>Current state: funktional
# <br>ToDo:
# 
# - Code optimieren ? (TR)
# - Skalierung der RMSE Box im Plot

# <center>
#     
# ### MWW01 - Grundwasserbewirtschaftung mit Computermodellen
# ### Thema 02: Konzeptionelle Modelle
# # Physikalische und Hydraulische Randbedingungen  </center>
# Dieses Notebook erläutert den Unterschied zwischen physikalischen und hydraulischen Randbedingungen anhand eines einfach Beispiels. Mit den zugehörigen Aufgaben kann das Verhalten der einzelnen Randbedingungen näher untersucht werden.
# 
# ### Lernziele:
# Nachdem Sie dieses Notebook erfolgreich bearbeitet haben, können Sie:
# * Die Unterschiede zwischen hydraulischen und physikalischen Randbedingungen für Grundwasserströmungsmodelle erläutern;
# * Disskutieren, mit welchen Möglichkeiten eine Grundwasserscheide für ein Grundwasserströmungsmodell berücksichtigt werden kann.
# 
# ### Ausgangssituation
# Betrachtet wird ein ungespannter, homogener und isotroper Grundwasserleiter, der durch zwei Standgewässer begrenzt wird. Der Grundwasserleiter wird zusätzlich durch Grundwasserneubildung gespeist. Die Grundwasserhydraulik kann durch eine 1D Vereinfachung der allgemeinen Grundwasserströmungsgleichung beschrieben werden, für die folgende analytische Lösung existiert (siehe nachfolgende Abbildung)
# 
# <br>
# <center>$\large h(x)=\sqrt{H_l^2-\frac{H_l^2-H_r^2}{L}x+\frac{q_n}{K}x(L-x)}$</center>
# <br>
# mit
# 
# * $h$ = Potentialhöhe
# * $x$ = Ortskoordinate entlang der Strömungsrichtung
# * $K$ = hydraulische Leitfähigkeit
# * $L$ = Distanz zwischen den beiden definierten Potentialhöhen
# * $H_l$ = definierte Potentialhöhe am linken Rand
# * $H_r$ = definierte Potentialhöhe am rechten Rand
# * $q_n$ = Grundwasserneubildung
# 
# <figure>
#   <img src="../FIG/MWW01_T02_1D_unconfied_BC_4.png" alt="Trulli" style="width:50%">
#   <figcaption>Abb.1 - Skizze der hydrogeologischen Situation.</figcaption>
# </figure>
# 
# <br> Im nachfolgenden Teil des Notebooks wird die analytische Lösung ermittelt und grafisch dargestellt.

# In[1]:


# Initialize librarys
from scipy.special import erfc, erf
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import json
from ipywidgets import *
import sys
# NOTE: Um einen relativen Import von JupyterQuiz zu ermoeglichen, sollte es in demselben Ordner
# sein wie dieses Script. Wenn das nicht der Fall ist, muss der absolute Pfad angegeben werden.
# Das Jupyter Quiz: https://github.com/jmshea/jupyterquiz.git von John Shea
sys.path.append("/home/grabow/Dropbox/SHK/Jupyter/HydroTeaching/jupyterquiz")
from jupyterquiz import display_quiz

# Read in dictionary
# NOTE: Auch hier sollte das Woerterbuch in demselben Ordner sein
with open("/home/grabow/Dropbox/SHK/Jupyter/HydroTeaching/dictionary.json", 'r') as f:
    language = json.load(f)


# In[2]:


# Choose language for plot and assign to strings
english = True    

if english:
    label = language["Berechnetes Druckpotential"]
    ylabel = language["Druckpotential"] + " (m)"
    title_ax1 = language["Potentialhöhe für die 1D ungespannte Grundwasserströmung"]
    title_ax2 = language["Lineare Regression"]
    gwn = language["GWN"]
else:
    label = "Berechnetes Druckpotential"
    ylabel = "Druckpotential" + " (m)"
    title_ax1 = "Potentialhöhe für die 1D ungespannte Grundwasserströmung"
    title_ax2 = "Lineare Regression"
    gwn = "GWN"


# In[3]:


def calculate_head(R, hr, hl, K, L, x):
    """
    Calculate the head after Depuit-Forchheimer-equation
    
    Keyword Arguments:
    R -- Recharge, m/s
    hr -- head at x=0, m
    hl -- head at x=L, m
    K -- hydraulic conductivity, m/s
    L -- distance to hl, m
    x -- points to get the head for, m
    """
    h=(hl**2-(hl**2-hr**2)/L*x+(R/K*x*(L-x)))**0.5
    return h


# Definition of the function
def head(hl, hr, L, R, K,y_scale):
    """
    Plot the head and correlation
    
    Keyword Arguments:
    R -- Recharge, mm/a
    hr -- head at x=0, m
    hl -- head at x=L, m
    K -- hydraulic conductivity, m/s
    L -- distance to hl, m
    y_scale -- scale of y, -
    """
    # Transform recharge units from mm/a to m/s
    R=R/1000/365.25/86400
    
    # Calculate head
    x = np.arange(0, L,L/1000)
    #global h
    h = calculate_head(R, hr, hl, K, L, x)
    
    # Prepare figure and subplots

    fig = plt.figure(figsize=(9,6))
    ax = fig.add_subplot(1, 1, 1)
    
    # Plot the head
    ax.plot(x,h, label=label)
    ax.set(xlabel='x (m)', ylabel=ylabel,title=title_ax1)
    ax.fill_between(x,0,h, facecolor='lightblue')
    
    # Plot the CH-boundaries
    ax.vlines(0, 0, hl, linewidth = 10, color='b')
    ax.vlines(L, 0, hr, linewidth = 10, color='b')
    
    # Make Water-Triangle
    y_range = abs((hl*(1-y_scale/100))-(hr*(1+y_scale/100)))
    h_arrow = (hl**2-(hl**2-hr**2)/L*(L*0.96)+(R/K*(L*0.96)*(L-(L*0.96))))**0.5  #water level at arrow
    ax.arrow(L*0.96,(h_arrow+(h_arrow*0.004)), 0, -0.01, fc="k", ec="k", head_width=(L*0.015), head_length=(y_range*0.03))
    #ax.hlines(y= h_arrow-(h_arrow*0.0005), xmin=L*0.95, xmax=L*0.97, colors='blue')   
    ax.hlines(y= h_arrow-(y_range*0.01), xmin=L*0.95, xmax=L*0.97, colors='blue') 
    #ax.hlines(y= h_arrow-(h_arrow*0.001), xmin=L*0.955, xmax=L*0.965, colors='blue')
    ax.hlines(y= h_arrow-(y_range*0.015), xmin=L*0.955, xmax=L*0.965, colors='blue')

    
    # Add arrows for recharge
    if R != 0:
        head_length=(R*1000*0.0005 * (86400*365.25))*y_range # DIVISION STATT MULTIPLIKATION?
        h_rch1 = (hl**2-(hl**2-hr**2)/L*(L*0.25)+(R/K*(L*0.25)*(L-(L*0.25))))**0.5  #water level at arrow for Recharge Posotion 1
        ax.arrow(L*0.25,(h_rch1+head_length), 0, -0.01, fc="k", ec="k", head_width=(L*0.03), head_length=head_length)
        h_rch2 = (hl**2-(hl**2-hr**2)/L*(L*0.50)+(R/K*(L*0.50)*(L-(L*0.50))))**0.5  #water level at arrow for Recharge Postition 2
        ax.arrow(L*0.50,(h_rch2+head_length), 0, -0.01, fc="k", ec="k", head_width=(L*0.03), head_length=head_length)
        h_rch3 = (hl**2-(hl**2-hr**2)/L*(L*0.75)+(R/K*(L*0.75)*(L-(L*0.75))))**0.5  #water level at arrow for Recharge Position 3
        ax.arrow(L*0.75,(h_rch3+head_length), 0, -0.01, fc="k", ec="k", head_width=(L*0.03), head_length=head_length)

    # Add Grundwasserscheide
    max_y = max(h)
    max_x = x[h.argmax()]
    R_min_ms=K*abs(hl**2-hr**2)/L**2
    if R>R_min_ms:
        plt.vlines(max_x,0,max_y, color="r")
        
    
    # Set y min und max
    y_min = hl*(1-y_scale/100)
    y_max = hr*(1+y_scale/100)
    
    # Add flow arrows
    grad_h = np.diff(h)
    num_arows = 10

    max_arrow_length = 100
    start_arrow = 0 + 30 + max_arrow_length
    end_arrow = L - max_arrow_length
    steps = (end_arrow - start_arrow) / num_arows
    numeric_x = np.concatenate((np.arange(max_x,start_arrow, -steps), np.arange(max_x, end_arrow, steps)))
    numeric_x = np.delete(numeric_x, np.argwhere(numeric_x == max_x))
    
    sample_x = [np.absolute(x-x_val).argmin() for x_val in numeric_x]
    max_grad = np.max(abs(grad_h[sample_x]))
    
    
    puffer = (y_max - y_min) * 0.02

    for idx in sample_x:
        arrow_length = max_arrow_length * grad_h[idx] / max_grad
        x_pos = x[idx]
        y_positions = np.linspace(y_min + puffer, h[idx]  - puffer, num_arows)
        for y_pos in y_positions:
            #print(f"plot arrow at {x_pos}, {y_pos}")
            plt.arrow(x_pos, y_pos, -arrow_length, 0, width=0.1, head_width=0.2,head_length=y_scale*3.,length_includes_head=True)

    plt.ylim(y_min, y_max)
    plt.xlim(-50,L+50)
    plt.text(L, (hr*(1+y_scale/100))-0.1*y_range, gwn+r': {:.3e} m/s '.format(R), horizontalalignment='right', bbox=dict(boxstyle="square", facecolor='azure'),fontsize=12)
    ax.grid()
    
    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.show()


# In[4]:


# Start interactive plot

interactive_plot = interact(head,
         y_scale = widgets.BoundedFloatText(value=5, min=1, max=100, step=1, description='Skal. y-Achse:', disabled=False),
         hl=widgets.BoundedFloatText(value=100, min=0, max=1000, step=1, description='H_l:', disabled=False),
         hr=widgets.BoundedFloatText(value=102, min=0, max=1000, step=1, description='H_r:', disabled=False),
         L= widgets.BoundedFloatText(value=2500,min=0, max=20000,step=100, description='L:' , disabled=False),
         R=widgets.FloatSlider(value=0,min=-500, max=500, step=10,description=gwn+":", disabled=False),
         K=widgets.FloatLogSlider(value=0.0001,base=10,min=-6, max=-2, step=0.1,description='K:',readout=True,readout_format='.2e'))


# In[5]:


display_quiz("/home/grabow/Dropbox/SHK/Jupyter/HydroTeaching/questions_unconfined_flow.json")


# <hr>
# &copy; 2022 | Thomas Reimann
# <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img style="float: right" alt="Creative Commons Lizenzvertrag" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a>
