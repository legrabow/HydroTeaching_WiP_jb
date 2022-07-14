#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from ipywidgets import *
import time


# In[2]:


# Download real time 
url = 'https://data.geobasis-bb.de/geofachdaten/Wasser/Grundwasser/gw_basis_mn.zip'
local_path = 'tmp/'

r = requests.get(url)
z = zipfile.ZipFile(io.BytesIO(r.content))
z.extractall(path=local_path) # extract to folder
filenames = [y for y in sorted(z.namelist()) for ending in ['dbf', 'prj', 'shp', 'shx'] if y.endswith(ending)] 


# In[2]:


from ipyleaflet import Map, basemaps, WidgetControl, Polyline, Marker
from ipywidgets import IntSlider, ColorPicker, jslink



def ymax_conf(Q, K, i, b):
    ymax = Q/(2.*K*np.abs(i)*b)
    return ymax

def x0_conf(Q, K, i, b):
    x0 = -Q/(2.*np.pi*K*i*b)
    return x0


def TSL_conf(Q, K, i, b, dist_iso, plot_coor):
    h2 = h1 - i * (x1 - x2)
    L = x1 - x2

    ymax = ymax_conf(Q, K, i, b)
    #x0 = x0_conf(Q, K, i, b)
    y = np.arange(-ymax+0.1, ymax-0.1, 0.01)
    x = y/(np.tan(2*np.pi*K*i*b*y/Q))
    
    
    h_iso = isohypsen_conf(h1, h2, x1, x2, K, Q, b, x_iso_tilted, y_iso_tilted)

    tsl_normal = np.vstack([x,y])
    tsl_tilted = rotation_matrix.dot(tsl_normal)
    
    final_x = tsl_tilted[0,:] + lon_brunnen
    final_y = tsl_tilted[1,:] + lat_brunnen
    return final_x, final_y


def isohypsen_conf(h1, h2, x1, x2, K, Q, b, x, y):
    h_r1 = 0
    h_iso = (h2-h1)*(x-x1)/(x2-x1) +h1 + Q * np.log(np.sqrt(x**2+y**2) / r1) / (2 * np.pi * K * b) + h_r1
    return h_iso

def update_lines(change):
    Q = change.new
    
    x,y = TSL_conf(Q, K, i, b, dist_iso, plot_coor)
    locations = [[sub_x, sub_y] for sub_x, sub_y in zip(x,y)]
    line.locations = locations
    m.remove_layer(line)
    m.add_layer(line)

    
    

K=6E-4
i=0.0017
dist_iso=1
plot_coor=False
continuous_update=False

    
    
theta = np.radians(29)
c, s = np.cos(theta), np.sin(theta)
rotation_matrix = np.array(((c, -s), (s, c)))
inv_theta = np.radians(-29)
c, s = np.cos(inv_theta), np.sin(inv_theta)
inv_rotation_matrix = np.array(((c, -s), (s, c)))

readout_format='.3f'
y_arrow_end = rotation_matrix.dot([[0],[2000]])
y_arrow_start = rotation_matrix.dot([[0],[-2000]])
x_arrow_start = rotation_matrix.dot([[-4000],[0]])
x_arrow_end = rotation_matrix.dot([[2000],[0]])


xscale = 1. * 2000 / (0.48 + 0.35)
yscale = 1. * 2000 / (0.125 + 0.265)
x_min = -4 * xscale
x_max = 1.18 * xscale
y_min = -1.07 * yscale
y_max = 1 * yscale
x_iso, y_iso = np.meshgrid(np.arange(x_min, x_max, 10), np.arange(y_min, y_max, 10))
iso_normal = np.stack([x_iso,y_iso], axis=2)
iso_tilted= np.tensordot(inv_rotation_matrix, iso_normal, axes=([1],[2]))
x_iso_tilted = iso_tilted[0,:,:]
y_iso_tilted = iso_tilted[1,:,:]

r1 = 10000


lon_brunnen = 13.71639
lat_brunnen = 51.88775

x1 = -200
x2 = 400
h1 = 61
b = 100

#layout = widgets.Layout(width='auto', height='40px')

x,y = TSL_conf(0.0001, K, i, b, dist_iso, plot_coor)

line = Polyline(
    locations=[[sub_y, sub_x] for sub_x, sub_y in zip(x,y)],
    color="green" ,
    fill=False
)

m = Map(center = (lat_brunnen,lon_brunnen), zoom =10)
m.add_layer(line)
marker = Marker(location=(lat_brunnen,lon_brunnen), draggable=False)
m.add_layer(marker);

#zoom_slider = IntSlider(description='Q:', min=0, max=0.1, value=0.05)
#zoom_slider.observe(recalculate, names='value')
#jslink((zoom_slider, 'value'), (m, 'zoom'))
#widget_control1 = WidgetControl(widget=zoom_slider, position='topright')
#m.add_control(widget_control1)


m

        


# In[ ]:




