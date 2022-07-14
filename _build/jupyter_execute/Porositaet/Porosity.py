#!/usr/bin/env python
# coding: utf-8

# # Porosity

# In[1]:


import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import numpy as np
import json
from ipywidgets import *
from IPython.display import display,clear_output
from matplotlib.patches import Patch


# In[2]:


def plot_box(portion_solid, portion_water):
    """
    Plot a box with solid, water and air in it
    
    Keyword Arguments:
    portion_solid -- portion of total height to plot for solid
    portion_water -- portion of total height to plot for water
    """
    
    
    # Prepare plot
    fig = plt.figure(figsize=(7,7))
    ax = plt.subplot(111)
    
    # Increase axes width and equalize plot
    ax.set_aspect("equal")
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(3)
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)

    # Define dummy boundaries for plot
    sandbox_height = 1
    sandbox_width = 1
    ax.set_xlim(0, sandbox_width)
    ax.set_ylim(0, sandbox_height)
    
    # Plot solid
    ax.add_patch(Rectangle((0,0), sandbox_width, portion_solid,
                 color = "gold"))
    
    # Plot water
    ax.add_patch(Rectangle((0,portion_solid), sandbox_width, portion_water,
                 color = "skyblue", alpha=0.8))
    
    # The rest will be air...
    
    # Plot legend
    legend_elements = [Patch(facecolor='white', edgecolor='black',
                         label='Air'),
                   Patch(facecolor='skyblue', edgecolor='black',
                         label='Water'),
                       Patch(facecolor='gold', edgecolor='black',
                         label='Solid')]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.1, 0.6),  prop={'size': 15})
    
    
    plt.show()


# In[3]:


# Define the sliders for the portions of solid and water
portion_solid = widgets.FloatSlider(value=0.3, min=0, max=1.0, step=0.01, description='Portion of Solid in %',
                                   continuous_update = False, style= {'description_width': 'initial'})

portion_water = widgets.FloatSlider(value=0.4, min=0, max=1.0, step=0.01, description='Portion of Water in %',
                                   continuous_update = False, style= {'description_width': 'initial'})

# Define the field for the question and the answer
question = widgets.FloatText(
    value=0,
    description='What will the total porosity be?',
    style= {'description_width': 'initial'})

answer = widgets.Valid(
    value=None,
    description='',
    disabled = True,
    style= {'description_width': 'initial'}
)

# Tell not to show the answer before the question has not been answered
answered = False

# Couple the portion of solid and water such that they are never greater than one together
def change_solid(change):
    """
    Adjust the portion of solid once the portion of water and solid togehter is greater than one
    
    Keyword Arguments:
    change -- ipywidgets listening object when the water portion is adjusted
    """
    if change["new"] + portion_solid.value > 1:
        portion_solid.value = 1 - change["new"]
def change_water(change):
    """
    djust the portion of water once the portion of water and solid togehter is greater than one
    
    Keyword Arguments:
    change -- ipywidgets listening object when the solid portion is adjusted
    """
    if change["new"] + portion_water.value > 1:
        portion_water.value = 1 - change["new"]
portion_solid.observe(change_water, names='value')
portion_water.observe(change_solid, names='value')
        
# Couple the question input with the answer output
def check_answer(change):
    """
    Check the inputed answer against the calculated answer and update the answer field
    
    Keyword Arguments:
    change -- ipywidgets listening object when the answer is inputed
    """
    global answered
    tot_poro = 1 - portion_solid.value
    if not answered:
        display(answer)
        answered = True
    if float(change["new"]) == tot_poro:
        answer.description = "Right!"
        answer.value = True
    else:
        answer.description = f"Wrong! Try again..."
        answer.value = False
question.observe(check_answer, names='value')

# Start the plot and the question/answer field



# In[4]:


interact(plot_box,
         portion_solid = portion_solid,
         portion_water = portion_water)
display(question)

