#!/usr/bin/env python
# coding: utf-8

# In[1]:


hl=100
hr=102
L=2500
R=0
K=0.0001
x = np.arange(0, L,L/1000)

fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(2, 1, 1)
(ln,) = ax.plot(np.arange(0, L,L/1000),calculate_head(R, hr, hl, K, L, x), animated=True)
plt.show(block=False)
plt.pause(0.1)
bg = fig.canvas.copy_from_bbox(fig.bbox)
# draw the animated artist, this uses a cached renderer
ax.draw_artist(ln)
# show the result to the screen, this pushes the updated RGBA buffer from the
# renderer to the GUI framework so you can see it
fig.canvas.blit(fig.bbox)

# Definition of the function
def head_faster(hl, hr, L, R, K,y_scale, english, add_observation_data):
    x = np.arange(0, L,L/1000)
    h = calculate_head(R, hr, hl, K, L, x)
    
    # PLOT FIGURE

    if english:
        label = language["Berechnetes Druckpotential"]
        label_syn = language["Beobachtetes Druckpotential"]
        ylabel = language["Druckpotential"] + " (m)"
        title = language["Potentialhöhe für die 1D ungespannte Grundwasserströmung"]
        title_ax2 = language["Lineare Regression"]
    else:
        label = "Berechnetes Druckpotential"
        label_syn = "Beobachtetes Druckpotential"
        ylabel = "Druckpotential" + " (m)"
        title = "Potentialhöhe für die 1D ungespannte Grundwasserströmung"
        title_ax2 = "Lineare Regression"

        
    fig.canvas.restore_region(bg)
    # update the artist, neither the canvas state nor the screen have changed
    ln.set_ydata(calculate_head(R, hr, hl, K, L, x))
    # re-render the artist, updating the canvas state, but not the screen
    ax.draw_artist(ln)
    
    
    #ax.set(xlabel='x (m)', ylabel=ylabel,title=title)
    #ax.fill_between(x,0,h, facecolor='lightblue')
    
    # BOUNDARY CONDITIONS hl, hr
    #ax.vlines(0, 0, hl, linewidth = 10, color='b')
    #ax.vlines(L, 0, hr, linewidth = 10, color='b')
    
    # MAKE 'WATER'-TRIANGLE
    y_range = abs((hl*(1-y_scale/100))-(hr*(1+y_scale/100)))
    h_arrow = (hl**2-(hl**2-hr**2)/L*(L*0.96)+(R/K*(L*0.96)*(L-(L*0.96))))**0.5  #water level at arrow
    ax.arrow(L*0.96,(h_arrow+(h_arrow*0.004)), 0, -0.01, fc="k", ec="k", head_width=(L*0.015), head_length=(y_range*0.03))
    #ax.hlines(y= h_arrow-(h_arrow*0.0005), xmin=L*0.95, xmax=L*0.97, colors='blue')   
    ax.hlines(y= h_arrow-(y_range*0.01), xmin=L*0.95, xmax=L*0.97, colors='blue') 
    #ax.hlines(y= h_arrow-(h_arrow*0.001), xmin=L*0.955, xmax=L*0.965, colors='blue')
    ax.hlines(y= h_arrow-(y_range*0.015), xmin=L*0.955, xmax=L*0.965, colors='blue')

    
    #MAKE DATA SYNTHETIC DATA
    x_syn = np.linspace(0, L,20)
    h_syn = calculate_head(hl=hl, hr=hr, K=2*10**(-5), R=250, x=x_syn,L=L)
    ax.scatter(x_syn,h_syn, label=label_syn)
    
    #ARROWS FOR RECHARGE 
    if R != 0:
        head_length=(R*86400*365.25*1000*0.0005)*y_range
        h_rch1 = (hl**2-(hl**2-hr**2)/L*(L*0.25)+(R/K*(L*0.25)*(L-(L*0.25))))**0.5  #water level at arrow for Recharge Posotion 1
        ax.arrow(L*0.25,(h_rch1+head_length), 0, -0.01, fc="k", ec="k", head_width=(L*0.03), head_length=head_length)
        h_rch2 = (hl**2-(hl**2-hr**2)/L*(L*0.50)+(R/K*(L*0.50)*(L-(L*0.50))))**0.5  #water level at arrow for Recharge Postition 2
        ax.arrow(L*0.50,(h_rch2+head_length), 0, -0.01, fc="k", ec="k", head_width=(L*0.03), head_length=head_length)
        h_rch3 = (hl**2-(hl**2-hr**2)/L*(L*0.75)+(R/K*(L*0.75)*(L-(L*0.75))))**0.5  #water level at arrow for Recharge Position 3
        ax.arrow(L*0.75,(h_rch3+head_length), 0, -0.01, fc="k", ec="k", head_width=(L*0.03), head_length=head_length)

    
    #Grundwasserscheide
    max_y = max(h)
    max_x = x[h.argmax()]
    R_min_ms=K*abs(hl**2-hr**2)/L**2
    if R>R_min_ms:
        plt.vlines(max_x,0,max_y, color="r")

    plt.ylim(hl*(1-y_scale/100),hr*(1+y_scale/100))
    plt.xlim(-50,L+50)
    plt.text(L, (hr*(1+y_scale/100))-0.1*y_range, r'GWN: {:.3e} m/s '.format(R), horizontalalignment='right', bbox=dict(boxstyle="square", facecolor='azure'),fontsize=12)
    ax.grid()
    plt.legend(loc="upper left")
    
    
    # RESIDUALS
    h_inter = np.interp(x_syn, x, h)
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.scatter(h_syn, h_inter)
    ax2.set(xlabel=label_syn + " (m)", ylabel=label + " (m)",title=title_ax2)
    
    # copy the image to the GUI state, but screen might not be changed yet
    fig.canvas.blit(fig.bbox)
    # flush any pending GUI events, re-painting the screen if needed
    fig.canvas.flush_events()
    #print('R: ',R, ' m/s')

