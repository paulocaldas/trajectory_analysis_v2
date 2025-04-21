# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:24:39 2020

@author: pcaldas
"""

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

def angle(v1,v2):
        ''' calculates the angle between two vectors between 0 and 2pi 
        v1 and v2 are arrays with the form np.array([a,b]) '''
    
        #dot_prod = np.dot(v1,v2) 
        #norm_prod = np.linalg.norm(v1) * np.linalg.norm(v2)
        #cos_tetha = dot_prod/norm_prod    
        
        signed_angle = math.atan2(v2[1],v2[0]) - math.atan2(v1[1],v1[0])
        
        return signed_angle # in rads
    
def angle_correlation(X,Y):
    
    XY_df = pd.DataFrame([X,Y]).T
    
    angles = []
      
    #moving step for pairwise comparison
    
    for step in np.arange(1, int(int(len(XY_df)) * 0.5)):

        XY_table = XY_df[::step]
        
        step_angles = []
        
        for i in np.arange(1, len(XY_table) - 1):
            
            a = XY_table.iloc[i].values - XY_table.iloc[i-1].values
            b = XY_table.iloc[i+1].values - XY_table.iloc[i].values
                
            angles.append(angle(a,b))
        
        angles.append(np.mean(step_angles))
    
    return angles

def plot_angles(angles, polar_ax = None, bins = 24, cmap = plt.cm.Blues):
    '''angles in radians'''
    
    n_bins = bins
    
    angle_bins  = np.linspace(0,2 * np.pi, n_bins + 1, endpoint = True) # divide angles (0 to 360 degrees) into 8 bins
    counts, bins = np.histogram(angles, angle_bins, density = True) #split my data into bins

    width = 2 * np.pi / n_bins # width of each slice/bins
    
    # by default create a new plot
    if polar_ax == None:
        plt.figure()
        axes = plt.subplot(projection='polar')
    
    # to add to existent axis set polar_ax = fig.add_subplot(1, 2, 2, projection="polar") 
    
    else:
        axes = polar_ax
    
    bars = axes.bar(angle_bins[:n_bins], counts, width = width, bottom = 0.1, color = cmap(counts), 
                  edgecolor = 'black', alpha = 0.8, align = 'edge', )
    
    # y axis corresponds to counts on a normal histogra, proportional to the number of counts in each bin
    axes.set_yticklabels([])
    #ax.set_yticks([0, np.max(counts)])
    
    # set the lable go clockwise and start from the top
    axes.set_theta_zero_location("N")
    axes.set_theta_direction(-1)  