# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 14:04:10 2020

@author: pcaldas
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# multicolor lines functions
def multicolored_lines(x,y, cmap = 'plasma', linewidth = 1.5, alpha = 1.0, bar = True):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    """

    #fig, ax = plt.subplots()
    lc = colorline(x, y, cmap = cmap, linewidth = linewidth, alpha = alpha)
	
    if bar == True: plt.colorbar(lc)
    
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(), y.max())
    
    plt.axis('auto')
    plt.show()

def colorline(x, y, cmap='viridis', linewidth = 1, alpha = 1, norm = plt.Normalize()):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]
    z = np.linspace(0.0, 1.0, len(x))

    segments = make_segments(x, y)
    lc = LineCollection(segments,  array = z, cmap = cmap, norm = norm,
                              linewidth = linewidth, alpha = alpha)
    ax = plt.gca()
    ax.add_collection(lc)

    return lc

def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments
