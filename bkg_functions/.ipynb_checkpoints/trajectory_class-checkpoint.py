# -*- coding: utf-8 -*-
"""
Created on Fri May 29 12:23:08 2020
@author: pcaldas
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import FormatStrFormatter
import math 

class Trajectory:
    
    ''' X,Y position over time: X and Y are 1D arrays (e.g pixels or microns)
    step_rate: interval between time points (e.g in seconds);
    clip_fit: % of the MSD curve to perform the fit/analysis '''
    

    # Initializer / Instance Attributes
    def __init__(self, X, Y, step_rate = 1, clip_fit = 0.5):
        
        self.X = X
        self.Y = Y
        self.step_rate = step_rate
        self.clip_fit = clip_fit
        
        self.n_points = len(self.X)
        
        # attributes computed from the msd analysis
        self.msd_curve = self.msd_analysis()[0:2]
        self.msd_fit   = self.msd_analysis()[2:4]
        self.msd_alpha = self.msd_analysis()[4]
        self.msd_param = self.msd_analysis()[5]
        self.msd_mode  = self.msd_analysis()[6]
        
    # instance method to compute total length of the track
    
    def distance(self):     
        """computes total lenght of the track"""
        
        XY_table = pd.DataFrame([self.X, self.Y]).T
        
        steps = []
    
        for frame in range(1,len(XY_table)):
            step_dist = np.linalg.norm(XY_table.iloc[frame] - (XY_table.iloc[frame-1]))
            steps.append(step_dist)
            
        return np.sum(steps)
    
    # instance method to compute displacement 
    
    def displacement(self):
        """computes euclidean distance from the initial position to the final position"""
        
        XY_table = pd.DataFrame([self.X, self.Y]).T
        return np.linalg.norm(XY_table.iloc[0] - XY_table.iloc[-1])
    
    # instance method to compute displacement 
    
    def directionality(self):
        """ratio between total distance and displacement"""
        return self.displacement() / self.distance()
    
    # instance method to compute total length of the track
    
    def step_velocity(self):
        """computes instant velocity of the track"""
        
        XY_table = pd.DataFrame([self.X, self.Y]).T
        
        velo = (XY_table - XY_table.shift(1)).dropna()
        instant_velo = np.sqrt((velo ** 2).sum(1)) / self.step_rate

        return instant_velo.mean()
    
    # instance method to compute MSD analysis for one trajectory
    
    def msd_analysis(self):
        """Compute MSD for one trajectory and fit a model to the curve
        depending on the type of behavior (directed, random or anamolous motion)
        returns: taus, msd_curve, taus_fit, msd_fit, alpha, msd_params, mode"""

        XY_table = pd.DataFrame([self.X, self.Y]).T
        
        msd_curve = []

        n_shifts = len(XY_table)

        for shift in range(1, n_shifts):
            diffs = abs(XY_table - XY_table.shift(-shift))
            msd = np.square(diffs.dropna()).sum(axis=1)
            msd_curve.append(msd.mean())  

        taus = np.array(range(1,n_shifts)) * self.step_rate

        # fit linear equation to find the scalling exponent (alpha)
        
        np.insert(taus, 0, 0)
        np.insert(msd_curve, 0,0)

        clip = int(len(msd_curve) * self.clip_fit)

        eq_linear = lambda t, alpha, d: alpha * t + d
        (alpha, d), cov_lin = curve_fit(eq_linear, np.log(taus[:clip]), np.log(msd_curve[:clip]),
                                        p0 = (0.5,0.5), bounds= ([0, -np.inf],[np.inf,np.inf]))

        # fit to a parabola with a diffusion (D) and velocity (V) component if alpha >> 1 (directed motion)
        
        if alpha > 1.2:
            
            mode = 'active transport'
            
            eq_parabola = lambda t, D, V: 4*D*t + (V*t)**2
            msd_params, cov_vel = curve_fit(eq_parabola, taus[:clip], msd_curve[:clip], p0 = (0.5,0.5),
                                                                 bounds= ([0, 0],[np.inf,np.inf]))
            # return attributes for active transport
            taus_fit = np.linspace(0, taus[-1], 500)
            msd_fit = eq_parabola(taus_fit, *msd_params)        
            
            msd_params = 'lateral diffusion = {:4.2} ; velocity = {:4.2}'.format(msd_params[0], msd_params[1])
            
        # fit to a linear equation with only a diffusion (D) component if alpha is around 1 (Brownian motion)
        elif alpha > 0.6:
            
            mode = 'random diffusion'
            
            eq_linear2 = lambda t, D: 4*D*t
            msd_params, cov_dif = curve_fit(eq_linear2, taus[:clip], msd_curve[:clip], p0 = 0.5, bounds= (0, np.inf))    
            
            # return attributes for Brownian diffusion
            taus_fit = np.linspace(0, taus[-1], 500)
            msd_fit = eq_linear2(taus_fit, *msd_params)
            
            msd_params = 'diffusion coefficient = {:4.2}'.format(msd_params[0])
            
        # fit to a linear equation with only a diffusion (D) component if alpha is around 1 (Brownian motion)
        
        else:
            mode = 'anommalous diffusion'
            
            eq_confined = lambda t,D,R: (R**2) * (1 - np.exp(-(4*D*t)/(R**2)))
            msd_params, cov_conf = curve_fit(eq_confined, taus[:clip], msd_curve[:clip],
                                            p0 = (0.5,0.5), bounds= ([0,0],[np.inf, np.inf]))    
            
            # return attributes for confned diffusion
            taus_fit = np.linspace(0, taus[-1], 500)
            msd_fit = eq_confined(taus_fit, *msd_params)
            
            msd_params = 'cage size = {:4.2} ; cage diffusion = {:4.2}'.format(msd_params[0], msd_params[1])
           
        return taus, msd_curve, taus_fit, msd_fit, alpha, msd_params, mode
    
    def msd_plot(self, colors = ['gray','red']):
        fig, ax = plt.subplots(figsize = (4,3), dpi = 120)
        ax.plot(*self.msd_curve, '--o', c = colors[0], alpha = 0.8, markeredgecolor = 'black', label = 'msd_curve')
        ax.plot(*self.msd_fit, '-', c = colors[1], label = 'msd_fit')
        ax.legend(frameon = False, loc = 0)
        ax.set_ylabel('$<MSD(\u03C4)>$'); ax.set_xlabel('\u03C4')
        ax.tick_params(direction = 'inout', width = 1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.show()
        
    def plot_trajectory(self, ax = None):
        ''' add to an existent axis (ax) if wanted'''
        
        if ax == None:
            fig, ax = plt.subplots(figsize = (4,3), dpi = 120)
            ax.plot(self.X, self.Y, '-o', lw = 2, color = 'gray', 
                     markeredgecolor = 'black', alpha = 0.8, label = self.msd_mode);
            ax.set_xlabel('X'); ax.set_ylabel('Y')
            ax.legend(frameon = False)
            
            ax.axis('off')
			#ax.spines['right'].set_visible(False)
            #ax.spines['top'].set_visible(False)   
            
            ax.tick_params(direction = 'inout', width = 1)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'));
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'));
            plt.show()
            
        else:
            ax.plot(self.X, self.Y, '--o', lw = 2, color = 'steelblue', 
                    markeredgecolor = 'black', alpha = 0.8, label = self.msd_mode)
            ax.set_xlabel('X'); ax.set_ylabel('Y')
            ax.legend(frameon = False)

    # under construction    
    def plot_angle_ori(self, cmap = plt.cm.Blues):
        
        
        def angle(v1,v2):
            ''' calculates the angle between two vectors between 0 and 2pi 
            v1 and v2 are arrays with the form np.array([a,b]) '''
        
            #dot_prod = np.dot(v1,v2) 
            #norm_prod = np.linalg.norm(v1) * np.linalg.norm(v2)
            #cos_tetha = dot_prod/norm_prod    
            
            signed_angle = math.atan2(v2[1],v2[0]) - math.atan2(v1[1],v1[0])
            
            return signed_angle # in rads
    
        
        def angle_correlation(X, Y):
        
            ''''computes angle correlation along the filament (i.e pairwise
            comparison as a funciton of lag time) and plot angle distribution'''
        
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

        def plot_angles(angles, polar_ax = None, bins = 24, cmap = cmap):
            '''plot angles distribution (in radians) on a polar axis'''
            
            n_bins = bins
            
            angle_bins  = np.linspace(0,2 * np.pi, n_bins + 1, endpoint = True) # divide angles (0 to 360 degrees) into 8 bins
            counts, bins = np.histogram(angles, angle_bins, density = True) #split my data into bins
        
            width = 2 * np.pi / n_bins # width of each slice/bins
            
            # by default create a new plot
            if polar_ax == None:
                plt.figure(figsize = (4,4), dpi = 120)
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
                
        angles_dist = angle_correlation(self.X, self.Y)
        plot_angles(angles_dist)
            
        return 
    
    def plot_all(self):
        self.plot_trajectory();
        self.msd_plot();
        self.plot_angle_ori();
        
def ComputeParameterDistribution(traj_table, traj_rate, param = 'msd_alpha', clip = -1):
    ''' returns a ditrbution of a chosen parameter for all tracks
    
    traj_table (dataframe): dataframe containing all tracks
    traj_rate (float): time interval of each step along the tracks (e.g frame rate)
    param (string): 'msd_alpha', 'msd_directed_velocity', 'length', 'displacement', 'directionality', 'step_velocity'
    clip: number of tracks to analyze; set to -1 to analyze all tracks
    ----
    '''
    
    param_dist = [] #empty list to add values

    for track in traj_table.TRACK_ID.unique()[:clip]:
        
        print('Analyzing Track_ID ', track, end='\r')
        
        track_to_analyze = traj_table[traj_table.TRACK_ID == track]
        track_info = Trajectory(track_to_analyze.POSITION_X, track_to_analyze.POSITION_Y, step_rate = traj_rate)
    
        # get MSD alpha parameter from each track
        if param == 'msd_alpha':
            param_dist.append(track_info.msd_alpha)
        
        if param == 'msd_directed_velocity':
            
            if track_info.msd_mode == 'active transport':
                msd_vel = float(track_info.msd_param.split('velocity = ')[-1])
                param_dist.append(msd_vel)
            
        if param == 'length':
            param_dist.append(track_info.distance())
            
        if param == 'displacement':
            param_dist.append(track_info.displacement())
            
        if param == 'directionality':
            param_dist.append(track_info.directionality())
            
        if param == 'step_velocity':
            param_dist.append(track_info.step_velocity())
            
    return param_dist