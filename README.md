# trajectory_analysis_v2
compute multiple parameters from a 2D spatio-temporal trajectories

# trajectory_analysis
compute multiple parameters from a 2D trajectory based on MSD analysis

#### use `Trajectory` to construct an object from X,Y coordinates

`my_traj = Trajectory(X,Y, step_rate = 1, clip_fit = 0.5)`

- X,Y are coordenates over time (e.g pixels or microns);
- step_rate: interval between time points (e.g in seconds);
- clip_fit: % of the MSD curve to perform the fit/analysis;

#### then several properties/methods are available with a simple command

- `my_traj.distance()` 		total length of the trajectory
- `my_traj.displacement()` 	total displacement of the trajectory (euclidian between init and final position)
- `my_traj.directionlity()`	between 0 (completely random) and 1 (directed motion); ration between displacement and total distance
- `my_traj.step_velocity()` 	average of step velocities along the trajectory
- `my_traj.plot_trajectory()`	show 2D plot

- `my_traj.msd_alpha` 		returns the scalling factor (α) from MSD ~ t^α
- `my_traj.msd_mode` 		returns the type of motion of the trajectory based on MSD scalling factor value
- `my_traj.msd_param` 		returns the parameters from MSD fit, according to the type of motion:
	- if `α > 1`: directed motion, fit to `4Dt + (Vt)^2`; msd_param = [D,V]
	- if `α ~ 1`: Browniam motion, fit to `4Dt` ; msd_param = [D]
	- if `α < 1`: anommalous diffusion, fit to `(R^2) (1 - np.exp((4Dt)/(R^2)))` ; msd_param = [D, R]

- `my_traj.msd_curve` 		returns X and Y values to plot msd curve
- `my_traj.msd_fit`		returns X and Y values to plot msd fit
- `my_traj.msd_plot()`		plot msd curve and msd fit

- `my_traj.plot_angle_ori()`	plot angle distribution for directionality analysis
