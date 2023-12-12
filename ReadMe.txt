1- Run "main.m" file in MATLAB.
2- Change the dialogbox defualt values if you want to run for a different environment.
	Number of rays: The number of optical rays that will be used for simulation.
	Refraction index: The refraction index of the prism.
	Simulation, on=1, off=0: If it set 1, the process of simulation will be shown graphically (but in expense of lower speed).
	        Only 1/1000 of simulated rays are shown graphically.
	Reflection: The reflection coefficient of absorbing surfaces.
	Apex_0: The smallest apex angle which is used for the starting of simulation. (Simulation will be performed for prisms 
		with apex angles between Apex_0 and 180 degrees)
	Angular distribution: It defines the initial angular distribution for the radiation source. 1= Lambertian, according to 
		the Lambert's cosine law (for natural termal radiation), 2= Uniformly distributed radiation, 3= directional radiation.
		In this case you will asked for the polar and the azimuthal angles for the ray directions.
3- Click the "ok" button to start the simulation.

"Ray_tracer", "tra" and "unit_cell" are subroutines.