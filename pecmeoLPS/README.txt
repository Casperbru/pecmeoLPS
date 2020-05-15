****************************
LPS SIMULATION AND ORBIT DETERMINATION WITH ISL
Created by Casper Brunt 2020
c.c.brunt@student.tudelft.nl

****************************

In this repository the tools are found to simulate satellite constellations and perform position estimations on a constellation.
The tools enable you to make use of the capability of the constellation to use inter-satellite linking to enhance the positioning
accuracy. The tool is developed for the master thesis "Assessment of the LPS Ephemeris Accuracy using Inter-Satellite Linking".

In general, the main file is the only file that has to be modified to perfrom different simulations. 
The source code is provided with comments explaining the steps set.

The main file:
- PECMEORUN: this is the fail where all the settings for the simulation and estimation can be set, here it 
		is also determined how the results are processed and saved. Here all the other functions are
		accessed.

There are four separate tool files: 
- constellationLoader: in this file the data of the constellations is stored and the initial conditions are
			produced.

- toolsGeneral: contains functions that call other functions, creates the orbits and observations simulation,
		calls all the estimation functions, used for processing of data

- toolsSimulation: function where the simulation is completely done, from constellation orbits to observations

- toolsNavigation: this is the file that contains all the navigation estimation algorithms, within this file
			the DOP values, SPP estimations and KOD estimations are computed

There is one verification file
- verificationRun: this is a very minimal file. It only test SPP by using data from a numerical example