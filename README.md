## Description
Project for Robot Control class (MSc Robotics Engineering). The project models and simulates a legged robot whose legs are delta robots. The aim is to control the system to do follow a certain number of waipoints. 

## Folder Contents:
1. *dynamics* folder: contains functions for computing dynamics matrices.
2. *kinematics* folder: contains functions for homogeneous transformations. 
3. *utility* folder: contains functions for computing direct and inverse kinematics, and for 
                     exporting matrices to Simulink.
4. *init_parameters* fodler: contains all files needed to initialize parameters. 

- compute_deltino_dyn_parallel.m: this file has been used for computing kinematics, differential 
kinematics, dynamics, constraints of the robot. It uses Parallel Computing Toolbox to speed up the computations. 
- main.m: this file must be called before everything else, since it initializes and prepares the workspace for the Simulink scheme. 
- plot_deltino_*.m: after Simulink simulation, is possible to plot the robot using the two files.  
- deltino.slx: main Simulink scheme. 
