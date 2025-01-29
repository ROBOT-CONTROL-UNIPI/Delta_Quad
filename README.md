# Legged Robot with Delta-Leg Kinematics

## Overview
This project was developed for the **Robot Control** course (MSc Robotics Engineering). It focuses on modeling and simulating a **legged robot** with **delta robots as legs**, forming closed kinematic chains instead of the conventional serial-leg configurations.

The kinematics and dynamics of the robot were modeled from scratch using MATLAB **Symbolic Toolbox**. The derived equations were then implemented in Simulink with quasi-velocities formulation. 

## Features
- **Kinematic and dynamic modeling** of the robot
- **Implementation in Simulink**, including Dahl model for tangential contact dynamics
- **Gait generation system** based on the Double Support Triangle (DST) method
- **Control and guidance systems** for waypoint tracking

## Getting Started
To run the simulations, ensure you have MATLAB with the required toolboxes installed.

```bash
# Clone this repository
git clone https://github.com/yourusername/your-repo.git
cd your-repo
```
TODO

## Authors
- [ngazzanelli](https://github.com/ngazzanelli)
- [AlbertoNobili](https://github.com/AlbertoNobili)

## References
- B. Chen, S. Li and D. Huang, "Quadruped robot crawl gait planning based on DST," Proceedings of the 33rd Chinese Control Conference, Nanjing, China, 2014, pp. 8578-8582, doi: 10.1109/ChiCC.2014.6896440

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
