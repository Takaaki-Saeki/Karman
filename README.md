# Karman

https://raw.github.com/wiki/Takaaki-Saeki/Karman/gif/pressure.gif

## Overview
This project is Karman vortex simulation around a square prism. Poisson equation is solved for the analysis of the pressure field, and Kawamura-Kuwahara scheme is used for the analysis of the velocity field. Visualization of pressure field is conducted by using python3.6 + matplotlib (details below).

## details
### elements
|item            |value                |
|----------------|---------------------|
|max step        |5000                 |
|output interval |250                  |
|x region        | [-10, 30]           |
|y region        | [-10, 10]           |
|dx              | 0.1                 |
|dy              | 0.1                 |
|Re              | 70                  |
|CFL             | 0.2                 |
|corr coeff      | 1.0                 |

### boundary conditions
|boundary type   |  conditions          |
|----------------|----------------------|
|inflow          |u=1, v=0, p=0         |
|outflow         |u,v: linear-extrapolation, p=0|
|prism surface   |u=v=0, dp/dn = 0      |


## Requirement
### main.f90
gfortran (gcc 8.2.0)

### pressure_output.py
python 3.6  
numpy  
matplotlib  
pandas  
imagemagick







