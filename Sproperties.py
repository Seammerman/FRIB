# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 10:27:52 2016

solenoid configuration file

All information regarding: size, shape, turn, current, beam properties etc

the goal of this file will be to make the properties easily accessible through
simple reference of class or function
@author: samme
"""

## solenoid constants (length, width, currnet, etc..)
'''
using classes to store information regarding solenoid and beam properties should
allow for flexability and simple control over parameters without altering functional 
default properties.

**note for future: impliment a function that takes a configuration file
                    config file and function will update class imformation
'''

class Inner_coil():
    turns = 312
    layers = 20
    Rmax = (8.6/2 * 1/100)
    Rmin = (5.4/2 * 1/100)
    WT = (0.8 * 1/1000)
   
class Outer_coil():
    turns = 454
    layers = 24
    Rmax = (11.24/2 * 1/100)
    Rmin = (8.6/2 * 1/100)
    WT = (0.55 * 1/1000)
        
class solenoid():
    length = (25 * 1/100)
    Rmax = Outer_coil.Rmax
    Rmin = Inner_coil.Rmin
    current = 110
    max_B = 9

class Beam():
    T = (0.5 * 10**6)                   # energy per nucleon .5 MeV/u
    charge = 33                         # beam charge
    mass = (938 * 10**6)                # eV per nucleon based on nuetron/proton
    N = 238                             # number of nucleons
    c = 299792458                       # speed of light
    Cpc = (2 * mass)**.5 * (T)**.5      # classical momentum
    JP = 8 * 10 **-6                    # J/p momentum value
    Vz = Cpc / mass * c                 # Velocity in the z direction, p = mv ; p/m = v
    dz = 1/1000                         # 1 mm differential length
    
    # to be used for relativistic energies
    beta = Vz / c                       # relativity value beta = (v/c)
    gamma = 1 / (1 - (beta)**2)**(0.5)  # relativity value gamma