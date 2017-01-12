# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 11:11:06 2016

Solenoid calculation's file

the purpose of this file is to store functions and classes related to constructing
fields: Br, Bz, B, K, theta, etc

@author: samme
"""
import numpy as np
import scipy
import Sproperties as SP
import Stools as ST
# asign members of SP and ST to relavent names
    # SP objects
Beam = SP.Beam
solenoid = SP.solenoid
Inner_coil = SP.Inner_coil
Outer_coil = SP.Outer_coil
    #ST objects
tools = ST.tools

class calcs():
    # this class is responsible for calculatiing stuff (Br, Bz, B, K, etc..)

    def magnetic(coil):
        r,z = tools.grid()
        # selecting properties
        if coil == 'Inner':
            layers = Inner_coil.layers
            turns = Inner_coil.turns
            WT = Inner_coil.WT
            Rstart = Inner_coil.Rmin
            
        elif coil == 'Outer':
            layers = Outer_coil.layers
            turns = Outer_coil.turns
            WT = Outer_coil.WT
            Rstart = Outer_coil.Rmin
            
        else:
            print('which coil?')
            print('Inner or Outer')
            coil = input()
            calcs.magnetic(coil)
        # initializing zero matrices of correct size
        Br = np.zeros(np.shape(r))
        Bz = np.zeros(np.shape(z))
        
        print(coil + ' coil')
        print('')
        
        tstart = time.time() # keeping track of working time
        
        for ii in range(layers):
            Rmax = Rstart + WT * ii
            for jj in range(turns): 
                z_temp = z + (solenoid.length / 2 - WT * jj  * (solenoid.length / 0.25))
                Br_temp, Bz_temp = tools.calculator(Rmax,r,z_temp)
                Br = Br + Br_temp
                Bz = Bz + Bz_temp

        tstop = time.time()
        telapsed = tstop - tstart # change in time (final - initial)
        print('time elapsed = ' + str(np.round(telapsed, 1)))

        return Br, Bz
        
    def K(r,z):
        Q = Beam.charge
        N = Beam.N
        c = Beam.c
        Cpc = Beam.Cpc
        Bzfun = tools.interp_(Bz)
        Ko = (1/4 * (Q/N * c / Cpc)**2 * Bzfun(r,z)**2)
        return Ko - ((Beam.JP) * (1/r**2))**2
        
    def R_SODE():
        zmax = solenoid.length * 2 
        zgap = (2 * zmax) / (1/1000)
        R0 = Inner_coil.Rmin / 2
        dR0 = 1 * 1/1000
        vec0 = ([R0,dR0])
        zvec = np.linspace(-zmax,zmax, zgap + 1)
        ii = 0
        def dvec(vec,zvec):
            R,dR = vec
            d2Rdz2 = -calcs.K(R,zvec) * R
            return ([dR, d2Rdz2])
            
        R = scipy.integrate.odeint(dvec,vec0,zvec)
        
        return R,zvec
        
    def dthetadz(r,z):
        Bzfun = tools.interp_(Bz)
        return -1/2 * (Beam.charge/Beam.N) * (Beam.c/Beam.Cpc) * Bzfun(r,z) + (Beam.JP) * (1/r**2)
        
    def theta_ODE():
        # constants
        Q = Beam.charge
        N = Beam.N
        c = Beam.c
        Cpc = Beam.Cpc
        step_size = 1 * 1/1000
        # interpolation function for Bz
        Bzfun = tools.interp_(Bz)
        # grabbing R and z values from our R(z) solution
        R, z = calcs.R_SODE()
        # R has two columns of values with 0 having position values and 1 having velocity values
        r = R[:,0]
        theta = [0]
        for i in range(np.size(r)):
            theta = np.append(theta, theta[i] + calcs.dthetadz(r[i],z[i]) * step_size)
            
        theta = np.delete(theta, np.size(theta)-1)
        
        return theta, r, z
        
    def p_conservation(r,z):
        
        # calculate dtheta / dz 
        dth = calcs.dthetadz(r,z)
        # generate a function for the magnetic field, takes r and z as positional values
        Bzfun = tools.interp_(Bz)
        # formula constant. generalization of mew = integral(r Bz dr dtheta); no theta dependence -> dtheta = 2pi; r dependance of Bz << z - > mew = integral(r dr) = r**2/2
        def mew(r,z):
            return (2 * np.pi) * r**2 / 2 * Bzfun(r,z)
        # angular momentum, should be conserved
        L = Beam.gamma * Beam.mass * r**2 * dth + Beam.charge * mew(r,z) / (2 * np.pi)
        return L
