# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 21:41:15 2016
leap frong with yoshi's formula
@author: samme
"""
import numpy as np
import Sproperties as sp
import Stools as st
import Scalculator as sc


#magnetic field vector, piecewise inside of solenoid vs outside of solenoid
def Bvec(z):
    if np.abs(z) <= 0.125:
        magnetic_field_vec = np.array([0,0,9.0])
    else:
        magnetic_field_vec = np.array([0,0,0])
    return magnetic_field_vec
    
#equation of motion provided by Yoshi, R'' = -Kr
def K(r,z):
    sp.Beam.JP = 8*10**(-6)
    Q = sp.Beam.charge
    N = sp.Beam.N
    c = sp.Beam.c
    Cpc = sp.Beam.Cpc
    Bz = Bvec(z)[2]
    Ko = (1/4 * (Q/N * c / Cpc)**2 * Bz**2)
    return Ko - ((sp.Beam.JP) * (1/r**2))**2

def delTheta(r,z):
    sp.Beam.JP = 8*10**(-6)
    Q = sp.Beam.charge
    N = sp.Beam.N
    c = sp.Beam.c
    Cpc = sp.Beam.Cpc
    Bz = Bvec(z)[2]
    return -(1/2 * (Q/N * c / Cpc) * Bz) + ((sp.Beam.JP) * (1/r**2))

def storagematrix(Nsteps):
    velocitymatrix = np.zeros((Nsteps,3));
    positionmatrix = np.zeros((Nsteps,3));
    accelerationmatrix = np.zeros((Nsteps,3))
    return velocitymatrix, positionmatrix, accelerationmatrix

def velocitypush(r,z,ii, **kwargs):
    dz = st.Beam.dz
    velocity = kwargs.pop('vmat')
    Vminus = velocity[ii,:]
    Vplus = Vminus + np.array([-K(r,z)*r,0,0])*dz
    return Vplus, Vminus
    
def leapfrog():
    #constants
    rvec,zvec = st.tools.vectors()
    Vz = sp.Beam.Vz;
    dz = sp.Beam.dz
    dt = dz/Vz; #dz/dt = Vz, dt/dz = 1/Vz
    z0 = zvec[0]
    r0 = sp.Inner_coil.Rmin / 2.0
    
    
    #iteration information and initialization of matrices
    itersteps = np.size(zvec)
    itercounter = np.arange(np.size(zvec)-1)
    velocity, position, acceleration = storagematrix(itersteps)
    
    #establishing initial conditions
    position[0,:] = np.array([r0,0,z0])
    velocity[0,:] = np.array([0,-delTheta(r0,z0)*r0,1]) #the initial velocity Vz here is 1 because we are in terms of z, not t
    acceleration[0,:] = np.array([-K(r0,z0)*r0,0,0])
    
    #half-step backward for velocity to establish the leapfrog scheme
    velocity[0,:] = velocity[0,:] - acceleration[0,:]*dz/2

    for ii in itercounter:
        zpos = zvec[ii]
        rpos = position[ii,0]
        
        Vplus, Vminus = velocitypush(rpos,zpos,ii,vmat = velocity)
        Vavg = (Vplus+Vminus)/2
        thetadot = np.array([0,delTheta(rpos,zpos),0])
        position[ii+1,:] = position[ii,:] + (Vavg + thetadot)*dz
        velocity[ii+1,:] = Vplus + thetadot
        acceleration[ii,0] = -K(rpos,zpos)*rpos

    return position, velocity, acceleration
    
def LFmain():
    
    from matplotlib import pyplot as plt
    
    N = sp.Beam.N;
    mass = N*1.673*10**-27
    
    position, velocity, acceleration = leapfrog()
    
    radius = position[:,0]
    angle = position[:,1]
    z = position[:,2]
    
    momentum = mass*velocity
    angular_momentum = np.cross(position,momentum)
    
    plt.figure('r vs z')
    plt.plot(z,radius)
    plt.title('radius vs z')
    plt.ylabel('radius')
    plt.show()
    
    plt.figure('radius vs angle')
    #plt.subplot('111',projection='polar')
    plt.plot(radius,angle)
    plt.title('angle vs radius')
    plt.ylabel('theta(rad)')
    plt.show()
    
    plt.figure('angle vs z')
    plt.plot(z,angle)
    plt.title('angle vs z')
    plt.ylabel('theta')
    plt.show()
    
    plt.figure('L vs z')
    plt.plot(z,np.abs((angular_momentum[0,2]-angular_momentum[:,2])/angular_momentum[0,2]))
    plt.title('angular momentum vs z')
    plt.ylabel('delta L / L')
    plt.show()
    
    pass
    
