# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 14:20:34 2016

Boris method used for beam simulation

@author: samme
"""
import numpy as np
import Sproperties as sp
import Stools as st
import Scalculator as sc

class vectors():
    def Bvec(z):
        if np.abs(z) <= 0.125:
            mag_vec = np.array([0,0,9.0])
        else:
            mag_vec = np.array([0,0,0])
        return mag_vec
    def vvec(temp):
        vr, vtheta, vz = temp
        return np.array([vr, vtheta, vz])
    def rvec(temp):
        rr, rtheta, rz = temp
        return np.array([rr, rtheta, rz])
    

def storagematrix(Nsteps):
    velocitymatrix = np.zeros((Nsteps,3));
    positionmatrix = np.zeros((Nsteps,3));
    return velocitymatrix, positionmatrix
    
def vel(n,vrotation,**kwargs):
    '''
    need to redefine this function to work with sympy matricies
    '''
    vmat = kwargs.pop('vmat')
    velplus = vmat[n,:] + vrotation;
    
    vr = velplus[0];
    vtheta = velplus[1];
    vz = velplus[2];
    return (vr, vtheta, vz)
    
def pos(x,**kwargs):
    '''
    rvectors are stored in a 1-D list, position, this will set the first entry
        to the proper r0 vector
    the current list position, x, is evaluated based on r = r0 + vt + at^2/2
    x - index, current position in the interation list
    *argv - in principle should be a pass-through for velocity and a
    '''
    pmat = kwargs.pop('pmat')
    vmat = kwargs.pop('vmat')
    dt = sp.Beam.dz / sp.Beam.Vz;
    posplus = pmat[0,:] + vmat[x,:]*dt;
    rr = posplus[0];
    rtheta = posplus[1];
    rz = posplus[2];
    return (rr, rtheta, rz)
    
def fun(n,r,z,**kwargs):
    '''
    n - index for velocity vector
    r,z - positions for magnetic field vector
    '''
    vmat = kwargs.pop('vmat')
    pmat = kwargs.pop('pmat')
    N = sp.Beam.N;
    mass = N*1.673*10**-27;
    charge = sp.Beam.charge*1.6*10**-19;
    alpha = charge / (2*mass);
    dt = sp.Beam.dz/sp.Beam.Vz;

    # rotation vectors
    magneticfield = vectors.Bvec(pmat[n,2])
    t = alpha*(magneticfield)*dt;
    s = 2*t / (1+(np.dot(t,t)));
    
    # transformation vectors
    vminus = vmat[n,:]
    vprime = vminus + np.cross(vminus,t)
    
    #change in velocity vector    
    vrotation = np.cross(vprime,s)
    
    return vrotation

def boris():
    #constants for reference
    sp.Beam.dz = 1*10**-3
    rvec,zvec = sc.tools.vectors();
    itersteps = np.size(zvec);
    itercounter = np.arange(np.size(zvec)-1);
    dt = sp.Beam.dz/sp.Beam.Vz;
    z0 = zvec[0]
    r0 = sp.Inner_coil.Rmin/2.0
    theta0 = (8*10**(-6)*sp.Beam.Vz/r0**2)
    #loading container lists
    velocitymatrix, positionmatrix = storagematrix(itersteps);
    
    #populate the 0th list position of each* container matrix
    velocitymatrix[0,:] = vectors.vvec((0, -theta0*r0, sp.Beam.Vz));
    positionmatrix[0,:] = vectors.rvec((r0,0,z0));
    
    #half step forward in velocity
    positionmatrix[0,:] = positionmatrix[0,:]
    for ii in itercounter:
        zpos = zvec[ii]
        rpos = positionmatrix[ii,0]
        vrotation = fun(ii,rpos,zpos,vmat = velocitymatrix, pmat = positionmatrix)
        
        velocitymatrix[ii+1,:] = (velocitymatrix[ii,:] + vrotation)
        positionmatrix[ii+1,:] = (positionmatrix[ii,:] + velocitymatrix[ii,:]*dt)
    return positionmatrix, velocitymatrix
    
def BMmain():
    from matplotlib import pyplot as plt
    
    N = sp.Beam.N;
    mass = N*1.673*10**-27
    
    position, velocity = boris()
    radius = position[:,0]
    angle = position[:,1]
    z = position[:,2]
    
    momentum = mass*velocity
    angular_momentum = np.cross(position,momentum)
    
    plt.figure('r vs z')
    plt.plot(z,radius)
    plt.title('radius vs z')
    plt.ylabel('radius')
    
    
    plt.figure('theta vs angle')
    plt.plot(radius,angle)
    plt.title('angle vs radius')
    plt.ylabel('theta(rad)')
    
    
    plt.figure('angle vs z')
    plt.plot(z,angle)
    plt.title('angle vs z')
    plt.ylabel('theta')
    
    
    plt.figure('L vs z')
    plt.plot(z,np.abs((angular_momentum[0,2]-angular_momentum[:,2])/angular_momentum[0,2]))
    plt.title('angular momentum vs z')
    plt.ylabel('delta L / L')
    plt.show()

    pass