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
    def Bvec(temp):
        Br, By, Bz = [float(temp[N]) for N in range(3)]
        return np.array([Br, By, Bz])
    def vvec(temp):
        vr, vtheta, vz = temp
        return np.array([vr, vtheta, vz])
    def rvec(temp):
        rr, rtheta, rz = temp
        return np.array([rr, rtheta, rz])
    
# importation and interpolation of the Magnetic Fields
Brfield, Bzfield, Bfield = st.data_manager.filemanager('mfield2','load');
Brinterp = st.tools.interp_(Brfield);
Bzinterp = st.tools.interp_(Bzfield);

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
    mass = 1.673*10**-27;
    charge = sp.Beam.charge*1.6*10**-19;
    alpha = charge / (2*N*mass);
    dt = sp.Beam.dz/sp.Beam.Vz;

    # rotation vectors
    t = alpha*(vectors.Bvec((0,0,Bzinterp[n])))*dt;
    s = 2*t / (1+(np.vdot(t,t)));
    
    # transformation vectors
    vminus = vmat[n,:]
    vprime = vminus + np.cross(vminus,t)
    
    #change in velocity vector    
    vrotation = np.cross(vprime,s)
    
    return t, s, vrotation

def boris():
    #constants for reference
    sp.Beam.dz = 1*10**-3
    rvec,zvec = sc.tools.vectors();
    itersteps = np.size(zvec);
    itercounter = np.arange(np.size(zvec)-1);
    dt = sp.Beam.dz/sp.Beam.Vz;
    z0 = zvec[0]
    r0 = sp.Inner_coil.Rmin/2.0
    #storage matrices for various vectors
    vptest, vrottest = storagematrix(itersteps)
    ttest, stest = storagematrix(itersteps)
    
    #loading container lists
    velocitymatrix, positionmatrix = storagematrix(itersteps);
    
    #populate the 0th list position of each* container matrix
    velocitymatrix[0,:] = vectors.vvec((1.0, -8*10**-6/r0, sp.Beam.Vz));
    positionmatrix[0,:] = vectors.rvec((r0,0,z0));
    
    #half step forward in velocity
    positionmatrix[0,:] = positionmatrix[0,:]
    for ii in itercounter:
        zpos = zvec[ii]
        rpos = positionmatrix[ii,0]
        t, s, vrotation = fun(ii,rpos,zpos,vmat = velocitymatrix, pmat = positionmatrix)
        ttest[ii+1,:]  = t
        stest[ii+1,:] = s
        vrottest[ii+1,:] = vrotation
        velocitymatrix[ii+1,:] = (velocitymatrix[ii,:] + vrotation)
        positionmatrix[ii+1,:] = (positionmatrix[ii,:] + velocitymatrix[ii,:]*dt)
    return positionmatrix, velocitymatrix