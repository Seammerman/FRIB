# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 10:34:59 2016

Solenoid "tool" reference file

The purpose of this file is to provide a discrete set of scripts to be used for
calculating region geometry as well as neccesary fields and interpolation objects

@author: samme
"""
import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import Sproperties as SP
# asign Sproperties objects to simple references
Beam = SP.Beam
solenoid = SP.solenoid
Inner_coil = SP.Inner_coil
Outer_coil = SP.Outer_coil

class tools():
    # this class will hold definitions for detector grids, interpolation, etc..
    
    def vectors():
        # spatial restrictions
        z_max = solenoid.length * (1 + 1)
        r_max = Inner_coil.Rmin
        
        zgap = (2 * z_max) / (Beam.dz)
        rgap = (2 * r_max) / (Beam.dz)
      
        # position vectors
        zvec = np.linspace(-z_max, z_max, zgap + 1)
        rvec = np.linspace(-r_max, r_max, rgap + 1)
        return rvec, zvec
          
    def grid():
        # grids
        rvec,zvec = tools.vectors()
        rgrid, zgrid = np.meshgrid(rvec,zvec)
        return rgrid, zgrid
        
    def interp_(field):
        from scipy import interpolate as interpolate
        r,z = tools.vectors()
        return interpolate.interp2d(r,z,field, kind='cubic')

    def calculator(Rmax,r,z):
        # initial values
        
        Br = np.zeros(np.shape(r))
        Bz = np.zeros(np.shape(z))
        r = np.abs(r)
        #constants        
        I = solenoid.current
        mew = 4 * np.pi * 10**(-7)
        C = mew * I / np.pi
        
        # special variables with (eliptic integrals)
        a = Rmax
        alphatwo = a**2 + r**2 + z**2 - 2 * a * r
        betatwo = a**2 + r**2 + z**2 + 2 * a * r        
        ksq = 1 - (alphatwo / betatwo)
        
        E = scipy.special.ellipe(ksq)
        K = scipy.special.ellipk(ksq)
        
        # fields
        Br = (C * z) / (2 * alphatwo * (betatwo**.5) * (r)) * ((a**2 + r**2 + z**2) * E - alphatwo * K)
        Bz = C / (2 * alphatwo * (betatwo**.5)) * ((a**2 - r**2 - z**2) * E + alphatwo * K)
        Bz_replacement = (mew * I * a**2) / (2*(a**2+z**2)**(3/2))
        
        # fixing
        np.place(Br, np.isnan(Br),0)
        np.place(Br, np.isinf(Br),0)
        np.place(Bz, np.isnan(Bz),0)
        np.place(Bz, np.isinf(Bz),0)
    
        np.place(Bz, r==0, Bz_replacement[np.where(r==0)])
        return Br, Bz
        
class data_manager():
    
    def plot3D(field, stride):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        r,z = tools.grid()
        
        Axes3D.plot_surface(ax,r,z,field, rstride = stride, cstride = stride)
        
        ax.set_xlabel('Radius')
        ax.set_ylabel('Z-Axis')
        ax.set_zlabel('field')
        ax.set_title('Surface Plot for some field')
        
        plt.show()
        
    def plot2D(x,y,*title):
        import matplotlib.patches as patches
        fig = plt.figure()
        plt.plot(x,y)
        plt.xlabel('z-axis')
        plt.ylabel('variable of interest')
        plt.title(title)
        '''
        rectangle is added if:
         - the z axis extends beyond the ends of the solenoid 
         - the y values are all positive
        '''
        if np.max(x) > 0.25 and 0 < np.all(y) :
            ax1 = fig.add_subplot(111)
            ax1.add_patch(patches.Rectangle(
                (-solenoid.length / 2 ,0),
                solenoid.length,
                Inner_coil.Rmin ,
                fill = False))
        plt.show()
        
    def plotpolar(theta,r):
        ax = plt.subplot(111, projection ='polar')
        ax.plot(theta,r)
        ax.set_rmax(0.027)
        ax.grid(True)
        ax.set_title('polar plot')
        plt.show()
        
    def filemanager(filename, option, *fieldlist):
        '''
        filename - specify the 'filename' by the text preceeding the .
        option - 'save' and 'load' including the ''
        *fieldlist - to be used when saving a set of fields. ie: (Br, Bz, B)
        '''
        import pickle
        
        if option == 'save':
            f = open(filename + '.pckl', 'wb')
            pickle.dump(fieldlist,f)
        elif option == 'load':
            f = open(filename + '.pckl', 'rb')
            objtuple = pickle.load(f)
            objlist = [x for t in objtuple for x in t]
            if np.size(objlist)>1: #conditional for any set of fields of size 3
                Br  = objlist[0]
                Bz = objlist[1]
                B = objlist[2]
                return Br, Bz, B
            else: #conditional for any set of fields of size 1
                K = objlist
                return K
        else:
            print('Invalid option')
            print('Choose load or save')
            option = input()
            tools.filemanager(filename,option, *fieldlist)
        

