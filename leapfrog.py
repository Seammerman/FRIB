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

#position, velocity, and acceleration vectors (if needed)
class vectors():
    def position():
        return [[r, theta, z]]
    def velocity():
        return [[vr, vtheta, vz]]
    def acceleration():
        return[[ar, atheta, az]]

class onetimecalc():
    '''
    the purpose of this function is to generate maps or values that only need to be done once
    -any maps or values that are calculated from this function will be saved for later use
    '''
    def __init__(self):
        rvector, zvector = st.tools.vectors()
        self.r = rvector
        self.z = zvector
        self.rgrid, self.zgrid = np.meshgrid(self.r,self.z)
    def kfield(self): # generate a map of K values
        Kmap = sc.calcs.K(self.rgrid,self.zgrid);
    # save the Kmap data to be loaded and used without redundent calculation
        st.data_manager.filemanager('kmap','save', Kmap)
    
#interpolation functions for b-field and K-field
class field():
    # import pre-calculated fields
    def __init__(self):
        self.Brmap, self.Bzmap, self.Bmap = st.data_manager.filemanager('mfield2','load')
    
        try: #attampt to load map data, if exception then generate map data
            self.Kmap = st.data_manager.filemanager('kmap','load')
        except FileNotFoundError:
            test = onetimecalc()
            test.kfield()
            print("Oops! needed to do some calculations first")
        
    def magnetic(self):
        return [[st.tools.interp_(self.Brmap),st.tools.interp_(self.Bzmap), st.tools.interp_(self.Bmap)]]
    def force():
        return st.tools.interp_(self.Kmap)