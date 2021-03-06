#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 15:42:58 2020

@author: paulg

Program for simple visualisation of the piecewise 3D function, expressed as a combination of heavisides,
associated with the metric deformation generated by a cylinder of constant electromagnetic energy density
moving at the speed of light c.
"""

#import sys
#reload(sys)
#sys.setdefaultencoding('utf-8') #Uncomment if there is an encoding problem

from mpl_toolkits import mplot3d

from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

###Default constants###

L=20
R=5

"""
CONDITION FUNCTIONS
"""

def v1(L,R,z,t):
    return -t +L - 2*z - np.sqrt((2*z)**2 - 8*L*(z+t)) 

def v2(L,R,z,t):
    return -t +L - 2*z + np.sqrt((2*z)**2 - 8*L*(z+t)) 

def z1(L,R,z,t):
    return L - np.sqrt(L*(2*t+L))

def zR(L,R,z,t):
    return -0.5*z + 0.5*(R**2)/z

def zRL(L,R,z,t):
    return 0.5*(L-z)-0.5*(R**2)/(L-z)

def va(L,R,z,t):
    return - np.sqrt(t**2 - R**2) 

def vb(L,R,z,t):
    return np.sqrt(t**2-R**2)

"""
ANTIDERIVATIVES (Multiply return by 0 to remove afunction for test purposes)
"""

def prodsup(v,L,R,z,t): #v represents the integration variable
    return 0.25*v*(2*L - 2*z-v)

def ct(v,L,R,z,t):
    return 0.5*v*t*np.ones(np.shape(R))

def Rplusz(v,L,R,z,t):
    return 0.25*(v*np.sqrt(R**2 + v**2)+ R**2 * np.arcsinh(v/R))

def prodinf(v,L,R,z,t):
    return 0.25*v*(v + 2*z)*np.ones(np.shape(R))

def zer(v,L,R,z,t):
    return -0.25*np.sign(v)*v**2*np.ones(np.shape(R))

"""
COMPLETE SOLUTION
"""

def Solution(L,R,z,t):
    CondGen=[t*np.ones(np.shape(L))>=0,L-z>=0] # General conditions list. Solution is 0 if not verified
    BMaxGen=[t*np.ones(np.shape(L)),0.5*L*np.ones(np.shape(z)) - 0.5*z] # General upper bounds
    BminGen=[-t*np.ones(np.shape(L))] # General lower bounds
    
    FoncList=[prodsup,ct,ct,ct,Rplusz,Rplusz,Rplusz,prodinf,zer,zer] #Function list. Functions are repeated if its associated conditions are written as sum of Heaviside functions
    BMaxLoc=[[],[va(L,R,z,t),L-t-z],[L-t-z],[L-t-z],[vb(L,R,z,t),zRL(L,R,z,t)],[vb(L,R,z,t),zRL(L,R,z,t),-z*np.ones(np.shape(L)),zR(L,R,z,t)],[vb(L,R,z,t),zRL(L,R,z,t),-z*np.ones(np.shape(L))],[-0.5*z*np.ones(np.shape(L))],[np.zeros(np.shape(z))*np.ones(np.shape(L))],[]] #Associated upper bounds
    BminLoc=[[zRL(L,R,z,t),L-z-t],[(-t-z)*np.ones(np.shape(L))],[vb(L,R,z,t),(-t-z)*np.ones(np.shape(L))],[(-t-z)*np.ones(np.shape(L))],[va(L,R,z,t),-z*np.ones(np.shape(L))],[va(L,R,z,t)],[va(L,R,z,t),zR(L,R,z,t)],[zR(L,R,z,t) ,(-t-z)*np.ones(np.shape(L))],[],[np.zeros(np.shape(z))*np.ones(np.shape(L)),-0.5*z*np.ones(np.shape(L))]] #Associated lower bounds

    
    for k,i in enumerate(BMaxLoc): # General bounds are added to local bounds list
        i.extend(BMaxGen)
        BMaxLoc[k]=np.array(i) # Lists turned into array to apply numpy functions
    for l,j in enumerate(BminLoc):
        j.extend(BminGen)
        BminLoc[l]=np.array(j)
        
#-----  np.zeros(np.shape(t),dtype=bool) can be added in some local conditions to negate some parts of the integration (especially useful for studying the effect of a condition on a function)
    CondLoc=[[np.ones(np.shape(t), dtype=bool)*np.ones(np.shape(L), dtype=bool)],[t>=R],[t>=R],[t<R],[np.ones(np.shape(t), dtype=bool)*np.ones(np.shape(L), dtype=bool),t>=R],[z*np.ones(np.shape(L))>=0,t>=R],[z*np.ones(np.shape(L))<0,t>=R],[z*np.ones(np.shape(L))<0],[z*np.ones(np.shape(L))>=0],[np.ones(np.shape(t), dtype=bool)*np.ones(np.shape(L), dtype=bool)]]

    for i,j in enumerate(CondLoc):
        j.extend(CondGen)
        CondLoc[i]=np.array(j)
    

    
    Bornes=[]
    for k,l in enumerate(CondLoc):
        BMax=np.nanmin(BMaxLoc[k],axis=0) #Taking the upper bounds' minimum, excluding NaNs
        Bmin=np.nanmax(BminLoc[k],axis=0) #Taking the lower bounds' maximum, excluding NaNs
        

        CondFailed=np.logical_not(np.all(l,axis=0)) # Set as True where conditions are not met
        
        BMax[CondFailed]=0 # Bounds set to 0 so that the function gives 0 on these points
        Bmin[CondFailed]=0
        
        minGtMax=np.greater_equal(Bmin,BMax) # Returns True where lower bound is >= upper bound
        
        BMax[minGtMax]=0 # Bounds set to 0 so that the function gives 0 on these points
        Bmin[minGtMax]=0
        
        Bornes.append([BMax,Bmin])


    Sol=np.zeros(np.shape(z))*np.ones(np.shape(L)) #Initializing the solution
    for i,j in enumerate(FoncList):
        Sol+=j(Bornes[i][0],L,R,z,t)-j(Bornes[i][1],L,R,z,t) #Adds to the solution F(Up_bound)-F(Low-bound) for each function.

        
    return [Sol,Bornes] #Boundary list is also returned in order to easily debug the progrtam for small arrays
   
    
