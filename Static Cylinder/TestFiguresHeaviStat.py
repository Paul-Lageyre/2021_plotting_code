#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 17:52:44 2020

@author: paulg
"""

from VisuHeaviStat import *
from cycler import cycler

########################      plot parameters     ########################



default_cycler = (cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#17becf'])+
cycler(linestyle=['-',(0,(3,1)),':',(0,(4,1,1,1)),(0,(3,5,3,5,1,5)),(0,(3,5,1,5,1,5)),(0,(30,5,1,5)),(0,(10,5,3,5))]))

plt.rcParams['figure.figsize']=[18.53,9.55] # used to set the figures to fullsize (important for tight_layout)
plt.rcParams['font.size']=35
plt.rc('lines',lw=8)
plt.rc('axes', prop_cycle=default_cycler)
plt.rc('xtick.major', size=14.0, width=3.2)
plt.rc('xtick.minor', size=8.0, width=2.4)
plt.rc('ytick.major', size=14.0, width=3.2)
plt.rc('ytick.minor', size=8.0, width=2.4)
#############################################################################

c= 299792458.

Khi= 8*np.pi *6.06742e-11 / c**4 # Einstein's constant/2

l0= 1e-6 # Problem's scaling length

I=1e26  #Laser beam intensity in W/m²
V=np.pi*500 #Scaled up volume of the cylinder of height L and radius R

E = l0**3*V*I/c #Energy held in the cylinder

Domaine=[-200,200,0,100] # plotting domain: [zmin,zmax,c*tmin,c*tmax]
Npoints = 500 # Number of points per variable: Npoints² total points for 2 variables.

x = np.linspace(Domaine[0],Domaine[1],Npoints)
y = np.linspace(Domaine[2],Domaine[3],Npoints)

X, Y = np.meshgrid(x,y)


def tslice(t,N=4,lim=None, ratiomin=0.004, ratiomax=100):
    """
    Cross-section at time t of the solution(z,ct) centered in L/2 for various size ratios
    
    t: c*time of plot
    N: number of different size ratios
    lim: puts the same specified upper limit for all graphs, while by default each graph is scaled on its own
    ratiomin: lowest size ratio
    ratiomax: highest size ratio
    """
    fig = plt.figure()
    ax=[]

    ratios = np.logspace(np.log10(ratiomin),np.log10(ratiomax),N)
    for a,i in enumerate(ratios):
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))

        Z= Solution(L,R,x+L/2,t*np.ones(np.shape(x)))
#        x=x/L
        line, = plt.plot(x, Khi*I*(l0**2)*Z/c)
        if lim:
            line.axes.set_ylim(0,lim)
        line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        
        label=''
        if L<0.5: 
            if L<0.003:
                label+=u'L=%.1e'
            else:
                label+=u'L=%.3f'
        else:
            label+=u'L=%.1f'
        label+=u'µm, '
        
        if R<0.5: 
            if R<0.003:
                label+=u'R=%.1e'
            else:
                label+=u'R=%.3f'
        else:
            label+=u'R=%.1f'
        label+=u'µm'
        line.set_label(label%(L,R))
        ax.append(line)
    plt.xlabel(u'z- L/2 (µm)')
    plt.ylabel('h')

    plt.legend()
    plt.tight_layout(pad=1,rect=(0,0,1,.95))
    


def xslice(tmax,x0=0,N=4,lim=None, centered=True, ratiomin=0.004, ratiomax=100):
    """
    Cross-section at distance x0 from z=L/2 of the solution(z,t) for various size ratios
    
    tmax: maximum value of c*time in plot
    x0: position on axis z from z=L/2
    N: number of different size ratios
    lim: puts the same specified upper limit for all graphs, while by default each graph is scaled on its own
    centered: graphs are given with z=L/2 as 0. If False, 0 is at z=0
    ratiomin: lowest size ratio
    ratiomax: highest size ratio
    """
    fig = plt.figure()
    ax=[]

    ratios = np.logspace(np.log10(ratiomin),np.log10(ratiomax),N)
    yloc = np.logspace(np.log10(0.3),np.log10(tmax),Npoints)
    for a,i in enumerate(ratios):
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))
        if centered:
            Z= Solution(L,R,L/2+x0*np.ones(np.shape(yloc)),yloc)
        else:
            Z= Solution(L,R,x0*np.ones(np.shape(yloc)),yloc)
        line, = plt.semilogx(yloc, Khi*I*(l0**2)*Z/c)
        if lim:
            line.axes.set_ylim(0,lim)
        line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-4,4))
        
        label=''
        if L<0.5: 
            if L<0.003:
                label+=u'L=%.1e'
            else:
                label+=u'L=%.3f'
        else:
            label+=u'L=%.1f'
        label+=u'µm, '
        
        if R<0.5: 
            if R<0.003:
                label+=u'R=%.1e'
            else:
                label+=u'R=%.3f'
        else:
            label+=u'R=%.1f'
        label+=u'µm'
        line.set_label(label%(L,R))
        
        ax.append(line)
    
    plt.xlabel(u"ct (µm)")
    plt.ylabel('h')
    plt.legend()
    plt.tight_layout(pad=1,rect=(0,0,1,.95))
    


def maxStat(rmin,rmax,N): 
    """
    Maximum of solution(z,ct) for a range of aspect ratios
    
    rmin: minimum power of 10 for the aspect ratio
    rmax: maximum power of 10 for the apect ratio
    N: number of points
    """
    rs=np.logspace(rmin,rmax,N)
    maxlist=[]
    for i in rs:
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))
        Z=Solution(L,R,[L/2],[10000])[0]
        maxlist.append(Khi*I*(l0**2)*Z/c)
    fig=plt.figure(figsize=(16,9))
    line,=plt.semilogx(rs,maxlist)
    line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    plt.xlabel("Aspect ratio L/R")
    plt.ylabel('Maximum of h')
    plt.tight_layout(pad=1,rect=(0,0,1,.95))

    
def tFWHM(t,rmin,rmax,N): 
    """
    Gives the full width at half maximum of the gravitational perturbation as a function of the aspect ratio for a given time and 
    
    t: c*time of evaluation
    rmin: minimum power of 10 for the aspect ratio
    rmax: maximum power of 10 for the apect ratio
    N: number of points
    """
    rs=np.logspace(rmin,rmax,N)
    fmax=[]
    for i in rs:
        L=np.cbrt(V*i**2 /np.pi)
        R=np.cbrt(V/(i*np.pi))
        absZ= np.abs(np.array(Solution(L,R,x+L/2,t*np.ones(np.shape(y))))-np.float(Solution(L,R,[L/2],[t])[0])/2)
        fmax.append(np.abs(x[absZ.argmin()]))
    fig=plt.figure(figsize=(16,9))
    line,= plt.semilogx(rs,l0*np.array(fmax))
    line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    plt.xlabel("Aspect ratio L/R")
    plt.ylabel('FWHM (m)')
    plt.tight_layout(pad=1,rect=(0,0,1,.95))

    
def surfplot(L,R,N,domain=Domaine):
    """
    Gives the 3D and 2D colorimetric plot of the deformation on a preset domain.
    
    L: cylinder length
    R: cylinder base radius
    N: number of points for each axis, total point number is N²
    domain: domain on (z,ct) upon which the solution is ploted
    """
    if np.shape(domain) != (4,):
        raise TypeError('domain must be a (4,) shape callable list, array or tuple')
    else: 
        fig = plt.figure(figsize=(16,9))
        
        fig2, ax2 = plt.subplots(figsize=(16,9))
        
        ax2.ticklabel_format(axis='both', style='sci', scilimits=(-4,4))
        ax2.set_xlabel(u"z (µm)")
        ax2.set_ylabel(u'ct (µm)')
        
        ax = fig.gca(projection='3d')
        ax.ticklabel_format(axis='both', style='sci', scilimits=(-4,4))
        ax.set_xlabel(u"z (µm)",labelpad=40)
        ax.set_ylabel(u'ct (µm)',labelpad=40)
        ax.set_zlabel('h',labelpad=20)
        ax.tick_params(axis='z', pad=10)
        
        xs = np.linspace(domain[0],domain[1],N)
        ys = np.linspace(domain[2],domain[3],N)
        
        Xs, Ys = np.meshgrid(xs,ys)
        
        Z= Solution(L,R,Xs,Ys)
        Z= Khi*(l0**2)*I*Z/c
        
        
        surf1=ax.plot_surface(Xs,Ys,Z , cmap= cm.RdBu) #3D Plot of the solution
        surfp=ax2.contourf(Xs,Ys,Z,100,extend='both') #Colormap plot of the solution
        
        cbar=fig.colorbar(surf1)
        cbar2=fig2.colorbar(surfp) 
        
        cbar.ax.set_ylabel('h', labelpad=40, rotation=270)
        cbar2.ax.set_ylabel('h', labelpad=40, rotation=270)
        
        fig.tight_layout(pad=1,rect=(0,0,1,.95))
        fig2.tight_layout(pad=1,rect=(0,0,1,.95))

