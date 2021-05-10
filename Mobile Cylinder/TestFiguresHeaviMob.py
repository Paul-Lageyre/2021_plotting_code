#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Created on Wed Sep 23 17:52:44 2020

@author: paulg
"""

from VisuHeaviMob import *
from cycler import cycler
import matplotlib.colors as colors
import matplotlib.ticker as ticker

########################     plot parameters     ########################



default_cycler = (cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#17becf'])+
cycler(linestyle=['-',(0,(3,1)),':',(0,(4,1,1,1)),(0,(3,5,3,5,1,5)),(0,(3,5,1,5,1,5)),(0,(30,5,1,5)),(0,(10,5,3,5))]))

plt.rcParams['figure.figsize']=[16,9] # used to set the figures to fullsize (important for tight_layout)
plt.rcParams['font.size']=35
plt.rcParams['grid.linewidth']=1.0
plt.rc('lines',lw=8)
plt.rc('axes', prop_cycle=default_cycler)
plt.rc('xtick.major', size=14.0, width=3.2)
plt.rc('xtick.minor', size=8.0, width=2.4)
plt.rc('ytick.major', size=14.0, width=3.2)
plt.rc('ytick.minor', size=8.0, width=2.4)
#############################################################################


c= 299792458.

Khi= 4*np.pi *6.06742e-11 / c**4 # Einstein's constant/2

l0= 1e-6 # Problem's scaling length

I=1e26 #Laser beam intensity in W/m²
V=np.pi*500 #Scaled up volume of the cylinder of height L and radius R

E = l0**3*V*I/c #Energy held in the cylinder

Domaine=[-150,100,0,100] # plotting domain: [zmin,zmax,c*tmin,c*tmax]
Npoints = 100 # Number of points per variable: Npoints² total points for 2 variables.

x = np.linspace(Domaine[0],Domaine[1],Npoints)
y = np.linspace(Domaine[2],Domaine[3],Npoints)

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
        
        if -2*L> Domaine[0] and L/2<Domaine[1]:
            x0=np.linspace(-2*L,L/2,100)
            x=np.concatenate((np.linspace(Domaine[0],-2*L,50)[:-2],x0,np.linspace(L/2,Domaine[1],50)[1:]))
        else:
            x=np.linspace(Domaine[0],Domaine[1],200)
#        
#        x=np.linspace(-10*L,L/2,100) # other plot setting
        ax.append(fig.add_subplot(int(np.ceil(np.sqrt(N))),int(np.ceil(np.sqrt(N))),a+1))
        Z= Solution(L,R,x+L/2,t*np.ones(np.shape(x)))[0]
#        x=x/L #for normalization purposes
        ax[a].plot(x, Khi*I*(l0**2)*Z/c)
        if lim:
            ax[a].set_ylim(0,lim)
        ax[a].ticklabel_format(axis='both', style='sci', scilimits=(-4,4))
        
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
        ax[a].set_title(label%(L,R), fontsize=plt.rcParams['font.size'])
        
        ax[a].set_xlabel(u'Z- L/2 (µm)')
        ax[a].set_ylabel('h')
    
    plt.tight_layout(pad=1,rect=(0,0,1,.95))


def xslice(tmax,x0=0,N=4,lim=None,centered=True, ratiomin=0.004, ratiomax=100):
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
    yloc = np.logspace(np.log10(0.1),np.log10(tmax),Npoints)
    for a,i in enumerate(ratios):
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))
        if centered:
            Z= Solution(L,R,L/2+x0*np.ones(np.shape(yloc)),yloc)[0]
        else:
            Z= Solution(L,R,x0*np.ones(np.shape(yloc)),yloc)[0]
        line, = plt.semilogx(yloc, Khi*I*(l0**2)*Z/c)
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
    
    plt.xlabel(u"ct (µm)")
    plt.ylabel('h')
    plt.legend()
    plt.tight_layout(pad=1,rect=(0,0,1,.95))
    
    
def xsliceP(tmax,x0=0,N=4,lim=None,centered=True, ratiomin=1e-5, ratiomax=0.1, P=1e15):
    """
    Cross-section at distance x0 from z=L/2 and at set power P of the solution(z,t) for various size ratios
    
    tmax: maximum value of c*time in plot
    x0: position on axis z from z=L/2
    N: number of different size ratios
    lim: puts the same specified upper limit for all graphs, while by default each graph is scaled on its own
    centered: graphs are given with z=L/2 as 0. If False, 0 is at z=0
    ratiomin: lowest size ratio
    ratiomax: highest size ratio
    P: set power taken
    """
    fig = plt.figure()
    ax=[]

    ratios = np.logspace(np.log10(ratiomin),np.log10(ratiomax),N)
    yloc = np.logspace(np.log10(10),np.log10(tmax),Npoints)
    
    ccycle=iter(plt.rcParams['axes.prop_cycle'])
    for a,i in enumerate(ratios):
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))
        if centered:
            Z= P*Solution(L,R,L/2+x0*np.ones(np.shape(yloc)),yloc)[0] /R**2
        else:
            Z= P*Solution(L,R,x0*np.ones(np.shape(yloc)),yloc)[0]/R**2
        line, = plt.semilogx(yloc, Khi*Z/c)
        stepcycle=next(ccycle)
        vcolor=stepcycle['color']
        vstyle=stepcycle['linestyle']
        line.axes.axvline(R**2/(2*L),c=vcolor,ls=vstyle)
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

    plt.xlabel(u"ct (µm)")
    plt.ylabel('h')
    plt.legend()
    plt.tight_layout(pad=1,rect=(0,0,1,.95))


def maxMob(t,rmin,rmax,N):
    """
    Maximum of solution(z,ct) for a range of aspect ratios at set time t
    A figure needs to be set beforehand for this to plot
    
    t: c*time of evaluation
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
        x=np.linspace(-2*L,L/2,100)
        Z=Solution(L,R,np.array([0]),np.array([t]))[0]
        maxlist.append(Khi*I*(l0**2)*Z/c)
    line, = plt.semilogx(rs,maxlist)
    line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    
    t_label=u'ct='
    if t>9999:
        t_label+=u'{:.1e}'
    else:
        t_label+=u'{:.0f}'
    t_label+=u'µm'
    line.set_label(t_label.format(t))
    
    plt.xlabel("Rapport d'aspect L/R")
    plt.ylabel('Maximum de h')
    plt.legend()


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
        fig = plt.figure()
        
        fig2, ax2 = plt.subplots()

        ax2.ticklabel_format(axis='both', style='sci', scilimits=(-4,4))
        ax2.set_xlabel(u'Z- L (µm)')
        ax2.set_ylabel(u'ct (µm)')

        ax = fig.gca(projection='3d')
        ax.ticklabel_format(axis='both', style='sci', scilimits=(-4,4))
        ax.set_xlabel(u'Z (µm)',labelpad=20)
        ax.set_ylabel(u'ct (µm)',labelpad=20)
        ax.set_zlabel('h',labelpad=20)
        ax.tick_params(axis='z', pad=10)

        xs = np.linspace(domain[0],domain[1],N)
        ys = np.linspace(domain[2],domain[3],N)

        Xs, Ys = np.meshgrid(xs,ys)

        Z= Solution(L,R,Xs,Ys)[0]


        surf1=ax.plot_surface(Xs,Ys,Khi*I*(l0**2)*Z/c , cmap= cm.RdBu) #3D Plot of the solution
        surfp=ax2.contourf(Xs,Ys,Khi*I*(l0**2)*Z/c,100)                #Colormap plot of the solution

        cbar=fig.colorbar(surf1)
        cbar2=fig2.colorbar(surfp)

        cbar.ax.set_ylabel('h', labelpad=20, rotation=270)
        cbar2.ax.set_ylabel('h', labelpad=20, rotation=270)


def maxPcte(t,N,Rmin=1,Rmax=1e5,Lmin=5,Lmax=5e3,NL=4,P=1e15,lim=None):
    """
    Maximum of deformation for a set time and a set power as a function of the cylinder radius R, for various cylinder lengths L
    t: c*time of evaluation
    N: number of points for R
    Rmin: minimum value of R
    Rmax: maximum value of R
    Lmin: minimum value of L
    Lmax: maximum value of L
    NL: number of different values of L
    P: set power taken
    lim: puts the same specified upper limit for all graphs, while by default each graph is scaled on its own
    """
    Rs=np.logspace(np.log10(Rmin),np.log10(Rmax),N)
    Ls=np.logspace(np.log10(Lmin),np.log10(Lmax),NL)

    fig, ax = plt.subplots(nrows=int(np.ceil(np.sqrt(NL))), ncols=int(np.ceil(np.sqrt(NL))), sharex=True, sharey=True)
    for i in range(len(ax)): ax[i,0].set_ylabel('h')
    for i in range(len(ax[-1])): ax[-1,i].set_xlabel(u'R (µm)')
    ax=ax.flatten()
    for a,i in enumerate(Ls):

        Z=np.zeros(N)
        for b,j in enumerate(Rs):
            Iloc=P/(np.pi*j**2)
            Z[b]= Iloc*Solution(i,j,np.array([0]),np.array([t]))[0]
            
        ax[a].loglog(Rs, Khi*Z/c)
        if i<t:
            ax[a].axvline(np.sqrt(2*i*t-i**2),ls='--',lw=plt.rcParams['lines.linewidth']*3/5)
        else:
            ax[a].axvline(i,lw=3,ls='--')
        if lim:
            ax[a].set_ylim(1e-4*lim,lim)
            
        label=''
        if i<0.5: 
            if i<0.003:
                label+=u'L=%.1e'
            else:
                label+=u'L=%.3f'
        elif i>9999:
            label+=u'L=%.1e'
        else:
            label+=u'L=%.1f'
        label+=u'µm'
        
        ax[a].set_title(label%(i))
        ax[a].yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=4))
#        ax[a].set_xlabel(u'R (µm)')
#        ax[a].set_ylabel('h')
        
    plt.tight_layout(pad=1,rect=(0,0,1,.95))


def maxPcteRt(N,Rmin=1,Rmax=1e5,Lmin=5,Lmax=5e3,NL=4,P=1e15,lim=None):
    """
    Maximum of deformation for a set power as a function of the cylinder radius R, for various cylinder lengths L.
    The c*time at which this maximum is evaluated is equal to double the Rayleigh length
    t: c*time of evaluation
    N: number of points for R
    Rmin: minimum value of R
    Rmax: maximum value of R
    Lmin: minimum value of L
    Lmax: maximum value of L
    NL: number of different values of L
    P: set power taken
    lim: puts the same specified upper limit for all graphs, while by default each graph is scaled on its own
    """
    Rs=np.logspace(np.log10(Rmin),np.log10(Rmax),N)
    Ls=np.logspace(np.log10(Lmin),np.log10(Lmax),NL)

    fig = plt.figure(figsize=(16,9)) #figsize=(16,9) pour avoir le plein écran pour l'application de tight_layout
    ax=[]
    for a,i in enumerate(Ls):

        ax.append(fig.add_subplot(int(np.ceil(np.sqrt(NL))),int(np.ceil(np.sqrt(NL))),a+1))

        Z=np.zeros(N)
        for b,j in enumerate(Rs):
            Iloc=P/(np.pi*j**2)
            t=np.pi*j**2
            Z[b]= Iloc*Solution(i,j,np.array([0]),np.array([t]))[0]

        ax[a].semilogx(Rs, Khi*Z/c)
        if lim:
            ax[a].set_ylim(1e-4*lim,lim)
        ax[a].ticklabel_format(axis='y', style='sci', scilimits=(-4,4))
        
        label=''
        if i<0.5: 
            if i<0.003:
                label+=u'L=%.1e'
            else:
                label+=u'L=%.3f'
        elif i>9999:
            label+=u'L=%.1e'
        else:
            label+=u'L=%.1f'
        label+=u'µm'
        
        ax[a].set_title(label%(i))
        
        ax[a].set_xlabel(u'R (µm)')
        ax[a].set_ylabel('h')
    plt.tight_layout(pad=1,rect=(0,0,1,.95))
