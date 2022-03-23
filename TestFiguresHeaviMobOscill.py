#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Created on Wed Sep 23 17:52:44 2020

@author: paulg
"""

from VisuHeaviMobOscill import *
from cycler import cycler
import matplotlib.colors as colors
import matplotlib.ticker as ticker

########################     plot parameters     ########################



default_cycler = (cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#17becf'])+
cycler(linestyle=['-',(0,(3,1)),':',(0,(4,1,1,1)),(0,(3,5,3,5,1,5)),(0,(3,5,1,5,1,5)),(0,(30,5,1,5)),(0,(10,5,3,5))]))

plt.rcParams['figure.figsize']=[18.53,9.55] # used to set the figures to fullsize (important for tight_layout)
plt.rcParams['font.size']=35
plt.rcParams['grid.linewidth']=1.0
plt.rc('lines',lw=8)
plt.rc('axes', prop_cycle=default_cycler, titlesize=u'medium')
plt.rc('xtick.major', size=14.0, width=3.2)
plt.rc('xtick.minor', size=8.0, width=2.4)
plt.rc('ytick.major', size=14.0, width=3.2)
plt.rc('ytick.minor', size=8.0, width=2.4)
#############################################################################


c= 299792458.

Khi= 8*np.pi *6.06742e-11 / c**4 # Einstein's constant/2

l0= 1e-6 # Problem's scaling length

wavelength=5

k0=2.*np.pi/wavelength

I=1e26 #Laser beam intensity in W/m²
V=np.pi*500 #Scaled up volume of the cylinder of height L and radius R

E = l0**3*V*I/c #Energy held in the cylinder

Domaine=[-100,100,0,100] # plotting domain: [zmin,zmax,c*tmin,c*tmax]
Npoints = 100 # Number of points per variable: Npoints² total points for 2 variables.

x = np.linspace(Domaine[0],Domaine[1],Npoints)
y = np.linspace(Domaine[2],Domaine[3],Npoints)

def tslice(t,N=4,lim=None, lim_inf=None, ratiomin=0.004, ratiomax=100, Npoints=5000):
    """
    Cross-section at time t of the solution(z,ct) centered in L/2 for various size ratios

    t: c*time of plot
    N: number of different size ratios
    lim: puts the same specified upper limit for all graphs, while by default each graph is scaled on its own
    ratiomin: lowest size ratio
    ratiomax: highest size ratio
    Npoints: number of points for the graph
    """
    if lim: fig, ax = plt.subplots(nrows=int(np.ceil(np.sqrt(N))), ncols=int(np.ceil(np.sqrt(N))), sharex=True, sharey=True)
    else: fig, ax = plt.subplots(nrows=int(np.ceil(np.sqrt(N))), ncols=int(np.ceil(np.sqrt(N))), sharex=True)

    for i in range(len(ax)): ax[i,0].set_ylabel('h')
    for i in range(len(ax[-1])): ax[-1,i].set_xlabel(u'Z- L/2 (µm)')
    ax=ax.flatten()

    ratios = np.logspace(np.log10(ratiomin),np.log10(ratiomax),N)
    for a,i in enumerate(ratios):
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))

        if -2*L> Domaine[0] and L/2<Domaine[1]:
            x0=np.linspace(-2*L,L/2,Npoints/2)
            x=np.concatenate((np.linspace(Domaine[0],-2*L,Npoints/4)[:-2],x0,np.linspace(L/2,Domaine[1],Npoints/4)[1:]))
        else:
            x=np.linspace(Domaine[0],Domaine[1],Npoints)
#
#        x=np.linspace(-10*L,L/2,100) # other plot setting
        Z= Solution(L,R,x+L/2,t*np.ones(np.shape(x)),k0)[0]
#        x=x/L #for normalization purposes
        ax[a].plot(x, Khi*I*(l0**2)*Z/c)
        if lim:
            if type(lim_inf)==int or type(lim_inf)==float :
                ax[a].set_ylim(lim_inf,lim)
            else:
                ax[a].set_ylim(-lim,lim)
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


#        ax[a].set_xlabel(u'Z- L/2 (µm)')
#        ax[a].set_ylabel('h')

    plt.tight_layout(pad=1,rect=(0,0,1,.95))


def tslice_kprop(t, Nk, N=4,lim=None, ratiomin=0.004, ratiomax=5000):
    """
    Cross-section at time t of the solution(z,ct) centered in L/2 for various size ratios, with number of oscillations
    staying constant to cylinder length

    t: c*time of plot
    Nk: number of periods in the light cylinder (ex: wavelength=5 for L=20--> Nk=4)
    N: number of different size ratios
    lim: puts the same specified upper limit for all graphs, while by default each graph is scaled on its own
    ratiomin: lowest size ratio
    ratiomax: highest size ratio
    """
    if lim: fig, ax = plt.subplots(nrows=int(np.ceil(np.sqrt(N))), ncols=int(np.ceil(np.sqrt(N))), sharey=True)
    else: fig, ax = plt.subplots(nrows=int(np.ceil(np.sqrt(N))), ncols=int(np.ceil(np.sqrt(N))))

    for i in range(len(ax)): ax[i,0].set_ylabel('h')
    for i in range(len(ax[-1])): ax[-1,i].set_xlabel(u'Z- L/2 (µm)')
    ax=ax.flatten()

    ratios = np.logspace(np.log10(ratiomin),np.log10(ratiomax),N)
    for a,i in enumerate(ratios):
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))

        k0=2.*np.pi*Nk/L

        x=np.linspace(-2*L,L/1.8,1000)
#        Amplitude Max
#        x=np.linspace(-10*L,L/2,100) # other plot setting
        Z= Solution(L,R,x+L/2,t*np.ones(np.shape(x)),k0)[0]
#        x=x/L #for normalization purposes
        ax[a].plot(x, Khi*I*(l0**2)*Z/c)
        if lim:
            ax[a].set_ylim(-lim,lim)
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


#        ax[a].set_xlabel(u'Z- L/2 (µm)')
#        ax[a].set_ylabel('h')

    plt.tight_layout(pad=1,rect=(0,0,1,.95))


def xslice(tmax,x0=0,N=4,lim=None,centered=True, ratiomin=0.004, ratiomax=100, Nk=False):
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
    yloc = np.logspace(np.log10(1),np.log10(tmax),Npoints)
    for a,i in enumerate(ratios):
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))
        x=np.linspace(-L,L/2,500)
        X, Yloc=np.meshgrid(x,yloc)
        if Nk:
            k0=2.*np.pi*Nk/L
        if centered:
            Z= np.max(Solution(L,R,X+L/2+x0,Yloc,k0)[0],axis=1)
        else:
            Z= np.max(Solution(L,R,X+x0,Yloc,k0)[0],axis=1)
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


def xsliceMax(tmax,N=4,lim=None, ratiomin=0.004, ratiomax=100, Nk=False):
    """
    Maximum amplitude along z of the solution(z,t) for various size ratios

    tmax: maximum value of c*time in plot
    N: number of different size ratios
    lim: puts the same specified upper limit for all graphs, while by default each graph is scaled on its own
    ratiomin: lowest size ratio
    ratiomax: highest size ratio
    Nk: number of periods in the electromagnetic source. Spatial frequency for all aspect ratio is k0 if declared False.
    """
    fig = plt.figure()
    ax=[]

    ratios = np.logspace(np.log10(ratiomin),np.log10(ratiomax),N)
    yloc = np.logspace(np.log10(1),np.log10(tmax),Npoints)
    for a,i in enumerate(ratios):
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))
        x=np.linspace(-L,L,500)
        X, Yloc=np.meshgrid(x,yloc)
        if Nk:
            k0=2.*np.pi*Nk/L
        Z= np.max(Solution(L,R,X,Yloc,k0)[0],axis=1)
        Zm= np.min(Solution(L,R,X,Yloc,k0)[0],axis=1)
        line, = plt.semilogx(yloc, Khi*I*(l0**2)*(Z-Zm)/c)
        if lim:
            line.axes.set_ylim(-lim,lim)
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
            Z= P*Solution(L,R,L/2+x0*np.ones(np.shape(yloc)),yloc,k0)[0] /R**2
        else:
            Z= P*Solution(L,R,x0*np.ones(np.shape(yloc)),yloc,k0)[0]/R**2
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


def xslicePMax(tmax,N=4,lim=None, ratiomin=1e-5, ratiomax=0.1, P=1e15, Nk=False):
    """
    Cross-section at distance x0 from z=L/2 and at set power P of the solution(z,t) for various size ratios

    tmax: maximum value of c*time in plot
    N: number of different size ratios
    lim: puts the same specified upper limit for all graphs, while by default each graph is scaled on its own
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
        x=np.linspace(-L,L,500)
        X, Yloc=np.meshgrid(x,yloc)
        if Nk:
            k0=2.*np.pi*Nk/L
        Z= np.max(Solution(L,R,X,Yloc,k0)[0],axis=1)
        Zm= np.min(Solution(L,R,X,Yloc,k0)[0],axis=1)
        line, = plt.semilogx(yloc, Khi*P*(Z-Zm)/R**2 /c)
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



def maxStat(t,rmin,rmax,N):
    """
    Maximum of solution(z,ct) for a range of aspect ratios at set time t
    A figure needs to be set beforehand for this to plot

    t: c*time of evaluationnp.cbrt(V/(r*np.pi))
    rmin: minimum power of 10 for the aspect ratio
    rmax: maximum power of 10 for the aspect ratio
    N: number of points
    """
    rs=np.logspace(rmin,rmax,N)
    maxlist=[]
    minlist=[]
    for i in rs:
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))
        x=np.linspace(-L,L,500)
        Z=np.max(Solution(L,R,x,t*np.ones(np.shape(x)),k0)[0])
        Zm=np.min(Solution(L,R,x,t*np.ones(np.shape(x)),k0)[0])
        maxlist.append(Khi*I*(l0**2)*Z/c)
        minlist.append(Khi*I*(l0**2)*Zm/c)
    line, = plt.semilogx(rs,maxlist)
    plt.semilogx(rs,minlist)
    line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    line.set_label('ct={:.1e}m'.format(l0*t))
    plt.xlabel("Aspect ratio L/R")
    plt.ylabel('Maximum of h')
    plt.legend()

def maxStatKvar(t,rmin=-4,rmax=4,Npoints=500,wlmin=-1,wlmax=1,N=4):
    """
    Maximum of solution(z,ct) for a range of aspect ratios at set time t,
    and for several wavelenghts.

    t: c*time of evaluationnp.cbrt(V/(r*np.pi))
    rmin: minimum power of 10 for the aspect ratio
    rmax: maximum power of 10 for the aspect ratio
    Npoints: number of points
    wmin: minimum power of 10 for the wavelength
    wmax: maximum power of 10 for the wavelength
    N: number of wavelengths chosen
    """
    fig, ax=plt.subplots()
    rstart=np.logspace(rmin,rmax,Npoints)
    ws=np.logspace(wlmin,wlmax,N)
    ccycle=iter(plt.rcParams['axes.prop_cycle'])
    Ls=np.cbrt(V*rstart**2 /np.pi)
    for j in ws:
        maxlist=[]
        minlist=[]
        stepcycle=next(ccycle)
        vcolor=stepcycle['color']
        vstyle=stepcycle['linestyle']
        kl=2.*np.pi/j
        rs=np.sqrt(np.pi/V * (Ls-Ls%(j/2))**3)
        rs=np.delete(rs,np.where(rs==0))
        for i in rs:
            r=i
            L=np.cbrt(V*r**2 /np.pi)
            R=np.cbrt(V/(r*np.pi))
            x=np.linspace(-L,L,500)
            Z=np.max(Solution(L,R,x,t*np.ones(np.shape(x)),kl)[0])
            Zm=np.min(Solution(L,R,x,t*np.ones(np.shape(x)),kl)[0])
            maxlist.append(Khi*I*(l0**2)*Z/c)
            minlist.append(Khi*I*(l0**2)*Zm/c)
        line, = plt.semilogx(rs,maxlist,ls=vstyle,c=vcolor)
        plt.semilogx(rs,minlist,ls=vstyle,c=vcolor)
        ax.axvline(np.sqrt(np.pi*(j/2)**3/V),ls=vstyle,c=vcolor)
        line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        line.set_label(r'$\lambda$={:.1e}m'.format(l0*j))
        plt.legend()
    plt.xlabel("Aspect ratio L/R")
    plt.ylabel('Maximum of h')


def maxStat_kprop(t,rmin,rmax,N,Nk):
    """
    Maximum of solution(z,ct) for a range of aspect ratios at set time t, with number of oscillations
    staying constant to cylinder length.
    A figure needs to be set beforehand for this to plot

    t: c*time of evaluation
    rmin: minimum power of 10 for the aspect ratio
    rmax: maximum power of 10 for the aspect ratio
    N: number of points
    Nk: number of periods in the light cylinder (ex: wavelength=5 for L=20--> Nk=4)
    """
    rs=np.logspace(rmin,rmax,N)
    maxlist=[]
    for i in rs:
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))
        x=np.linspace(-L,L,500)
        k0=2.*np.pi*Nk/L
        Z=np.max(Solution(L,R,x,t*np.ones(np.shape(x)),k0)[0])
        Zm=np.min(Solution(L,R,x,t*np.ones(np.shape(x)),k0)[0])
        maxlist.append(Khi*I*(l0**2)*(Z-Zm)/c)
    line, = plt.semilogx(rs,maxlist)
    line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    line.set_label(r'ct={:.0f}$\mu$m, L/$\lambda$={}'.format(t,Nk))
    plt.xlabel("Aspect ratio L/R")
    plt.ylabel('Maximum of h')
    plt.legend()


def maxStat_kpropLvar(t,Lmin,Lmax,N,Nk,R=50):
    """
    Maximum of solution(z,ct) for a range of cylinder lengths at set time t, with number of oscillations
    staying constant to cylinder length.
    A figure needs to be set beforehand for this to plot

    t: c*time of evaluation
    Lmin: minimum power of 10 for the aspect ratio
    Lmax: maximum power of 10 for the aspect ratio
    N: number of points
    Nk: number of periods in the light cylinder (ex: wavelength=5 for L=20--> Nk=4)
    R: cylinder radius
    """
    Ls=np.logspace(Lmin,Lmax,N)
    maxlist=[]
    for i in Ls:
        x=np.linspace(-i,i,500)
        k0=2.*np.pi*Nk/i
        Z=np.max(Solution(i,R,x,t*np.ones(np.shape(x)),k0)[0])
        Zm=np.min(Solution(i,R,x,t*np.ones(np.shape(x)),k0)[0])
        maxlist.append(Khi*I*(l0**2)*(Z-Zm)/c)
    line, = plt.semilogx(Ls,maxlist)
    line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    line.set_label(r'ct={:.0f}$\mu$m, L/$\lambda$={}'.format(t,Nk))
    plt.xlabel("Cylinder length L")
    plt.ylabel('Maximum of h')
    plt.legend()

def maxk(t,wlmin,wlmax,N,scale='log',L=20,R=5):
    """
    Maximum of solution(z,ct) for a range of wavelengths dividing L at set time t.
    A figure needs to be set beforehand for this to plot

    t: c*time of evaluation
    wlmin: minimum power of 10 for the wavelength
    wlmax: maximum power of 10 for the wavelength
    N: number of points
    scale: sets the scale for the 'y' axis which presents the results
    L: cylinder length
    R: cylinder radius
    """
    if scale!='linear' and scale!='log':
        raise TypeError("'scale' must be str, options available are 'linear' or 'log'")
    wls=np.logspace(wlmin,wlmax,N)
    wls=2*L/np.floor(2*L/wls) #transforms the array into an array where wavelengths are dividers of L/2
    x=np.linspace(-L,L,3500)
    maxlist=[]
    for i in wls:
        k0=2.*np.pi/i
        Z=np.max(Solution(L,R,x,t*np.ones(np.shape(x)),k0)[0])
        Zm=np.min(Solution(L,R,x,t*np.ones(np.shape(x)),k0)[0])
        maxlist.append(Khi*I*(l0**2)*(Z-Zm)/c)

    if scale=='linear':
        line, = plt.semilogx(wls,maxlist)
        line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    elif scale=='log':
        line, = plt.loglog(wls,maxlist)
#    line.axes.axvline(2*L)

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
    plt.xlabel("Wavelength")
    plt.ylabel('Maximum of h')
    plt.legend()

def maxkStat(t,wlmin=-4,wlmax=4,Npoints=500,rmin=-4,rmax=4,N=4):
    """
    Maximum of solution(z,ct) for a range of wavelengths at set time t,
    and for several aspect ratios.
    Needs to have function 'maxk' loaded to work properly.

    t: c*time of evaluation
    wlmin: minimum power of 10 for the wavelength
    wlmax: maximum power of 10 for the wavelength
    Npoints: number of points

    rmin: minimum power of 10 for the aspect ratio
    rmax: maximum power of 10 for the aspect ratio
    N: number of different aspect ratios
    """
    fig, ax=plt.subplots()
    rs=np.logspace(rmin,rmax,N)

    ccycle=iter(plt.rcParams['axes.prop_cycle'])
    for i in rs:
        r=i
        L=np.cbrt(V*r**2 /np.pi)
        R=np.cbrt(V/(r*np.pi))
        maxk(t,wlmin,wlmax,Npoints,scale='log',L=L,R=R)

        stepcycle=next(ccycle)
        vcolor=stepcycle['color']
#        vstyle=stepcycle['linestyle']
        ax.axvline(2*L,c=vcolor)



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
        ax.set_xlabel(u'Z (µm)',labelpad=40)
        ax.set_ylabel(u'ct (µm)',labelpad=40)
        ax.set_zlabel('h',labelpad=40)
        ax.tick_params(axis='z', pad=10)

        xs = np.linspace(domain[0],domain[1],N)
        ys = np.linspace(domain[2],domain[3],N)

        Xs, Ys = np.meshgrid(xs,ys)

        Z= Solution(L,R,Xs,Ys,k0)[0]


        surf1=ax.plot_surface(Xs,Ys,Khi*I*(l0**2)*Z/c , cmap= cm.terrain) #3D Plot of the solution
        surfp=ax2.contourf(Xs,Ys,Khi*I*(l0**2)*Z/c,100)                #Colormap plot of the solution

        cbar=fig.colorbar(surf1)
        cbar2=fig2.colorbar(surfp)

        cbar.ax.set_ylabel('h', labelpad=40, rotation=270)
        cbar2.ax.set_ylabel('h', labelpad=40, rotation=270)

        fig.tight_layout(pad=1,rect=(0,0,1,.95))
        fig2.tight_layout(pad=1,rect=(0,0,1,.95))


def maxPcte(t,N,Rmin=1,Rmax=1e5,Lmin=5,Lmax=5e3,NL=4,P=1e15,lim=None,Nk=False):
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
    Nk: if equal to a number sets the number of oscillations in the cylinder for each L. Disabled if False.
    """
    Rs=np.logspace(np.log10(Rmin),np.log10(Rmax),N)
    Ls=np.logspace(np.log10(Lmin),np.log10(Lmax),NL)

    if lim: fig, ax = plt.subplots(nrows=int(np.ceil(np.sqrt(NL))), ncols=int(np.ceil(np.sqrt(NL))), sharex=True, sharey=True)
    else: fig, ax = plt.subplots(nrows=int(np.ceil(np.sqrt(NL))), ncols=int(np.ceil(np.sqrt(NL))), sharex=True)

    for i in range(len(ax)): ax[i,0].set_ylabel('h')
    for i in range(len(ax[-1])): ax[-1,i].set_xlabel(u'R (µm)')
    ax=ax.flatten()
    for a,i in enumerate(Ls):

        x=np.linspace(-i,i,500)
        Z=np.zeros(N)
        Zm=np.zeros(N)
        if Nk:
            k0=2.*np.pi*Nk/i
        for b,j in enumerate(Rs):
            Iloc=P/(np.pi*j**2)
            Z[b]= Iloc*np.max(Solution(i,j,x,t*np.ones(np.shape(x)),k0)[0])
            Zm[b]= Iloc*np.min(Solution(i,j,x,t*np.ones(np.shape(x)),k0)[0])

        ax[a].loglog(Rs, Khi*(Z-Zm)/c)
#        ax[a].loglog(Rs, -Khi*Zm/c)
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


def maxPcteRt(N,Rmin=1,Rmax=1e5,Lmin=5,Lmax=5e3,NL=4,P=1e15,lim=None,Nk=False):
    """
    Maximum of deformation for a set power as a function of the cylinder radius R, for various cylinder lengths L.
    The c*time at which this maximum is evaluated is equal to double the Rayleigh length

    N: number of points for R
    Rmin: minimum value of R
    Rmax: maximum value of R
    Lmin: minimum value of L
    Lmax: maximum value of L
    NL: number of different values of L

    P: set power taken
    lim: puts the same specified upper limit for all graphs, while by default each graph is scaled on its own
    Nk: if equal to a number sets the number of oscillations in the cylinder for each L. Disabled if False.
    """
    Rs=np.logspace(np.log10(Rmin),np.log10(Rmax),N)
    Ls=np.logspace(np.log10(Lmin),np.log10(Lmax),NL)

    fig, ax = plt.subplots(nrows=int(np.ceil(np.sqrt(NL))), ncols=int(np.ceil(np.sqrt(NL))), sharex=True, sharey=True)
    for i in range(len(ax)): ax[i,0].set_ylabel('h')
    for i in range(len(ax[-1])): ax[-1,i].set_xlabel(u'R (µm)')
    ax=ax.flatten()
    for a,i in enumerate(Ls):

        x=np.linspace(-i,i,500)
        Z=np.zeros(N)
        Zm=np.zeros(N)
        if Nk:
            k0=2.*np.pi*Nk/i
        for b,j in enumerate(Rs):
            Iloc=P/(np.pi*j**2)
            t=k0*j**2
            Z[b]= Iloc*np.max(Solution(i,j,x,t*np.ones(np.shape(x)),k0)[0])
            Zm[b]= Iloc*np.min(Solution(i,j,x,t*np.ones(np.shape(x)),k0)[0])

        ax[a].semilogx(Rs, Khi*(Z-Zm)/c)
#        ax[a].semilogx(Rs, -Khi*Zm/c)
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
#        ax[a].yaxis.set_major_locator(ticker.MaxNLocator(3))

#        ax[a].set_xlabel(u'R (µm)')
#        ax[a].set_ylabel('h')
#    plt.suptitle('Nk={}'.format(Nk))
    plt.tight_layout(pad=1,rect=(0,0,1,.95))

def maxMaxPcteRt(Nkmin=.5,Nkmax=1e4,Npoints=1000,P=1e15,scale='lin'):
    """
    Maximum of deformation for a set power as a function of the number of optical cycles Nk.
    The c*time at which this maximum is evaluated is equal to double the Rayleigh length
    for a cylinder radius Re such that  Re>Le the cylinder length.

    Nkmin: minimum value of R
    Nkmax: maximum value of R
    Npoints: number of different values of Nk
    
    P: set power taken
    scale: defines if the y-axis has a linear or logarithmic scale. Set by default to linear.
    """
    fig=plt.figure()
    Le=500
    Re=1e4
    Iloc=P/(np.pi*Re**2)
    Nkunform=np.logspace(np.log10(Nkmin),np.log10(Nkmax),Npoints)
    Nks=Nkunform- Nkunform%.5 #ensures that every wavelength is a divider of 2*L
    Z=np.zeros(Npoints)
    Zm=np.zeros(Npoints)
    for b,i in enumerate(Nks):
        k0=2.*np.pi*i/Le
        t=Re**2 *k0
#        x=np.linspace(-Le,Le,500)
        x=np.linspace(-3*Le/i +Le/2,Le/2+ 2*Le/i,600) #evaluating the maximum over a range of 4 periods centered in the middle of the cylinder
        Z[b]= Iloc*np.max(Solution(Le,Re,x,t*np.ones(np.shape(x)),k0)[0])
        Zm[b]= Iloc*np.min(Solution(Le,Re,x,t*np.ones(np.shape(x)),k0)[0])

    if scale=='log':
        line, =plt.loglog(Nks,Khi*(Z-Zm)/c)
    else:
        line, =plt.semilogx(Nks,Khi*(Z-Zm)/c)
        line.axes.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

    plt.xlabel(r'L/$\lambda$')
    plt.ylabel('Maximum of h')
    plt.tight_layout(pad=1,rect=(0,0,1,.95))
