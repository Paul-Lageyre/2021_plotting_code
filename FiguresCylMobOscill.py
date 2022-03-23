#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 11:47:58 2021

@author: paulg
"""

from TestFiguresHeaviMobOscill import *

directory='/home/paulg/Documents/Redactions/Courbes/VisuHeaviMobOscill/01_09_21' #fill with the directory you want the graphs to save to

def saveplt(name,direc=directory):
    
    plt.savefig(direc + '/png/' + name + '.png', format='png')
    plt.savefig(direc + '/pdf/' + name + '.pdf', format='pdf')
    plt.close()
    

tslice(250)
saveplt('tslice250')

tslice(5000)
saveplt('tslice5000')

tslice(250,lim=5.5e-36)
saveplt('tslice250lim')

tslice(5000,lim=3.3e-35)
saveplt('tslice5000lim')


tslice_kprop(50,4,lim=4.4e-37)
saveplt('tslice_kprop50')

tslice_kprop(250,4,lim=2.1e-36)
saveplt('tslice_kprop250')


k0=2.*np.pi/1e7
tslice(50,lim=3.2e-36,lim_inf=0)
saveplt('tslice_kSmall50')

tslice(5000,lim=3.3e-35,lim_inf=0)
saveplt('tslice_kSmall5000')


xsliceMax(1e3,Nk=1)
saveplt('xsliceMaxNk1_1e3')

xsliceMax(1e5,Nk=1)
saveplt('xsliceMaxNk1_1e5')


plt.figure()
maxStat_kprop(1500,-4,4,500,.5)
maxStat_kprop(250,-4,4,500,.5)
maxStat_kprop(50,-4,4,500,.5)
maxStat_kprop(10,-4,4,500,.5)
a=plt.gca()
a.set_ylim(0,2.15e-35)
saveplt('maxStat_kprop_ctvar')

plt.figure()
maxStat_kprop(250,-3,4,500,.5)
maxStat_kprop(250,-3,4,500,1)
maxStat_kprop(250,-3,4,500,2.5)
maxStat_kprop(250,-3,4,500,20)
a=plt.gca()
a.set_ylim(0,8.75e-36)
saveplt('maxStat_kprop_Nkvar')

maxkStat(250,wlmin=-2,wlmax=3,rmin=-2,rmax=2)
saveplt('maxkStat_t250_closeup')

maxPcteRt(500,Rmin=1e-1,Nk=.5)
saveplt('maxPcteRtOscill')

maxMaxPcteRt()
saveplt('maxMaxPcteRt')

surfplot(20,5,500,domain=[-10,25,0,100])
saveplt('Ex2D_5mum')

a=plt.gca()
a.view_init(21.3,-103.5)
saveplt('Ex3D_5mum')
