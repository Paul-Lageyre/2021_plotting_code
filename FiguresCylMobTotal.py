#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 11:47:58 2021

@author: paulg
"""

from TestFiguresHeaviMobTotal import *

directory='/home/paulg/Documents/Redactions/Courbes/VisuHeaviMobTotal/06_09_21' #fill with the directory you want the graphs to save to

def saveplt(name,direc=directory):
    
    plt.savefig(direc + '/png/' + name + '.png', format='png')
    plt.savefig(direc + '/pdf/' + name + '.pdf', format='pdf')
    plt.close()
    

tslice(250,ratiomin=0.177, ratiomax=4)
saveplt('tslice250')

tslice(5000,ratiomin=0.177, ratiomax=4)
saveplt('tslice5000')

tslice(50,ratiomin=0.177, ratiomax=4,lim=4.4e-36,lim_inf=0)
saveplt('tslice50lim')

tslice(250,ratiomin=0.177, ratiomax=4, lim=1.1e-35, lim_inf=0)
saveplt('tslice250lim')

tslice(5000,ratiomin=0.177, ratiomax=4,lim=3.1e-35,lim_inf=0)
saveplt('tslice5000lim')



#k0=2.*np.pi/1e7
#tslice(50,lim=3.2e-36,lim_inf=0)
#saveplt('tslice_kSmall50')
#
#tslice(5000,lim=3.3e-35,lim_inf=0)
#saveplt('tslice_kSmall5000')


#xsliceMax(1e3,Nk=1)
#saveplt('xsliceMaxNk1_1e3')
#
#xsliceMax(1e5,Nk=1)
#saveplt('xsliceMaxNk1_1e5')


plt.figure()
maxStat_kprop(1500,-4,4,500,.5)
maxStat_kprop(250,-4,4,500,.5)
maxStat_kprop(50,-4,4,500,.5)
maxStat_kprop(10,-4,4,500,.5)
a=plt.gca()
a.set_ylim(0,2.7e-35)
saveplt('maxStat_kprop_ctvar')

plt.figure()
maxStat_kprop(250,-3,4,500,.5)
maxStat_kprop(250,-3,4,500,1)
maxStat_kprop(250,-3,4,500,2.5)
maxStat_kprop(250,-3,4,500,20)
a=plt.gca()
a.set_ylim(0,1.1e-35)
saveplt('maxStat_kprop_Nkvar')

maxkStat(250,wlmin=-4,wlmax=3,rmin=-2,rmax=2)
saveplt('maxkStat_t250_closeup')

maxPcteRt(500,Rmin=1e-1,Nk=.5)
saveplt('maxPcteRtOscill')

maxMaxPcteRt()
saveplt('maxMaxPcteRt')


k0=2.*np.pi/10

surfplot(20,5,500,k0,domain=[-30,25,0,100])
saveplt('Ex2D_5mum')

a=plt.gca()
a.view_init(26.9,-94.2)
saveplt('Ex3D_5mum')
