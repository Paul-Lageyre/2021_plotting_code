#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 11:47:58 2021

@author: paulg
"""

from TestFiguresHeaviMob import *

directory='/home/paulg/Documents/Redactions/Courbes/VisuHeaviMob/Corrected_05_2021' #fill with the directory you want the graphs to save to

def saveplt(name,direc=directory):
    
    plt.savefig(direc + '/png/' + name + '.png', format='png')
    plt.savefig(direc + '/pdf/' + name + '.pdf', format='pdf')
    plt.close()
    

tslice(50,lim=1.55e-36)
saveplt('tslice50')

tslice(250,lim=4.9e-36)
saveplt('tslice250')

tslice(5000,lim=1.65e-35)
saveplt('tslice5000')

xslice(50,lim=1.5e-36,centered=False)
saveplt('xslice50uncen')

xslice(5000,lim=1.6e-35,centered=False)
saveplt('xslice5000uncen')

xsliceP(1e8,centered=False,lim=1.2e-36)
saveplt('xsliceP1e8uncen')

plt.figure()
maxMob(50,-10,4,500)
maxMob(250,-10,4,500)
maxMob(1500,-10,4,500)
maxMob(5000,-10,4,500)
a=plt.gca()
a.set_ylim(0,2.2e-35)
saveplt('maxMob50-250-1500-5000')

maxPcte(1e4,500,lim=6.2e-37)
saveplt('maxPcte10000')

maxPcteRt(500,NL=4)
saveplt('maxPcteRt')

surfplot(20,5,500,domain=[-70,30,0,100])
saveplt('Ex3DFlat')

a=plt.gca()
a.view_init(39,-47.8)
saveplt('Ex3DFront')

surfplot(20,5,500,domain=[-70,30,0,100])
plt.close()
a=plt.gca()
a.view_init(23.9,-108.7)
saveplt('Ex3DBack')