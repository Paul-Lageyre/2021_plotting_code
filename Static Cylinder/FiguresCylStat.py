#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 11:47:58 2021

@author: paulg
"""

from TestFiguresHeaviStat import *

directory='/home/paulg/Documents/Redactions/Courbes/VisuHeaviStat/Corrected_07_2021' #fill with the directory you want the graphs to save to

def saveplt(name,direc=directory):
    
    plt.savefig(direc + '/png/' + name + '.png', format='png')
    plt.savefig(direc + '/pdf/' + name + '.pdf', format='pdf')
    plt.close()
    

tslice(10,lim=1.56e-36)
saveplt('tslice10')

tslice(50,lim=1.56e-36)
saveplt('tslice50')

tslice(1000,lim=1.56e-36)
saveplt('tslice1000')

xslice(150)
saveplt('xslice150')


maxStat(-4,4,500)
a=plt.gca()
a.set_ylim(0,1.66e-36)
saveplt('maxStat')

tFWHM(1000,-4,2,1000)
saveplt('tFWHM')

surfplot(20,5,500,domain=[-50,70,0,100])
saveplt('Ex3DFlat')

a=plt.gca()
a.view_init(38.0,-90.0)
saveplt('Ex3DFront')

surfplot(20,5,500,domain=[-50,70,0,100])
plt.close()
a=plt.gca()
a.view_init(45.4,-50.3)
saveplt('Ex3DSide')