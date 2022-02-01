#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Tue Oct  5 15:45:14 2021

@author: francesco, updated 1st February 2022

Find the paramters for the GLO fitting the GDR for 127Sb, as the average of the 
ones calculated for Sn126 and Te128 (separate files, produced in Tin_gsf-fit.py 
and fitRoot_gSF_128Te_oneGauss.cpp and imported from files written by the codes).
Then, import and plot 126Sn and 128Te from Darmstadt and tin for reference.
'''

import numpy as np
import matplotlib.pyplot as plt
from systlib import load_known_gsf, GLO

#import GLO parameters for Te and Sn
Te128_par = np.loadtxt('data/generated/128Te_params_root')
Sn126_par = np.loadtxt('data/generated/126Sn_params')

#Calculate GLO parameters of 127Sb as the average of Te and Sn
Sb127_par = np.zeros_like(Te128_par)
Sb127_par[:,0] = (Te128_par[:,0] + Sn126_par[:,0])/2
Sb127_par[:,1] = (Te128_par[:,1] + Sn126_par[:,1])

np.savetxt('data/generated/127Sb_params', Sb127_par)

#load other nuclei
tins = [load_known_gsf(A,50,lab=lab) for A in (120,124) for lab in ('o','d')]
Te128 = load_known_gsf(128,52)

labs = ['E0', 'Gamma', 'sigma', 'T']
for i in range(4):
    print(labs[i] + ': %s +- %s'%(Sb127_par[i,0],Sb127_par[i,1]))

#plot the data
energies = np.linspace(0, 20, 1000)
plt.plot(energies, GLO(energies, Sb127_par[3,0], Sb127_par[0,0], Sb127_par[1,0], Sb127_par[2,0]), label = 'Sb127_GLO')
[tin.plot() for tin in tins]
Te128.plot()
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()