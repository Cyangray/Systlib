#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 12:58:32 2022

@author: francesco, updated 1st February 2022

Code plotting the MACS of 127Sb from MACS_whole.txt generated in make_nlds_gsfs_lists.py
and compares it to TALYS predictions, and other libraries
"""
import numpy as np
import matplotlib.pyplot as plt
#import sys
#sys.path.insert(1, '/home/francesco/Documents/Talys-calculations')
from readlib import readreaclib, nonsmoker, Z2Name, rate2MACS, readastro, xs2MACS
from systlib import ToLatex
import pandas as pd
import matplotlib

#inputs
Z= 51
nucleus = Z #isotope from which to extract strength function. one, or a series of names, or Zs
A = 126         #mass number of isotope
omp = [1,2]         #three options: 0, 1 and 2. 0 for no arguments, 1 for jlmomp-y and 2 for localomp-n
strength = [1,2,3,4,5,6,7,8]    #Which strength function to read - numbers from 1 to 8
nld = [1,2,3,4,5,6]
mass = [1,2,3]
errorbars = False
Sb127_mass = 125.907247480 #in a.u.

#constants
k_B = 8.617333262145e1 #Boltzmann konstant, keV*GK^-1

#useful arrays
T9range = np.array([0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
T9extrange = np.array([0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
keVrange = T9extrange*k_B

#Make REACLIB MACS
a = readreaclib(nucleus, A, reaclib_file = 'data/ncapture/reaclib')
ncrate_reaclib = nonsmoker(a, T9range)
MACS_reaclib = rate2MACS(ncrate_reaclib,Sb127_mass,T9range)

#Make TALYS uncertainty span
dummyastro = readastro(nucleus, A, nld[0], mass[0], strength[0], omp[0])
maxMACS = np.zeros_like(dummyastro[:,2])
minMACS = dummyastro[:,2].copy()
for ld in nld:
    for m in mass:
        for s in strength:
            for o in omp:
                currentastro = readastro(nucleus, A, ld, m, s, o)
                alpha = 0.2
                for i,MACS in enumerate(currentastro[:,2]):
                    if MACS > maxMACS[i]:
                        maxMACS[i] = MACS
                    if MACS < minMACS[i]:
                        minMACS[i] = MACS
                    
#Import Oslo data
MACS_mat = np.genfromtxt('data/generated/MACS_whole.txt', unpack = True).T

#Import other models
tendl_data = np.loadtxt('data/ncapture/tendl_ng_rates.tot')
rates_bruslib = np.loadtxt('data/ncapture/bruslib-126Sb')
other_rates_df = pd.read_csv('data/ncapture/MACS_from_libs.txt', sep='	', header = 1)
models = [other_rates_df[other_rates_df['modn']==j] for j in range(1,7)]
MACS_tendl = xs2MACS(tendl_data[:,0]*1e3,tendl_data[:,1])
kTs = MACS_mat[:,0]*k_B #keV

#plot figure
cmap = matplotlib.cm.get_cmap('YlGnBu')
fig, ax = plt.subplots()
ax.fill_between(kTs, minMACS, maxMACS, color = cmap(1/8), alpha = 1, label = 'TALYS uncertainty span')
if errorbars:
    lower_err = MACS_mat[:,1] - MACS_mat[:,3]
    upper_err = MACS_mat[:,-2] - MACS_mat[:,1]
    ax.errorbar(kTs, MACS_mat[:,1],yerr=[lower_err, upper_err],ecolor='k')
else:
    ax.fill_between(kTs, MACS_mat[:,2], MACS_mat[:,5], color = cmap(2/8), alpha= 1, label = r'Oslo data, $2\sigma$')
    ax.fill_between(kTs, MACS_mat[:,3], MACS_mat[:,4], color = cmap(3/8), alpha= 1, label = r'Oslo data, $1\sigma$')
    ax.plot(kTs, MACS_mat[:,1], color = 'b', linestyle ='-', label = 'Oslo data')

for i, model in enumerate(models):
    label = model['model'].iloc[0]
    label = label[:-1]
    matr = model[['kT','MACS']].to_numpy()
    if i==0:
        ax.plot(matr[:,0], matr[:,1], color = cmap(5/8), linestyle = '-.', label=label)
ax.plot(T9range*k_B, MACS_reaclib, color = 'k', linestyle = '--',label = 'JINA REACLIB')
ax.plot(keVrange, MACS_tendl, color = 'tab:brown', linestyle = '-.',label = 'TENDL')
ax.plot(rates_bruslib[:,0]*k_B, rate2MACS(rates_bruslib[:,1],Sb127_mass,rates_bruslib[:,0]), color = cmap(8/8), linestyle = ':', label = 'BRUSLIB')
ax.set_title('MACS for '+ ToLatex(str(A) + Z2Name(nucleus)) + '$(n,\gamma)$')
ax.set_xlabel(r'$k_B$ $T$ [keV]')
ax.set_ylabel('MACS [mb]')
#ax.grid()
ax.set_yscale('log')
ax.set_xlim([-2,110])
ax.set_ylim([50,14e3])
ax.legend(ncol=2)        
fig.show()