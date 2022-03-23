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
from readlib import readreaclib, nonsmoker, Z2Name, rate2MACS, readastro, xs2MACS, xs2rate
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
<<<<<<< HEAD
rates = False
=======
rates = True
>>>>>>> 1caa57bc5e6f56f0bd01a8f8f793865ef0b0a907
Sb127_mass = 125.907247480 #in a.u.

#constants
k_B = 8.617333262145e1 #Boltzmann konstant, keV*GK^-1

#useful arrays
T9range = np.array([0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
T9extrange = np.array([0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
keVrange = T9extrange*k_B



#Make TALYS uncertainty span
dummyastro = readastro(nucleus, A, nld[0], mass[0], strength[0], omp[0])
if rates:    
    val_index = 1
else:
    val_index = 2
max_val = np.zeros_like(dummyastro[:,val_index])
min_val = dummyastro[:,val_index].copy()
for ld in nld:
    for m in mass:
        for s in strength:
            for o in omp:
                currentastro = readastro(nucleus, A, ld, m, s, o)
                alpha = 0.2
                for i,val in enumerate(currentastro[:,val_index]):
                    if val > max_val[i]:
                        max_val[i] = val
                    if val < min_val[i]:
                        min_val[i] = val
                    
<<<<<<< HEAD
#Import Oslo data, statistical and systematic errors
if rates:
    Oslo_mat = np.genfromtxt('data/generated/ncrates_whole.txt', unpack = True).T
    x_Oslo = Oslo_mat[:,0]
    Oslo_stat = np.loadtxt('data/generated/ncrates_stats.txt')
else:
    Oslo_mat = np.genfromtxt('data/generated/MACS_whole.txt', unpack = True).T
    x_Oslo = Oslo_mat[:,0]*k_B #keV
    Oslo_stat = np.loadtxt('data/generated/MACSs_stats.txt')
=======
#Import Oslo data
if rates:
    Oslo_mat = np.genfromtxt('data/generated/ncrates_whole.txt', unpack = True).T
    x_Oslo = Oslo_mat[:,0]
else:
    Oslo_mat = np.genfromtxt('data/generated/MACS_whole.txt', unpack = True).T
    x_Oslo = Oslo_mat[:,0]*k_B #keV

>>>>>>> 1caa57bc5e6f56f0bd01a8f8f793865ef0b0a907

#Import other models
a = readreaclib(nucleus, A, reaclib_file = 'data/ncapture/reaclib')
val_reaclib = nonsmoker(a, T9range) #rates
tendl_data = np.loadtxt('data/ncapture/tendl_ng_rates.tot')
mat_bruslib = np.loadtxt('data/ncapture/bruslib-126Sb')
other_libs_df = pd.read_csv('data/ncapture/MACS_from_libs.txt', sep='	', header = 1)
models = [other_libs_df[other_libs_df['modn']==j] for j in range(1,7)]
if rates:
    x_reaclib = T9range
    x_TENDL = tendl_data[:,0]*1e3/k_B
    val_TENDL = xs2rate(tendl_data[:,0]*1e3,tendl_data[:,1],Sb127_mass,Ts=x_TENDL)
    x_bruslib = mat_bruslib[:,0]
    val_bruslib = mat_bruslib[:,1]
else:
    x_reaclib = T9range*k_B
    val_reaclib = rate2MACS(val_reaclib,Sb127_mass,T9range)
    x_TENDL = tendl_data[:,0]*1e3
    val_TENDL = xs2MACS(x_TENDL,tendl_data[:,1],Ts = x_TENDL/k_B)
    x_bruslib = mat_bruslib[:,0]*k_B
    val_bruslib = rate2MACS(mat_bruslib[:,1],Sb127_mass,mat_bruslib[:,0])

#plot figure
cmap = matplotlib.cm.get_cmap('YlGnBu')
fig, ax = plt.subplots()
ax.fill_between(x_Oslo, min_val, max_val, color = cmap(1/8), alpha = 1, label = 'TALYS uncertainty span')
if errorbars:
<<<<<<< HEAD
    lower_err = Oslo_mat[:,1] - Oslo_mat[:,3] - (Oslo_stat[:,1] - Oslo_stat[:,2])
    upper_err = Oslo_mat[:,-2] - Oslo_mat[:,1] + (Oslo_stat[:,3] - Oslo_stat[:,1])
    ax.errorbar(x_Oslo, Oslo_mat[:,1],yerr=[lower_err, upper_err],ecolor='k')
else:
    lowerline_2s = Oslo_mat[:,2] - (Oslo_stat[:,1] - Oslo_stat[:,2])
    upperline_2s = Oslo_mat[:,5] + (Oslo_stat[:,3] - Oslo_stat[:,1])
    lowerline_1s = Oslo_mat[:,3] - (Oslo_stat[:,1] - Oslo_stat[:,2])
    upperline_1s = Oslo_mat[:,4] + (Oslo_stat[:,3] - Oslo_stat[:,1])
    ax.fill_between(x_Oslo, lowerline_2s, upperline_2s, color = cmap(2/8), alpha= 1, label = r'Oslo data, $2\sigma$')
    ax.fill_between(x_Oslo, lowerline_1s, upperline_1s, color = cmap(3/8), alpha= 1, label = r'Oslo data, $1\sigma$')
=======
    lower_err = Oslo_mat[:,1] - Oslo_mat[:,3]
    upper_err = Oslo_mat[:,-2] - Oslo_mat[:,1]
    ax.errorbar(x_Oslo, Oslo_mat[:,1],yerr=[lower_err, upper_err],ecolor='k')
else:
    ax.fill_between(x_Oslo, Oslo_mat[:,2], Oslo_mat[:,5], color = cmap(2/8), alpha= 1, label = r'Oslo data, $2\sigma$')
    ax.fill_between(x_Oslo, Oslo_mat[:,3], Oslo_mat[:,4], color = cmap(3/8), alpha= 1, label = r'Oslo data, $1\sigma$')
>>>>>>> 1caa57bc5e6f56f0bd01a8f8f793865ef0b0a907
    ax.plot(x_Oslo, Oslo_mat[:,1], color = 'b', linestyle ='-', label = 'Oslo data')

for i, model in enumerate(models):
    label = model['model'].iloc[0]
    label = label[:-1]
    if rates:
        matr = model[['kT','ngamma']].to_numpy()
        x_libs = matr[:,0]/k_B
    else:
        matr = model[['kT','MACS']].to_numpy()
        x_libs = matr[:,0]
    if i==0:
        ax.plot(x_libs, matr[:,1], color = cmap(5/8), linestyle = '-.', label=label)

ax.plot(x_reaclib, val_reaclib, color = 'k', linestyle = '--',label = 'JINA REACLIB')
ax.plot(x_TENDL, val_TENDL, color = 'tab:brown', linestyle = '-.',label = 'TENDL')
ax.plot(x_bruslib, val_bruslib, color = cmap(8/8), linestyle = ':', label = 'BRUSLIB')

if rates:
    ax.set_title(r'rates for '+ ToLatex(str(A) + Z2Name(nucleus)) + '$(n,\gamma)$')
    ax.set_xlabel(r'$T$ [GK]')
    ax.set_ylabel(r'$(n,\gamma)$ rate [cm$^{3}$mol$^{-1}$s$^{-1}$]')
    #ax.grid()
    ax.set_yscale('log')
    ax.set_xlim([-0.1,110/k_B])
    ax.set_ylim([3e6,1.3e9])
    ax.legend(ncol=2)
else:
<<<<<<< HEAD
    #ax.set_title('MACS for '+ ToLatex(str(A) + Z2Name(nucleus)) + '$(n,\gamma)$')
=======
    ax.set_title('MACS for '+ ToLatex(str(A) + Z2Name(nucleus)) + '$(n,\gamma)$')
>>>>>>> 1caa57bc5e6f56f0bd01a8f8f793865ef0b0a907
    ax.set_xlabel(r'$k_B$ $T$ [keV]')
    ax.set_ylabel('MACS [mb]')
    #ax.grid()
    ax.set_yscale('log')
    ax.set_xlim([-2,110])
    ax.set_ylim([50,14e3])
    ax.legend(ncol=2)
fig.show()