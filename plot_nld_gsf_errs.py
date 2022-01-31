#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 17:15:38 2021

@author: francesco, updated 31st January 2022

Draw nld and gsf from the nld_whole.txt and gsf_whole.txt files produced from make_nlds_gsfs_lists.py
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
sys.path.insert(1, '/home/francesco/Documents/Talys-calculations')
from readlib import *
from systlib import *

#parameters. play with these
L1max = 23
L2max = 23
enbin = 22
errorbars = True
conf90=False
confsigmas = True

#Don't modify unless you know what you're doing
blist = np.linspace(0.40,1.4,50)
base_rho = 375720
Sn = 8.380
a0 = -0.7599
a1 = 0.1562
sig = 0.341
limit_points = np.linspace(0.0,1.0,21)

ranges = [[-1,6], [np.log10(5e-9),-6]]
database_path = 'Make_dataset/127Sb-database/'

#import experimental nld and gsf
nld_mat = np.genfromtxt('nld_whole.txt', unpack = True).T
gsf_mat = np.genfromtxt('gsf_whole.txt', unpack = True).T
#delete rows with nans
nld_mat = nld_mat[~np.isnan(nld_mat).any(axis=1)]
gsf_mat = gsf_mat[~np.isnan(gsf_mat).any(axis=1)]

#delete some points?
#delete last 3 rows in nld
nld_mat = np.delete(nld_mat,[-3,-2,-1],0)
#delete last rows in gsf
gsf_mat = np.delete(gsf_mat,[-1],0)

#import known levels
known_levs = import_ocl(database_path + '040-150288/L1-1_L2-5/rholev.cnt',a0,a1)

#import TALYS calculated GSFs
TALYS_strengths = [readstrength(51, 127, 1, 1, strength, 1) for strength in range(1,9)]

#start plotting
cmap = matplotlib.cm.get_cmap('Dark2')
fig,doubleaxs = plt.subplots(nrows = 1, ncols = 2)
fig0,ax0 = plt.subplots()
fig1,ax1 = plt.subplots()
singleaxs = [ax0,ax1]
fig.set_size_inches(20, 10)
chi2_lim = [6,10]
chi2_lim2 = [14,18]
rhos = [b*base_rho for b in blist]
for ax in (doubleaxs[0], ax0):
    ax.axvspan(known_levs[int(chi2_lim[0]),0], known_levs[int(chi2_lim[1]),0], alpha=0.2, color='red')
    ax.axvspan(known_levs[int(chi2_lim2[0]),0], known_levs[int(chi2_lim2[1]),0], alpha=0.2, color='red',label='Fitting intervals')
    ax.plot(known_levs[:,0],known_levs[:,1],'k-',label='Known levels')
    middlerho = base_rho*(0.53+1.28)/2
    errorrho = base_rho*1.28-middlerho
    ax.errorbar(Sn, middlerho,yerr=errorrho,ecolor='g',linestyle=None, elinewidth = 4, capsize = 5, label=r'$\rho$ at Sn')

#Plot experiment data
for doubleax, singleax, val_matrix in zip(doubleaxs, singleaxs, [nld_mat, gsf_mat]):
    for ax in (doubleax, singleax):
        ax.plot(val_matrix[:,0], val_matrix[:,1], 'bo', label = 'Oslo data')
        if errorbars:
            lower_err = val_matrix[:,1] - val_matrix[:,3]
            upper_err = val_matrix[:,-2] - val_matrix[:,1]
            ax.errorbar(val_matrix[:,0], val_matrix[:,1],yerr=[lower_err, upper_err],color = 'b', ecolor='b')
        else:
            if conf90:
                ax.fill_between(val_matrix[:,0], val_matrix[:,3], val_matrix[:,-2], color = 'c', alpha = 0.2, label='90% confidence')
            elif confsigmas:
                ax.fill_between(val_matrix[:,0], val_matrix[:,2], val_matrix[:,-1], color = 'c', alpha = 0.2, label=r'2$\sigma$ confidence')
                ax.fill_between(val_matrix[:,0], val_matrix[:,3], val_matrix[:,-2], color = 'b', alpha = 0.2, label=r'1$\sigma$ confidence')
            else:
                for conf in range(int(len(limit_points)/2 - 0.5)):
                    ax.fill_between(val_matrix[:,0], val_matrix[:,2+conf], val_matrix[:,-(conf+1)], color = cmap((conf*2 + 1)/len(limit_points)), alpha = 0.2, label=str(int((limit_points[-(conf+1)] - limit_points[conf])*100)) + '% confidence')
        
        ax.vlines(x=Sn, ymin=0, ymax = 1e6, colors='purple', ls='--', alpha = 0.5, label=r'$S_n$')
        ax.set_xlabel('Energy [MeV]')
        ax.set_yscale('log')

#plot TALYS strengths
for i, TALYS_strength in enumerate(TALYS_strengths):
    stl = '--'
    doubleaxs[1].plot(TALYS_strength[:,0],TALYS_strength[:,2], color = cmap(i/8), linestyle = stl, label = 'strength %d'%(i+1))
    ax1.plot(TALYS_strength[:,0],TALYS_strength[:,2], color = cmap(i/8), linestyle = stl, label = 'strength %d'%(i+1))

#import GLO parameters for Te and Sn and calculate GLO parameters of 127Sb 
Te128_par = np.loadtxt('128Te_params_root')
Sn126_par = np.loadtxt('126Sn_params')
Sb127_par = (Te128_par + Sn126_par)/2

#plot
energies = np.linspace(0, 20, 1000)
for axs in (singleaxs, doubleaxs):
    axs[0].set_ylabel(r'NLD [MeV$^{-1}$]')
    axs[1].set_ylabel(r'GSF [MeV$^{-3}$]')
    axs[0].set_ylim(10**ranges[0][0],10**ranges[0][1])
    axs[1].set_ylim(10**ranges[1][0],10**ranges[1][1])
    axs[1].set_xlim(0.3,9)
    for ax in axs:
        #ax.grid()
        ax.legend(loc = 'upper left', ncol = 2)

fig.show()
fig0.show()
fig1.show()