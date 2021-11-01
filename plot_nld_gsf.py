#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 14:24:41 2021

@author: francesco
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(1, '/home/francesco/Documents/Talys-calculations')
from readlib import *
from systlib import *


'''
Draw nld and gsf according to their chi2 as histograms.
must have produced chitables.npy in plot_nld.py
'''

#parameters. play with these
L1max = 23
L2max = 23
plot_gsf_hists = True
dens = False #choose true to plot the probability distribution for each energy, choose false to plot each energy with the highest probability set to 1
#Do you want to plot a specific energy as a probability distribution of the gsf?
plot_single_energy_hist = True
enbin=51
plotlist=[
    #[1,3],
    #[11,15],
    #[19,22],
    ]

#Don't modify unless you know what you're doing
blist = np.linspace(0.53,1.2804,50)
base_rho = 375720
rho_Sn_err = 90000
fermistart=37 #from where to plot the CT model fit
Sn = 8.380
a0 =  -0.7599
a1 =   0.1562

#Load the chitable file, and prepare to plot by initializing stuff
chitables = np.load('chitables.npy')


nlds = np.load('nlds.npy', allow_pickle = True)
gsfs = np.load('gsfs.npy', allow_pickle = True)

#weigh down? Increase the chi2 value, so that when inverted, the probability will be reduced as if it were normally distributed about Gg=105, sigma=25
for nld, gsf in zip(nlds,gsfs):
    ratio = normal_dist_weigh(105, 25, nld.Gg)
    nld.chi2 = gsf.chi2 =  nld.chi2/ratio

#initializing loop for plotting. First loop through the nlds, and then through the gsfs
data_matrices = []
rgba_colorss = []
ranges = [[-1,6], [np.log10(5e-9),-6]]
bins = 500
for lst, rng in zip([nlds,gsfs],ranges):
    data_matrix, rgba_colors = calc_data_matrix(lst, rng, bins, dens)

    data_matrices.append(data_matrix)
    rgba_colorss.append(rgba_colors)
    

fig,axs = plt.subplots(nrows = 1, ncols = 2)
fig.set_size_inches(20, 10)
labels = ['127Sb NLD histograms', '127Sb GSF histograms']
for ax, data_matrix, rgba_colors, label in zip(axs, data_matrices, rgba_colorss, labels):
    ax.scatter(data_matrix[:,0],10**data_matrix[:,1], marker = '.', c=rgba_colors, label = label)#, label=r'$\rho$=%d'%rhos[j])
    ax.vlines(x=Sn, ymin=0, ymax = 1e6, colors='purple', ls='--', alpha = 0.5, label='Sn')
    ax.set_xlabel('Energy [MeV]')
    ax.set_yscale('log')

#import and plot other nuclei
Sn124_d = load_known_gsf(124,'Sn',lab='d')
Sn124_o = load_known_gsf(124,'Sn',lab='o')
Te128 = load_known_gsf(128,'Te')
Sn124_o.plot(ax=axs[1],color='r')
Sn124_d.plot(ax=axs[1],color='r')
Te128.plot(ax=axs[1])

#import known levels
known_levs = import_ocl('../countingexp/80/053-199131/L1-5_L2-12_B-053_Gg-80/rholev.cnt',a0,a1)
axs[0].plot(known_levs[:,0],known_levs[:,1],'k-',label='Known levels')
rhos = [b*base_rho for b in blist]
axs[0].plot([Sn for x in range(len(rhos))], rhos, 'g.', label=r'$\rho$ at Sn')

#import GLO parameters for Te and Sn and calculate GLO parameters of 127Sb 
Te128_par = np.loadtxt('128Te_params')
Sn126_par = np.loadtxt('126Sn_params')
Sb127_par = (Te128_par + Sn126_par)/2

#plot

chi2_lim = [6,10]
chi2_lim2 = [14,18]
#chi2_lim_e = [known_levs[int(chi2_lim[0]),0], ]
#chi2_lim_e2 = [known_levs[int(chi2_lim2[0]),0], known_levs[int(chi2_lim2[1]),0]]
ax[0].axvspan(known_levs[int(chi2_lim[0]),0], known_levs[int(chi2_lim[1]),0], alpha=0.5, color='red')
ax[0].axvspan(known_levs[int(chi2_lim2[0]),0], known_levs[int(chi2_lim2[1]),0], alpha=0.5, color='red')

energies = np.linspace(0, 20, 1000)
axs[1].plot(energies, GLO(energies, Sb127_par[3], Sb127_par[0], Sb127_par[1], Sb127_par[2]), 'c-', label = 'Sb127 GLO')
axs[0].set_title('NLD')
axs[1].set_title(r'$\gamma$SF')  
axs[0].set_ylabel(r'NLD [MeV$^{-1}$]')
axs[1].set_ylabel(r'$\gamma$SF [MeV$^{-3}$]')
axs[0].set_ylim(10**ranges[0][0],10**ranges[0][1])
axs[1].set_ylim(10**ranges[1][0],1e-5)
for ax in axs:
    ax.grid()
    ax.legend()
fig.show()

#plot gsf-value distribution for one energy point
if plot_single_energy_hist:
    energies2 = gsfs[4].energies #all energy vectors are alike. Pick the 5th, but any would do
    yvals = []
    inv_chis = []
    fig2, ax = plt.subplots()
    for line in data_matrices[1]:
        if line[0]==energies2[enbin]:
            yvals.append(line[1])
            inv_chis.append(line[2])
    ax.plot([10**x for x in yvals],inv_chis)
    ax.set_title('prob distr at %s MeV'%"{:.2f}".format(energies2[enbin]))
    ax.set_xlim([4e-8,4e-7])
    ax.set_xscale('log')
    fig2.show()
            