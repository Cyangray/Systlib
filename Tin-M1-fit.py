#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Wed Sep 15 10:04:01 2021

@author: francesco, updated 21 March 2022

Import all data from available tin isotopes from Darmstadt and try to fit for 
the M1 strength with a lorentzian.
When all parameters for the different tins are obtained, fit for them
and try to predict the value of the parameters ofthe GLO fitting for 126Sn. These
will be used in order to guess the M1 strength of 127Sb
'''

import numpy as np
import matplotlib.pyplot as plt
import lmfit
from sklearn.linear_model import LinearRegression
from systlib import load_known_gsf, gauss, GLO, f_fit_total, f_residuals, SLO_simple
import scipy as sc

#Import the tins nuclei (objects already made in systlib.py) and initialize lists
As = np.array([114, 116, 118, 120, 124])
nuclei = []
for A in As:
    nuclei.append(load_known_gsf(A, 50, lab='d', nature = 'M1'))

#initialize lists
Es = []
gammas = []
sigmas = []
x_values_cont = np.linspace(0, 20, 1000)

#initialize 3x2 figure where all the plots will be drawn
fig1, axs1_u = plt.subplots(nrows = 3, ncols = 2)
fig1.set_size_inches(20, 10)
axs1 = axs1_u.flatten()
#loop for every isotope of tin
for i,nucleus in enumerate(nuclei):
    #scipy optimize
    res, pcov = sc.optimize.curve_fit(SLO_simple,nucleus.energies,nucleus.y, p0=[9,2.5,4])#, bounds = ([7,1,3],[11,6,8]))
    
    #save guessed parameters
    Es.append(res[0])
    gammas.append(res[1])
    sigmas.append(res[2])
    
    #Plot results
    axs1[i].plot(x_values_cont, SLO_simple(x_values_cont, res[0], res[1], res[2]), '--', label="SLO_fit")
    nucleus.plot(ax = axs1[i])
    axs1[i].grid()
    axs1[i].legend()
    axs1[i].set_title(nucleus.label)
    axs1[i].set_yscale('log')
    axs1[i].set_ylim([1e-9,1e-5])
    
fig1.show()

#Use the data found for the 6 tin isotopes above, to predict the parameters for Sb with a linear regression.
datas = [Es, gammas, sigmas]
labs = ['E0', 'Gamma', 'sigma']
colors = ['r','b','g']
As2 = np.append(As, 127)

#Plot the parameters as function of neutron number, and the fit
fig2, ax = plt.subplots()
coeffs = []
output = '# E, Gamma, sigma\n'

output_matr = None
print('126Sn data:')
for i,data in enumerate(datas):
    regr = LinearRegression()
    regr.fit(As.reshape((-1, 1)), data)
    predictions_As = regr.predict(As.reshape(-1, 1))
    predictions_As2 = regr.predict(As2.reshape(-1, 1))
    data2 = np.append(data,predictions_As2[-1])
    noise = np.std(data2 - predictions_As2)
    noises = np.array([np.random.normal(predictions_As2,noise) for j in range(10_000)])
    l,u = np.quantile(noises, [0.025, 0.975], axis = 0)
    
    
    datapoint = predictions_As2[-1]
    stderr = u[-1]-datapoint
    print(labs[i] + ': ' + str(datapoint) + ' +- %s'%(stderr))
    
    if output_matr is None:
        output_matr = np.array([datapoint,stderr]).reshape(1,-1)
    else:
        newline = np.array([datapoint,stderr]).reshape(1,-1)
        output_matr = np.concatenate((output_matr, newline))

    ax.plot(As,data,colors[i] + 'o', label=labs[i])
    ax.plot(As2,predictions_As2,'k--',alpha=0.5)
    ax.plot(As2,l,colors[i] + '--')
    ax.plot(As2,u,colors[i] + '--')
    
    ax.plot(127,datapoint,colors[i] + 'o',mfc=None)
    output += str(datapoint) + '\n'

np.savetxt('data/generated/127Sn_M1_params',output_matr,header = '# E, Gamma, sigma')

ax.grid()
ax.legend()
fig2.show()












