#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Wed Sep 15 10:04:01 2021

@author: francesco, updated 31 January 2022

Import all data from available tin isotopes, and try to fit for the GDR with a GLO (plus two
gaussians modelling the pygmy). Import also the parameters obtained by Toft et al. (2011)
When all parameters for the different tins are obtained, fit for them
and try to predict the value of the parameters ofthe GLO fitting for 126Sn. These
will be used in order to fit for the 127Sb GDR together with the parameters for the GLO
fitting the 128Te GDR.
'''

import numpy as np
import matplotlib.pyplot as plt
import lmfit
from sklearn.linear_model import LinearRegression
from systlib import *

#Import the tins nuclei (objects already made in systlib.py) and initialize lists
As = np.array([112, 114, 116, 118, 120, 124])
nuclei = []
for A in As:
    nuclei.append(load_known_gsf(A, 50, lab='d'))

#import value from Toft et al. 2011 for the GLO parametrization of GDR for the other tins
As_Toft = np.array([117, 119, 121, 122])
#               E0      gamma   sigma   T

Sn117_params = np.array([15.66,  5.02,   254,    0.46])
Sn119_params = np.array([15.53,  4.81,   253,    0.40])
Sn121_params = np.array([15.53,  4.81,   253,    0.25])
Sn122_params = np.array([15.59,  4.77,   256,    0.25])
Toft_params = [Sn117_params, Sn119_params, Sn121_params, Sn122_params]
#initialize lists
Es = []
gammas = []
sigmas = []
Ts = []

#initialize 3x2 figure where all the plots will be drawn
fig1, axs1_u = plt.subplots(nrows = 3, ncols = 2)
fig1.set_size_inches(20, 10)
axs1 = axs1_u.flatten()
#loop for every isotope of tin
for i,nucleus in enumerate(nuclei):
    #give initial guess for the parameters to be guessed my the minimization algorithm
    #values got from fitting from Maria's Master's degree
    T = 0.33
    #First gaussian
    E1 = 6.47
    C1 = 0.83e-7
    sigma1 = 0.45
    #second gaussian
    E2 = 8.17
    C2 = 2.71e-7
    sigma2 = 0.87
    #GLO
    E3 = 15.54
    gamma3 = 5.7
    sigma3 = 270
    
    x_values_cont = np.linspace(0, 20, 1000)
    
    #minimization algorithm
    p0 = [T, #GLO temperature parameter
          E1, C1, sigma1, #first gaussian
          E2, C2, sigma2, #second gaussian
          E3, gamma3, sigma3 #GLO
          ]
    p0_functions = [
        gauss,
        gauss,
        GLO
        ]
    p0_functions_names = [
        "gauss",
        "gauss",
        "GLO"
        ]
    
    #These two coefficients give the lower and upper bound where the algorithm
    #will look for guesses. 0.7 here means 70% and 1.5 150% of the suggested value above
    coeff1 = 0.7
    coeff2 = 1.5
    minimum = np.array(p0)*coeff1
    maximum = np.array(p0)*coeff2
    
    params = lmfit.Parameters()
    params.add_many(
        ("T",      p0[0], True, minimum[0]*0.9/coeff1, maximum[0]*1.1/coeff2),
        ("E1",     p0[1], True, minimum[1], maximum[1]),
        ("C1",     p0[2], True, minimum[2], maximum[2]),
        ("Sigma1", p0[3], True, minimum[3], maximum[3]),
        ("E2",     p0[4], True, minimum[4], maximum[4]),
        ("C2",     p0[5], True, minimum[5], maximum[5]),
        ("Sigma2", p0[6], True, minimum[6], maximum[6]),
        ("E3",     p0[7], True, minimum[7], maximum[7]),
        ("Gamma3", p0[8], True, minimum[8], maximum[8]),
        ("Sigma3", p0[9], True, minimum[9], maximum[9])
        )
    
    #minimization procedure
    results = lmfit.minimize(fcn=f_residuals, params=params, args=(nucleus,p0_functions), method="least_squares")
    print("Chi square: ", results.redchi)
    chosen_params = results.params
    print(chosen_params.pretty_print())
    
    #save guessed parameters
    Es.append(chosen_params['E3'])
    gammas.append(chosen_params['Gamma3'])
    sigmas.append(chosen_params['Sigma3'])
    Ts.append(chosen_params['T'])
    
    #Plot results
    axs1[i].plot(x_values_cont, gauss(x_values_cont, chosen_params['T'], chosen_params['E1'], chosen_params['C1'], chosen_params['Sigma1']), '--', label="gauss_fit1")
    axs1[i].plot(x_values_cont, gauss(x_values_cont, chosen_params['T'], chosen_params['E2'], chosen_params['C2'], chosen_params['Sigma2']), '--', label="gauss_fit2")
    axs1[i].plot(x_values_cont, GLO(x_values_cont, chosen_params['T'], chosen_params['E3'], chosen_params['Gamma3'], chosen_params['Sigma3']), '--', label="GLO_fit3")
    axs1[i].plot(x_values_cont, f_fit_total(chosen_params, E=x_values_cont, functions = p0_functions), '-', label="Tot_fit")
    nucleus.plot(ax = axs1[i])
    axs1[i].grid()
    axs1[i].legend()
    axs1[i].set_title(nucleus.label)
    axs1[i].set_yscale('log')
    axs1[i].set_ylim([1e-9,1e-5])
    
fig1.show()

#Use the data found for the 6 tin isotopes above, to predict the parameters for Sb with a linear regression.
As = np.concatenate((As, As_Toft))
Es = Es + [nuc[0] for nuc in Toft_params]
gammas = gammas + [nuc[1] for nuc in Toft_params]
sigmas = sigmas + [nuc[2] for nuc in Toft_params]
Ts = Ts + [nuc[3] for nuc in Toft_params]
datas = [Es, gammas, sigmas, Ts]
labs = ['E0', 'Gamma', 'sigma', 'T']
colors = ['r','b','g','y']
As2 = np.append(As, 126)

#Plot the parameters as function of neutron number, and the fit
fig2, ax = plt.subplots()
coeffs = []
output = '# E, Gamma, sigma, T\n'

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
    
    ax.plot(126,datapoint,colors[i] + 'o',mfc=None)
    output += str(datapoint) + '\n'

np.savetxt('126Sn_params',output_matr,header = '# E, Gamma, sigma, T')

ax.grid()
ax.legend()
fig2.show()












