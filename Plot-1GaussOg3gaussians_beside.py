#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 19:09:37 2022

@author: francesco updated 1st February 2022

Import data from calculated best fit gsf, and try to fit for the GDR with a GLO plus one GLO for the pygmy
+ a fit with three gaussians. Plot together. Calculate the TRK sum rule.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from systlib import gauss, GLO, upbend_draw, PDR_to_MeVbarn

#import list 
best_fits = np.load('data/generated/best_fits.npy', allow_pickle = True)
best_gsf = best_fits[1]
best_gsf.clean_nans()
best_gsf.delete_point(-1)
Sb127_par = np.loadtxt('data/generated/127Sb_params')

#initialize common stuff
x_values_cont = np.linspace(0, 20, 1000)

#Plot one Gaussian PDR and 1 GLO GDR plot
#Values from root fit
#Gaussian PDR
E1 = 6.64092
E1_unc = 0.05
C1 = 2.03646e-7
C1_unc = 0.03e-7
Sigma1 = 7.7655e-1
Sigma1_unc = 0.11e-1
#GLO
T = 4.0472e-1
E2 = 15.00
Gamma2 = 6.3997
Sigma2 = 3.11377e2
#Upbend
a_up3 = -1.54164
C3 = 8.21051e-8

#prepare functions
Gauss1 = gauss(x_values_cont, T, E1, C1, Sigma1)
GDR = GLO(x_values_cont, T, E2, Gamma2, Sigma2)
upbend1 = upbend_draw(x_values_cont, C3, a_up3)
funcs_sum = [Gauss1[i] + GDR[i] + upbend1[i] for i in range(len(x_values_cont))]

#Plot functions
cmap = matplotlib.cm.get_cmap('YlGnBu')
fig1, ax = plt.subplots(nrows = 1, ncols = 2, sharey=True)
ax[0].plot(x_values_cont, Gauss1, color = cmap(2/3), linestyle = '-.', label="Gaussian")
ax[0].plot(x_values_cont, GDR, color = 'r', linestyle = '--', label="GDR")
ax[0].plot(x_values_cont, upbend1, color = 'y', linestyle = ':', label="Upbend")
ax[0].plot(x_values_cont, funcs_sum, color = 'k', linestyle = '-', label="Total fit")
best_gsf.plot(ax = ax[0], color = 'b', style = '.', alpha = 0.7, label = 'Oslo data')
#ax[0].grid()
ax[0].legend()
ax[0].set_title(r'$\chi^2$-score = 262.84')
ax[0].set_ylabel(r'GSF [MeV$^{-3}$]')
ax[0].set_xlabel(r'$E_\gamma$ [MeV]')
plt.yscale('log')
xlims = [-0.5,16]
ax[0].set_xlim(xlims)
plt.ylim([2e-9,1.5e-5])

#TRK
Z=51
A=127
N=A-Z
TRK_sum = 60*N*Z/A #mb*MeV
TRK_fac = PDR_to_MeVbarn(C1,Sigma1,E1)/TRK_sum
TRK_fac_max = PDR_to_MeVbarn(C1+C1_unc,Sigma1+Sigma1_unc,E1+E1_unc)/TRK_sum
TRK_fac_min = PDR_to_MeVbarn(C1-C1_unc,Sigma1-Sigma1_unc,E1-E1_unc)/TRK_sum
print('TRK fraction: %s percent'%(TRK_fac*100))
print('TRK err up: %s percent'%((TRK_fac_max-TRK_fac)*100))
print('TRK err down: %s percent'%((-TRK_fac_min+TRK_fac)*100))
print('\n')

#Plot one gaussian "scissor", two Gaussian PDR and 1 GLO GDR plot
#Values from root fit
#manual_fit
T = 0.352386
#scissors?
E1 = 3.88954e+00
E1_unc = 0.02
C1 = 2.32492e-09
C1_unc = 0.41e-9
Sigma1 = 1.99011e-01
Sigma1_unc = 0.027
#First gaussian
E2 = 6.45851e+00
E2_unc = 0.07
C2 = 1.76986e-07
C2_unc = 0.17e-7
Sigma2 = 7.23647e-01
Sigma2_unc = 0.42
#small gaussian
E3 = 7.58739e+00
E3_unc = 0.08
C3 = 6.14106e-08
C3_unc = 1.96e-8
Sigma3 = 3.21445e-01
Sigma3_unc = 0.094
#GLO
E4 = 15.000
Gamma4 = 6.39970e+00
Sigma4 = 3.11377e+02
#Upbend
C5 = 3.24958e-08
a_up5 = -8.58557e-01

#prepare functions
Gauss1 = gauss(x_values_cont, T, E1, C1, Sigma1)
Gauss2 = gauss(x_values_cont, T, E2, C2, Sigma2)
Gauss3 = gauss(x_values_cont, T, E3, C3, Sigma3)
GDR = GLO(x_values_cont, T, E4, Gamma4, Sigma4)
upbend1 = upbend_draw(x_values_cont, C5, a_up5)
funcs_sum = [Gauss1[i] + Gauss2[i] + Gauss3[i] + GDR[i] + upbend1[i] for i in range(len(x_values_cont))]

#Plot functions
ax[1].plot(x_values_cont, Gauss1, color = cmap(1/3), linestyle = '-.', label="Gaussian 1")
ax[1].plot(x_values_cont, Gauss2, color = cmap(2/3), linestyle = '-.', label="Gaussian 2")
ax[1].plot(x_values_cont, Gauss3, color = cmap(3/3), linestyle = '-.', label="Gaussian 3")
ax[1].plot(x_values_cont, GDR, color = 'r', linestyle = '--', label="GDR")
ax[1].plot(x_values_cont, upbend1, color = 'y', linestyle = ':', label="Upbend")
ax[1].plot(x_values_cont, funcs_sum, color = 'k', linestyle = '-', label="Total fit")
best_gsf.plot(ax = ax[1], color = 'b', style = '.', alpha = 0.7, label = 'Oslo data')
#ax[1].grid()
ax[1].legend()
ax[1].set_title(r'$\chi^2$-score = 157.43')
ax[1].set_yscale('log')
ax[1].set_xlabel(r'$E_\gamma$ [MeV]')
ax[1].set_xlim(xlims)
fig1.show()

#TRK
Z=51
A=127
N=A-Z
TRK_sum = 60*N*Z/A #mb*MeV
TRK_fac = (PDR_to_MeVbarn(C2,Sigma2,E2)+PDR_to_MeVbarn(C3,Sigma3,E3))/(TRK_sum)
TRK_fac_max = (PDR_to_MeVbarn(C2+C2_unc,Sigma2+Sigma2_unc,E2+E2_unc) + PDR_to_MeVbarn(C3+C3_unc,Sigma3+Sigma3_unc,E3+E3_unc))/TRK_sum
TRK_fac_min = (PDR_to_MeVbarn(C2-C2_unc,Sigma2-Sigma2_unc,E2-E2_unc) + PDR_to_MeVbarn(C3-C3_unc,Sigma3-Sigma3_unc,E3-E3_unc))/TRK_sum
print('TRK fraction: %s percent'%(TRK_fac*100))
print('TRK err up: %s percent'%((TRK_fac_max-TRK_fac)*100))
print('TRK err down: %s percent'%((TRK_fac- TRK_fac_min)*100))