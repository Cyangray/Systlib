#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 19:09:37 2022

@author: francesco updated 1st February 2022

Import data from calculated best fit gsf, and try to fit for the GDR with a GLO plus one Gaussian for the pygmy.
Add a fit with three gaussians. Plot the M1 strength as a SLO, extrapolated in Tin-M1-fit.py from tin isotopes.
Plot together. Calculate the TRK sum rule.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from systlib import gauss, GLO, upbend_draw, PDR_to_MeVbarn, SLO_simple, load_known_gsf

#import list 
best_fits = np.load('data/generated/best_fits.npy', allow_pickle = True)
best_gsf = best_fits[1]
best_gsf.clean_nans()
best_gsf.delete_point(-1)
Sb127_par = np.loadtxt('data/generated/127Sb_params')

#initialize common stuff
x_values_cont = np.linspace(0, 20, 1000)

#124Sn - M1
Sn124_M1 = load_known_gsf(124, 50, lab='d', nature = 'M1')

'''
Output root
Total fit funtion: 
 FCN=248.513 FROM HESSE     STATUS=OK             73 CALLS        2523 TOTAL
                     EDM=8.34445e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           1.50000e+01   5.60080e-02   8.71171e-04** at limit **
   2  p1           6.39970e+00   9.11941e-02   6.54523e-04** at limit **
   3  p2           3.11377e+02   8.47656e+00   1.22562e-03** at limit **
   4  p3           3.07468e-01   1.08669e-02   7.95917e-06   8.26944e-03
   5  p4           1.64175e-07   8.96521e-09   6.20476e-12  -3.88976e+03
   6  p5           6.99015e-01   2.67845e-02   1.76525e-05   1.86288e-02
   7  p6           6.52077e+00   4.71600e-02   2.42893e-05  -1.09472e-02
   8  p7           5.48108e-08   1.35308e-08   2.30148e-12  -4.36609e+04
   9  p8           1.13566e+00   1.45359e-01   3.81926e-06   2.20421e-03
  10  p9           9.30250e+00     fixed    
  11  p10          3.29840e+00     fixed    
  12  p11          5.48340e+00     fixed    
  
par_array_totalfit[12] = {E_r1,Gamma_r1,sigma_r1,temp,scaling1,stdev1,E_pyg1,constant1_upb,constant2_upb,E_m1,gamma_m1,sigma_m1};
'''


#Plot one Gaussian PDR and 1 GLO GDR plot
#Values from root fit
#Gaussian PDR
E1 = 6.52077
E1_unc = 4.71600e-02
C1 = 1.64175e-07
C1_unc = 8.96521e-09
Sigma1 = 6.99015e-01
Sigma1_unc = 2.67845e-02
#GLO
T = 3.07468e-01
E2 = 15.00
Gamma2 = 6.3997
Sigma2 = 3.11377e2
#Upbend
a_up3 = -1.13566
C3 = 5.48108e-08

#load M1-strength
Sn127_M1_par = np.loadtxt('data/generated/127Sn_M1_params')
res = Sn127_M1_par[:,0]

#prepare functions
x_values_cont = np.linspace(0,20,1000)
Gauss1 = gauss(x_values_cont, T, E1, C1, Sigma1)
GDR = GLO(x_values_cont, T, E2, Gamma2, Sigma2)
upbend1 = upbend_draw(x_values_cont, C3, a_up3)
M1_strength = SLO_simple(x_values_cont, res[0], res[1], res[2])
funcs_sum = [Gauss1[i] + GDR[i] + upbend1[i] + M1_strength[i] for i in range(len(x_values_cont))]

#Plot functions
cmap = matplotlib.cm.get_cmap('YlGnBu')
fig1, ax = plt.subplots(nrows = 1, ncols = 2, sharey=True)
ax[0].plot(x_values_cont, Gauss1, color = cmap(2/3), linestyle = '-.', label="Gaussian")
ax[0].plot(x_values_cont, GDR, color = 'r', linestyle = '--', label="GDR")
ax[0].plot(x_values_cont, upbend1, color = 'y', linestyle = ':', label="Upbend")
ax[0].plot(x_values_cont, M1_strength, 'r:', label="M1")
ax[0].plot(x_values_cont, funcs_sum, color = 'k', linestyle = '-', label="Total fit")
#Sn124_M1.plot(ax=ax[0], color = 'g')


ax[0].text(0.9, 0.05, 'a)', fontsize='medium', verticalalignment='center', fontfamily='serif', transform = ax[0].transAxes)
best_gsf.plot(ax = ax[0], color = 'b', style = '.', alpha = 0.7, label = 'Oslo data')
#ax[0].grid()
ax[0].legend()
#ax[0].set_title(r'$\chi^2$-score = 248.51')
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


'''
 Total fit funtion: 
 FCN=155.044 FROM HESSE     STATUS=OK            166 CALLS        3656 TOTAL
                     EDM=8.8415e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           1.50000e+01   5.21893e-01   9.06208e-03** at limit **
   2  p1           6.39970e+00   1.27220e+00   6.30436e-03** at limit **
   3  p2           3.11377e+02   4.73058e+01   1.19063e-02** at limit **
   4  p3           1.95381e-01   6.61197e-02   5.88448e-05  -1.10695e-01
   5  p4           2.18090e-09   4.08686e-10   1.18320e-12  -2.90394e+06
   6  p5           1.92935e-01   2.65329e-02   4.41428e-05  -4.74228e-03
   7  p6           3.88342e+00   2.42882e-02   1.42021e-04   1.20457e-03
   8  p7           1.56570e-07   1.54766e-08   2.18575e-11  -1.56437e+05
   9  p8           6.88135e-01   4.37870e-02   6.89848e-05  -5.09858e-02
  10  p9           6.41351e+00   7.13760e-02   8.76235e-05   4.11245e-02
  11  p10          4.61689e-08   1.54126e-08   4.95794e-11  -2.09494e+04
  12  p11          2.80511e-01   8.99204e-02   1.38009e-04  -7.44727e-04
  13  p12          7.52234e+00   7.68390e-02   3.87887e-04   1.35112e-03
  14  p13          3.10500e-08   5.53959e-09   1.89571e-04  -2.34154e-02
  15  p14          6.48115e-01   1.55938e-01   7.32284e-05   7.29318e-02
  16  p15          9.30250e+00     fixed    
  17  p16          3.29840e+00     fixed    
  18  p17          5.48340e+00     fixed    
 The chi^2 value: 155.044
 Degrees of freedom: 28
 Reduced chi^2 = 5.5373

double par_array_totalfit[18] = {E_r1,Gamma_r1,sigma_r1,temp,C_1,sigma_1,E0_1,C_2,sigma_2,E0_2,C_3,sigma_3,E0_3,constant1_upb,constant2_upb,E_m1,gamma_m1,sigma_m1};
    // In terminal output:        p0     p1       p2      p3  p4   p5      p6  p7    p8     p9  p10   p11   p12     p13            p14
'''


#Plot one gaussian "scissor", two Gaussian PDR and 1 GLO GDR plot
#Values from root fit
#manual_fit
T = 1.95381e-01
#scissors?
E1 = 3.88342
E1_unc = 2.42882e-02
C1 = 2.18090e-09
C1_unc = 4.08686e-10
Sigma1 = 1.92935e-01
Sigma1_unc = 2.65329e-02
#First gaussian
E2 = 6.41351
E2_unc = 7.13760e-02
C2 = 1.56570e-07
C2_unc = 1.54766e-08
Sigma2 = 6.88135e-01
Sigma2_unc = 4.37870e-02
#small gaussian
E3 = 7.52234
E3_unc = 7.68390e-02
C3 = 4.61689e-08
C3_unc = 1.54126e-08
Sigma3 = 2.80511e-01
Sigma3_unc = 8.99204e-02
#GLO
E4 = 15.000
Gamma4 = 6.39970e+00
Sigma4 = 3.11377e+02
#Upbend
C5 = 3.10500e-08
a_up5 = -6.48115e-01

#prepare functions
Gauss1 = gauss(x_values_cont, T, E1, C1, Sigma1)
Gauss2 = gauss(x_values_cont, T, E2, C2, Sigma2)
Gauss3 = gauss(x_values_cont, T, E3, C3, Sigma3)
GDR = GLO(x_values_cont, T, E4, Gamma4, Sigma4)
upbend1 = upbend_draw(x_values_cont, C5, a_up5)
funcs_sum = [Gauss1[i] + Gauss2[i] + Gauss3[i] + GDR[i] + upbend1[i] + M1_strength[i] for i in range(len(x_values_cont))]

#Plot functions
ax[1].plot(x_values_cont, Gauss1, color = cmap(1/3), linestyle = '-.', label="Gaussian 1")
ax[1].plot(x_values_cont, Gauss2, color = cmap(2/3), linestyle = '-.', label="Gaussian 2")
ax[1].plot(x_values_cont, Gauss3, color = cmap(3/3), linestyle = '-.', label="Gaussian 3")
ax[1].plot(x_values_cont, GDR, color = 'r', linestyle = '--')#, label="GDR")
ax[1].plot(x_values_cont, upbend1, color = 'y', linestyle = ':')#, label="Upbend")
ax[1].plot(x_values_cont, SLO_simple(x_values_cont, res[0], res[1], res[2]), 'r:')#, label="M1")
#Sn124_M1.plot(ax=ax[1], color = 'g')
ax[1].plot(x_values_cont, funcs_sum, color = 'k', linestyle = '-')#, label="Total fit")
best_gsf.plot(ax = ax[1], color = 'b', style = '.', alpha = 0.7)#, label = 'Oslo data')
ax[1].text(0.9, 0.05, 'b)', fontsize='medium', verticalalignment='center', fontfamily='serif', transform = ax[1].transAxes)
ax[1].legend(loc = 'upper left')#ncol = 2)
#ax[1].set_title(r'$\chi^2$-score = 155.044')
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