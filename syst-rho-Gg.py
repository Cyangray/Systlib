#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Fri Aug 20 11:52:45 2021

@author: francesco, last modified 1st February 2022

Code to find the level density at separation energy and the <Gamma_gamma> value 
for 127Sb by comparing to neighbouring nuclei.
for rho(Sn): 
Use the experimental values for neighbouring stable nuclei present in Mughabghab
and RIPL, and fit the theoretical values to these from E&B2005.
then use the fit in order to guess the value for 127Sb.
for <Gamma_gamma>:
Use experimental values again from Mughabghab and/or RIPL and extrapolate these
values to find the predicted one for 127Sb by either taking the average of the
values, or doing a linear extrapolation by 
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as sc
from sklearn.linear_model import LinearRegression
from systlib import ToLatex
from readlib import Z2Name

#what to fit? 'rho' or 'Gg'
What_fit = 'Gg'
#if 'Gg': choose what to have on the x axis
Gg_axis = 'A'
#if 'Gg': calculate average? 'Mugh', 'RIPL' or 'both'
calc_avg = 'both'
#if 'Gg': regression? 'Mugh', 'RIPL' or 'both'
Gg_regr = 'Mugh'
#Sort by element? 50, 51, 52. False or 0 if not
element = 0
#Subtract dn from Sn? False or 0 if not
subtract_dn = 0
#Plot?
plot = 1
#if 'rho': fit?
fit127Sb = True
#if 'rho': if so, fit with RIPL or Mugh?
fit_with='Mugh'




#import resonance data from RIPL
RIPL = np.genfromtxt('data/nuclear/resonances0.dat', skip_header=125, max_rows=16)
#Matrix with the neutron pairing energy, to be subtracted from the separation 
#energy of neutron-even nuclei in order to be able to compare their NLD at Sn 
#to the neutron-odd nuclei
dn_matr = np.loadtxt('data/nuclear/doba_delta_n.dat')

#create dataframe with nuclear data.
#MughD0 + uncertainty are taken from Mughabghab. RhoMugh and uncertainty calculated with d2rho
#RhoRIPL + uncertainty are are calculated with d2rho from D0 values
#rhoEB are the Egidy & Bucurescu values calculated with Robin (RMI(FG), AFG, RMI reduction factor = 1)
#Spincutoff calculated with Robin
elnamesRIPL =        ['116Sn',   '117Sn', '118Sn',     '119Sn', '120Sn',                '122Sn',             '124Sn',              '121Sb',              '123Sb',                          '122Te',   '123Te',   '124Te', '125Te',     '126Te',              '128Te',              '130Te']
MughD0 =             [np.NaN,    61,        700,       95,        1485,      np.NaN,    2965,      np.NaN,    5935,      np.NaN,    10.4,      np.NaN,    24,        np.NaN,    np.NaN,    131.7,     22.3,      220,       40,        617,       np.NaN,    2250,      np.NaN,    3390,     np.NaN]
MughD0uncertainty =  [np.NaN,    7,         100,       14,        130,       np.NaN,    145,       np.NaN,    650,       np.NaN,    1.1,       np.NaN,    1.7,       np.NaN,    np.NaN,    13.5,      2.8,       30,        2.8,       170,       np.NaN,    656,       np.NaN,    540,      np.NaN]
Sn =                 [9563.5,    6943.1,    9326,      6483.5,    9104,      6170,      8815,      5946,      8489,      5733,      9254,      6806,      8960,      6468,      8380,      9840,      6929,      9424.48,   6569,      9113.69,   6228,      8783,      6082,      8419.5,   5929]
rhoEB =              [0.78017E6, 0.17588E6, 0.71936E6, 0.12671E6, 0.56166E6, 0.84778E5, 0.32540E6, 0.46971E5, 0.14148E6, 0.20386E5, 0.58657E7, 0.83427E6, 0.29103E7, 0.39089E6, 0.37572E6, 0.30539E7, 0.49992E6, 0.16622E7, 0.26359E6, 0.81737E6, 0.12012E6, 0.31496E6, 0.45946E5, 0.9481E5, 0.14535E5]
rhoRIPL =            [np.NaN,    1.407E+05, 5.601E+05, 9.005E+04, 3.860E+05, 5.136E+04, np.NaN,    2.245E+04, np.NaN,    1.393E+04, np.NaN,    1.035E+06, np.NaN,    4.718E+05, np.NaN,    np.NaN,    4.663E+05, 1.498E+06, 3.145E+05, 8.984E+05, 1.115E+05, np.NaN,    4.958E+04, np.NaN,   5.349E+04]
rhoRIPLUncertainty = [np.NaN,    1.563E+04, 6.427E+04, 1.930E+03, 8.579E+04, 8.218E+03, np.NaN,    1.138E+03, np.NaN,    2.787E+03, np.NaN,    1.593E+05, np.NaN,    5.897E+04, np.NaN,    np.NaN,    4.471E+04, 1.671E+05, 4.289E+04, 6.268E+04, 2.962E+04, np.NaN,    1.215E+04, np.NaN,   1.783E+04]
rhoMugh=             [np.NaN,    1.248E+05, 5.601E+05, 9.005E+04, 3.657E+05, 4.323E+04, np.NaN,    2.241E+04, np.NaN,    1.174E+04, np.NaN,    1.294E+06, np.NaN,    4.718E+05, np.NaN,    np.NaN,    5.169E+05, 1.686E+06, 3.145E+05, 9.658E+05, 1.157E+05, np.NaN,    3.327E+04, np.NaN,   2.367E+04]
rhoMughUncertainty = [np.NaN,    1.230E+04, 6.427E+04, 1.286E+04, 5.390E+04, 3.785E+03, np.NaN,    1.096E+03, np.NaN,    1.286E+03, np.NaN,    1.369E+05, np.NaN,    3.342E+04, np.NaN,    np.NaN,    5.299E+04, 2.117E+05, 4.289E+04, 6.760E+04, 3.187E+04, np.NaN,    9.701E+03, np.NaN,   3.770E+03]
Zs =                 [50,        50,        50,        50,        50,        50,        50,        50,        50,        50,        51,        51,        51,        51,        51,        52,        52,        52,        52,        52,        52,        52,        52,        52,       52]
As =                 [116,       117,       118,       119,       120,       121,       122,       123,       124,       125,       121,       122,       123,       124,       127,       122,       123,       124,       125,       126,       127,       128,       129,       130,      131]
elnames =            ['116Sn',   '117Sn',   '118Sn',   '119Sn',   '120Sn',   '121Sn',   '122Sn',   '123Sn',   '124Sn',   '125Sn',   '121Sb',   '122Sb',   '123Sb',   '124Sb',   '127Sb',   '122Te',   '123Te',   '124Te',   '125Te',   '126Te',   '127Te',   '128Te',   '129Te',   '130Te',  '131Te']
spincutoff =         [5.760,     5.581,     5.769,     5.569,     5.819,     5.621,     5.894,     5.720,     5.999,     5.859,     6.133,     5.932,     6.196,     5.994,     6.445,     6.020,     5.791,     6.060,     5.839,     6.144,     5.931,     6.263,     6.077,     6.434,    6.294]
GgMugh =             [53,        117,       45,        92,        36,        np.NaN,    49,        np.NaN,    50,        np.NaN,    90,        np.NaN,    97.4,      np.NaN,    np.NaN,    73.4,      122,       59,        106,       53,        np.NaN,    37,        np.NaN,    50,       np.NaN]
GgMughUncertainty =  [3,         20,        7,         14,        3,         np.NaN,    5,         np.NaN,    5,         np.NaN,    5,         np.NaN,    13,        np.NaN,    np.NaN,    4.3,       8,         6,         6,         15,        np.NaN,    5,         np.NaN,    15,       np.NaN]
dns = []
for i in range(len(Zs)):
    found=False
    for line in dn_matr:
        if line[0]==As[i] and line[1]==Zs[i]:
            dns.append(line[3])
            found=True
            break
    if found==False:
        dns.append(0)

d = {"Z": pd.Series(Zs, index=elnames),
     "A": pd.Series(As, index=elnames),
     "Sn": pd.Series(Sn, index=elnames),
     "D0_RIPL": pd.Series(RIPL[:,5]*1000, index=elnamesRIPL),
     "D0_RIPL_err": pd.Series(RIPL[:,6]*1000, index=elnamesRIPL),
     "D0_Mugh": pd.Series(MughD0, index=elnames),
     "D0_Mugh_err": pd.Series(MughD0uncertainty, index=elnames),
     "rho_EB": pd.Series(rhoEB, index=elnames),
     "rho_Mugh": pd.Series(rhoMugh, index=elnames),
     "rho_Mugh_err": pd.Series(rhoMughUncertainty, index=elnames),
     "rho_RIPL": pd.Series(rhoRIPL, index=elnames),
     "rho_RIPL_err": pd.Series(rhoRIPLUncertainty, index=elnames),
     "dn": pd.Series(dns, index=elnames),
     "Gg_RIPL": pd.Series(RIPL[:,9], index=elnamesRIPL),
     "Gg_RIPL_err": pd.Series(RIPL[:,10], index=elnamesRIPL),
     "Gg_Mugh": pd.Series(GgMugh, index=elnames),
     "Gg_Mugh_err": pd.Series(GgMughUncertainty, index=elnames)
     }

df_whole = pd.DataFrame(d)
df_whole.sort_values(by=['Z','A'],inplace=True)

Sb127_serie = df_whole.loc['127Sb',:].copy(deep=True)
if subtract_dn:
    df_whole['Sn'] -= df_whole['dn']
if element:
    df = df_whole.loc[df_whole['Z'] == element]
else:
    df = df_whole

if What_fit == 'rho':
    #Follows the algorithm from Kullmann et al. 2019 for chi squared minimization. 
    #Fit the theoretical values to the experimental (up to a scale factor b).
    #From the fit, guess the theoretical prediction of 127Sb by multiplying the E&B value of rho, with b
    if fit127Sb:
        def chi2sum(b):
            chi2=0
            for index, row in df.iterrows():
                if (not np.isnan(row['rho_EB'])) and (not np.isnan(row['rho_' + fit_with])):
                    chi2+=(b*row['rho_EB'] - row['rho_'+ fit_with])**2/row['rho_'+ fit_with + '_err']**2
            return chi2
        
        result = sc.minimize(chi2sum, x0=0.63)
        bmin = result.x[0]
        chimin = result.fun
        if element:
            print('fitting by ' + Z2Name(element) + ' with ' + fit_with)
        print('Original, E&B value of 127Sb: %f'%(df_whole.at['127Sb','rho_EB']))
        print('b (scale factor): %s'%bmin)
        print('predicted value of rho for 127Sb: %s'%(df_whole.at['127Sb','rho_EB']*bmin))
        print('chi2 score: %f'%chimin)
    else:
        bmin = 1
        
    #Plot rho as function of separation energy
    if plot:
        df = df.append(Sb127_serie)
        df.reset_index(inplace=True)
        xaxis = df['Sn']
        fig, ax = plt.subplots()
        ax.plot(xaxis, df['rho_EB']*bmin, 'go', alpha=0.5, label='E&B')
        ax.plot(xaxis, df['rho_Mugh'], 'ro', alpha=0.5, label='Mughabghab')
        ax.errorbar(xaxis, df['rho_Mugh'], yerr=df['rho_Mugh_err'], color='r', capsize=2, ls='none')
        ax.plot(xaxis, df['rho_RIPL'], 'bo', alpha=0.5, label='RIPL')
        ax.errorbar(xaxis, df['rho_RIPL'], yerr=df['rho_RIPL_err'], color='b', capsize=2, ls='none')
        ax.set_ylabel(r'$\rho$($S_n$) [MeV$^{-1}$]')
        ax.set_xlim(np.nanmin(df.loc[:,'Sn'])-300,np.nanmax(df.loc[:,'Sn'])+500)
        plt.yscale('log')
        plt.grid()
        plt.legend()
        for index, row in df.iterrows():
            ax.annotate(r'$^{%d}$%s'%(row['A'],Z2Name(row['Z'])), xy = (row['Sn']+50, row['rho_EB']*7e-1))#, xytext=(row['Sn']+40, row['rho_EB']))
        plt.show()
        
        
        
elif What_fit == 'Gg':
    #extrapolate the Gg value og 127Sb from a model buildt on different nuclei
    if calc_avg:
        subdf = df.loc[df['A']%2 == 1]
        row_name2 = ''
        if calc_avg == 'both':
            row_name = 'Gg_Mugh'
            row_name2 = 'Gg_RIPL'
        elif calc_avg == 'Mugh' or calc_avg == 'RIPL':
            row_name = 'Gg_' + calc_avg
        else:
            print('Something went wrong')
        
        if calc_avg == 'both':
            mean = (subdf[row_name].sum(skipna=True) + subdf[row_name2].sum(skipna = True))/(subdf[row_name].count() + subdf[row_name2].count())
        else:
            mean = subdf[row_name].mean(skipna = True)
        print('calculated mean with %s: %f'%(calc_avg, mean))
        
    if Gg_regr:
        subdf = df.loc[df['A']%2 == 1]
        row_name2 = ''
        if Gg_regr == 'both':
            row_name = 'Gg_Mugh'
            row_name2 = 'Gg_RIPL'
        elif Gg_regr == 'Mugh' or Gg_regr == 'RIPL':
            row_name = 'Gg_' + Gg_regr
        else:
            print('Something went wrong')
        
        val_matr = subdf[['A',row_name]].values
        val_matr = val_matr[~np.isnan(val_matr).any(axis=1)]
        if Gg_regr == 'both':
            val_matr2 = subdf[['A',row_name2]].values
            val_matr2 = val_matr2[~np.isnan(val_matr2).any(axis=1)]
            val_matr = np.vstack((val_matr,val_matr2))
        
        X = val_matr[:,0].reshape(-1,1)
        Y = val_matr[:,1].reshape(-1,1)
        linear_regressor = LinearRegression()  # create object for the class
        linear_regressor.fit(X, Y)  # perform linear regression
        Y_pred = linear_regressor.predict([[127]])  # make predictions
        print('calculated prediction with linear regression with %s: %f'%(Gg_regr, Y_pred[0,0]))
                
    #Plot Gg as function of separation energy
    if plot:
        df.reset_index(inplace=True)
        xaxis = df[Gg_axis]
        fig, ax = plt.subplots()
        #fig.set_size_inches(20, 10)
        ax.plot(xaxis, df['Gg_Mugh'], 'bo', alpha=0.7, label='Mugh')
        ax.errorbar(xaxis, df['Gg_Mugh'], yerr=df['Gg_Mugh_err'], color='b', ecolor = 'b', capsize=2, ls='none')
        ax.plot(xaxis, df['Gg_RIPL'], 'gv', alpha=0.4, label='RIPL')
        ax.errorbar(xaxis, df['Gg_RIPL'], yerr=df['Gg_RIPL_err'], color='g', ecolor = 'g', capsize=2, ls='none')
        ax.set_ylabel(r'$\langle\Gamma_\gamma\rangle$ [meV]')
        if Gg_regr:
            Y_pred0 = linear_regressor.predict([[117]])
            ax.plot([117,127], [Y_pred0[0,0],Y_pred[0,0]],'k--')
        ax.plot(127,105,'k^', label = r'$^{127}$Sb, pred.')
        ax.errorbar(127, 105, yerr=25, color='k', capsize=2, ls='none')
        
        #plt.yscale('log')
        #plt.grid()
        plt.legend()
        if Gg_axis == 'A':
            for index, row in df.iterrows():
                ax.annotate(ToLatex(row['index']), xy = (row['A'], row['Gg_Mugh']), xytext=(row['A'] + 0.2, row['Gg_Mugh']))
                ax.annotate(r'$^{127}$Sb', xy = (127, 105), xytext=(127.2, 105))
            ax.set_xlabel('A')
            
        elif Gg_axis == 'Sn':
            for index, row in df.iterrows():
                ax.annotate(ToLatex(row['index']), xy = (row['Sn'], row['Gg_Mugh']), xytext=(row['Sn'] + 40, row['Gg_Mugh']))
                ax.annotate(r'$^{127}$Sb', xy = (8383, 105), xytext=(8423, 105))
            ax.set_xlabel(r'$S_n$ [MeV]')
        ax.set_xlim(115.5,131.9)
        plt.show()
