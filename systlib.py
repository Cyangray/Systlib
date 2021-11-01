#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 15:48:24 2021

@author: francesco
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/francesco/Documents/Talys-calculations')
from readlib import *
import lmfit
from sklearn.linear_model import LinearRegression


#convertion factor?
k_C = 1.439996 #e^2*fm*MeV
hc = 197.3269804 #MeV*fm
dBdE2f = 16*np.pi/27*hc**(-3)*k_C
sqrt2pi = np.sqrt(2*np.pi)
const = 1/(3*np.pi**2*hc**2*10) #mb

def sigma2f(sigma, E_g):
    return const*sigma/E_g

def SLO(E, T, E0, Gamma0, sigma0):
	# Special Lorentzian, adapted from Kopecky & Uhl (1989) eq. (2.1)
	funct = const * sigma0 * E * Gamma0**2 / ( (E**2 - E0**2)**2 + E**2 * Gamma0**2 )
	return funct

def GLO(E, T, E0, Gamma0, sigma0):
    # General Lorentzian, adapted from Kopecky & Uhl (1989) eq. (2.4)
    Gamma = Gamma0 * (E**2 + 4* np.pi**2 * T**2) / E0**2
    param1 = (E*Gamma)/( (E**2 - E0**2)**2 + E**2 * Gamma**2 )
    param2 = 0.7*Gamma0*4*np.pi**2 *T**2 /E0**5
    funct = const * (param1 + param2)*sigma0*Gamma0
    return funct

def gauss(E, T, E0, C, sigma):
    return C*np.exp(-(E-E0)**2/(2*sigma**2))/(sqrt2pi*sigma)

def f_fit_total(par, E, functions):
    """ make_fit currently supports up to 5 GLO/SLOs. To add more, continue the pattern. This 
    was, suprisingly, the best solution I found atm working with the scipy-syntax. """
    EX=[]
    Gamma=[]
    sigma=[]
    for i,function in enumerate(functions):
        if function==gauss:
            Ename = 'E' + str(i+1)
            Cname = 'C' + str(i+1)
            Sname = 'Sigma' + str(i+1)
            EX.append(par[Ename])
            Gamma.append(par[Cname])
            sigma.append(par[Sname])
            
        elif function==GLO:
            Ename = 'E' + str(i+1)
            Gname = 'Gamma' + str(i+1)
            Sname = 'Sigma' + str(i+1)
            EX.append(par[Ename])
            Gamma.append(par[Gname])
            sigma.append(par[Sname])
            
    T = par["T"]
    EX = np.array(EX)
    Gamma = np.array(Gamma)
    sigma = np.array(sigma)
    '''
    EX = np.array([
        par["E1"], 
        par["E2"],
        par["E3"]
        ])
    Gamma = np.array([
                      par["C1"], 
                      par["C2"],
                      par["Gamma3"]
                      ])
    sigma = np.array([
                      par["Sigma1"], 
                      par["Sigma2"],
                      par["Sigma3"]
                      ])
    '''
    # Calculating the sum of the singular GLOs
    output = np.zeros(len(E))
    for i in range(len(functions)):
        output += functions[i](E, T, EX[i], Gamma[i], sigma[i])
    return output

def f_residuals(parameters, nucleus, functions):
    E = nucleus.energies
    return (f_fit_total(parameters, E, functions) - nucleus.y)**2/(nucleus.yerr**2)

def chisquared(theo,exp,err,DoF=0,method = 'linear'):
    if len(theo) == len(exp) and len(exp) == len(err):
        chi2=0
        if method == 'linear':
            for i in range(len(theo)):
                chi2+=((theo[i] - exp[i])/err[i])**2
            return chi2/(len(theo)-DoF)
        elif method == 'log':
            for i in range(len(theo)):
                #chi2+=((np.log(theo[i]) - exp[i])/np.log(err[i]))**2
                expmax = exp[i] + err[i]/2
                expmin = exp[i] - err[i]/2
                chi2+=(np.log(theo[i]/exp[i])/np.log(expmax/expmin))**2
                #((np.log(theo[i]) - exp[i])/np.log(err[i]))**2
            return chi2/(len(theo)-DoF)
        else:
            print('Method not known')
    else:
        print('Lenths of arrays not matching')

def import_ocl(path,a0,a1,fermi=False):
    'import data generated with the oslo method software, and convert it to something legible'
    raw_matrix = np.loadtxt(path)
    channels = 61
    if fermi:
        polished_matrix = np.zeros((782,2))
        limit = 782
    else:
        polished_matrix = np.zeros((channels,3))
        limit=channels
    for i, el in enumerate(raw_matrix):
        if i<limit:
            polished_matrix[i,0] = a0 + a1*i
            if el == 0:
                polished_matrix[i,1] = np.nan
            else:
                polished_matrix[i,1] = el
        else:
            if el == 0:
                polished_matrix[i-channels,2] = np.nan
            else:
                polished_matrix[i-channels,2] = el
    return polished_matrix
    #return polished_matrix[~np.isnan(polished_matrix).any(axis=1)]

def normal_dist_weigh(mean, sigma, value):
    #gives the ratio between a gaussian evaluated at the mean, and at a specific value given the standard deviation
    return np.exp(-0.5*(1-((value - mean)/sigma)**2))

def calc_data_matrix(lst, rng, bins, dens):
    '''
    Function transforming a list of gsfs or nlds to a matrix containing the data
    for drawing the histogram graph, and the rgba_colors matrix to pass to a plt.scatter
    function to give the histogram varying transparency.

    Parameters
    ----------
    lst : list, ndarray
        a list containing either nld or gsf objects
    rng : tuple, list of two elements
        range for the histogram
    bins : int
        Number of bins to divide the range into
    dens : boolean
        False: distribution is normalized. True: each energy histogram is normalized
        in such a way that the highest bin equals one

    Returns
    -------
    data_matrix : ndarray
        Matrix. 1st column: energies (middle bin). 2nd column: gsf value. 3rd
        column: associated inverse chi2 value
    rgba_colors : ndarray
        First three columns: r,g,b values for a colour. 4th column: normalized inverse
        chi2 value
    '''
    
    energies = lst[4].energies #all energy vectors are alike. Pick the 5th, but any would do
    energy_array = np.array([])
    point_array = np.array([])
    inv_chi2_array = np.array([])
    #for every energy, make a histogram with numpy. Normalize each histogram so that
    for i, energy in enumerate(energies):
        y_stripe = []
        energy_stripe = []
        inv_chi2_stripe = []
        for graph in lst:
            if not np.isnan(graph.y[i]):
                y_stripe.append(graph.y[i])
                inv_chi2_stripe.append(1/graph.chi2)
        y_stripe = np.log10(y_stripe)
        hist, bin_edges = np.histogram(y_stripe, weights = inv_chi2_stripe, bins = bins, range = rng, density = dens)
        energy_stripe = np.array([energy for x in range(bins)])
        stripe_binned = [bin_edges[i] + (bin_edges[1]-bin_edges[0])/2 for i in range(hist.size)]
        if not dens:
            if hist.max() == 0:
                pass
            else:    
                hist = hist/np.nanmax(hist)
        point_array = np.concatenate([point_array, stripe_binned])
        inv_chi2_array = np.concatenate([inv_chi2_array, hist])
        energy_array = np.concatenate([energy_array, energy_stripe])
    if dens:
        inv_chi2_array = inv_chi2_array/np.nanmax(inv_chi2_array)
        
        
    data_matrix = np.stack((energy_array,point_array,inv_chi2_array)).T
    data_matrix = data_matrix[~np.isnan(data_matrix).any(axis=1)]
    data_matrix = data_matrix[~np.any(data_matrix == 0, axis=1)]
    rgba_colors = np.zeros((data_matrix[:,0].size,4))
    r,g,b,a = to_rgba('b')
    rgba_colors[:,0] = r
    rgba_colors[:,1] = g
    rgba_colors[:,2] = b
    rgba_colors[:,3] = data_matrix[:,2]
    return data_matrix, rgba_colors


class gsf:
    def __init__(self, path, label='', energycol = 0, xscol = 1, errcol = None, is_sigma = True, a0 = 0, a1 = 0, is_ocl = False):
        if is_ocl:
            self.rawmat = import_ocl(path,a0,a1)
            errcol = 2
        else:
            self.rawmat = np.loadtxt(path)
        self.path = path
        self.energies = self.rawmat[:,energycol]
        self.label = label
        if is_sigma:
            self.y = np.array([sigma2f(self.rawmat[i,xscol], self.energies[i]) for i in range(len(self.energies))])
        else:
            self.y = self.rawmat[:,xscol]
        
        if isinstance(errcol, int):
            self.error = True
            if is_sigma:
                self.yerr = np.array([sigma2f(self.rawmat[i,errcol], self.energies[i]) for i in range(len(self.energies))])
            else:
                self.yerr = self.rawmat[:,errcol]
            self.yerrplot = self.yerr
        elif isinstance(errcol, tuple) or isinstance(errcol, list):
            self.error = True
            if is_sigma:
                self.yerrdown = np.array([sigma2f(self.rawmat[i,errcol[0]], self.energies[i]) for i in range(len(self.energies))])
                self.yerrup = np.array([sigma2f(self.rawmat[i,errcol[1]], self.energies[i]) for i in range(len(self.energies))])
            else:
                self.yerrdown = self.rawmat[:,errcol[0]]
                self.yerrup = self.rawmat[:,errcol[1]]
            self.yerr = self.yerrup + self.yerrdown
            self.yerrplot = [self.yerrdown, self.yerrup]
        else:
            self.error = False
            
    
    def plot(self, ax = 0, alpha = 1, ploterrors = True, color = '', style = '.'):
        if ax:
            if self.error and ploterrors:
                ax.errorbar(self.energies, self.y, fmt = color+style, yerr = self.yerrplot, alpha = alpha, label = self.label)
            else:
                ax.plot(self.energies, self.y, color+style, alpha = alpha, label = self.label)
                
                
        else:
            if self.error and ploterrors:
                plt.errorbar(self.energies, self.y, fmt = color+style, yerr = self.yerrplot, alpha = alpha, label = self.label)
            else:
                plt.plot(self.energies, self.y, color+style, alpha = alpha, label = self.label)
        
class nld:
    def __init__(self, path, label='', energycol = 0, ldcol = 1, errcol = 100, a0 = 0, a1 = 0, is_ocl = False):
        if is_ocl:
            self.rawmat = import_ocl(path,a0,a1)
            errcol = 2
        else:
            self.rawmat = np.loadtxt(path)
        self.path = path
        self.energies = self.rawmat[:,energycol]
        self.label = label
        self.y = self.rawmat[:,ldcol]
        try:
            self.yerr = self.rawmat[:,errcol]
            self.error = True
        except:
            self.error = False
    
    def plot(self, ax = 0, alpha = 1, ploterrors = True, color = '', style = 'o'):
        if ax:
            if self.error and ploterrors:
                ax.errorbar(self.energies, self.y, yerr = self.yerr, alpha = alpha, label = self.label)
            else:
                ax.plot(self.energies, self.y, color+style, alpha = alpha, label = self.label)
        else:
            if self.error and ploterrors:
                plt.errorbar(self.energies, self.y, yerr = self.yerr, alpha = alpha, label = self.label)
            else:
                plt.plot(self.energies, self.y, color+style, alpha = alpha, label = self.label)


def load_known_gsf(A,Z,lab=''):
    '''
    function to load the gsf of a nucleus
    TODO: get path as input and override the had-coded one

    Parameters
    ----------
    A : int
        mass number
    Z : int,str
        element number, or element name
    lab : str, optional
        'oslo' if data from Oslo, 'darmstadt' if data from Darmstadt. 'o' and 'd' also work. The default is ''.

    Returns
    -------
    gsf object
        gsf object containing the loaded data

    '''
    Oslo = False
    Darmstadt = False
    nucleus = str(int(A)) + Z2Name(Z)
    if (lab=='oslo') or (lab=='Oslo') or (lab=='o') or (lab=='ocl') or (lab=='OCL'):
        nucleus = nucleus + '_o'
        Oslo = True
    elif (lab=='darmstadt') or (lab=='Darmstadt') or (lab=='d'):
        nucleus = nucleus +'_d'
        Darmstadt = True
    
    if nucleus == '129I':
        return gsf('129I-gn.txt', label = 'I129', energycol = 1, xscol = 0)
    elif (Z2Name(Z) == 'Sn') and (A in (112,114,116,118,120,124)) and Darmstadt:
        return gsf('../Tin/Darmstadt/' + str(int(A)) + 'Sn_Total_GSF_Darmstadt.dat', label = nucleus, energycol = 0, xscol = 1, errcol = 2, is_sigma = False)
    elif (Z2Name(Z) == 'Sn') and (A in (120,124)) and Oslo:
        return gsf('../Tin/Oslo/' + str(int(A)) + 'Sn_GSF.txt', label = 'Sn' + str(int(A)) + '_o', energycol = 0, xscol = 1, errcol = [3,2], is_sigma = False)
    elif (Z2Name(Z) == 'Te') and (A == 128):
        Te128_n = gsf('128Te-gn3.txt', label = 'Te128_n', energycol = 3, xscol = 0, errcol = 1)
        Te128_2n = gsf('128Te-gn3-2n.txt', label = 'Te128_2n', energycol = 3, xscol = 0, errcol = 1)
        Te128 = gsf('128Te-gn3-2n.txt', label='Te128', energycol = 3, xscol = 0, errcol = 1)
        #unite Te128_n and Te128_2n into Te128
        for i, energy in enumerate(Te128_2n.energies):
            index = np.where(Te128_n.energies == energy)[0]
            if len(index) == 0:
                Te128.energies[i] = np.nan
                Te128.y[i] = np.nan
                Te128.yerr[i] = np.nan
            else:
                Te128.energies[i] = energy
                Te128.y[i] += Te128_n.y[index[0]]
                Te128.yerr[i] += Te128_n.yerr[index[0]]
        Te128.energies = Te128.energies[~np.isnan(Te128.energies)]
        Te128.y = Te128.y[~np.isnan(Te128.y)]
        Te128.yerr = Te128.yerr[~np.isnan(Te128.yerr)]
        
        #for E<Emin for the new set, put Te128_n
        Te128_Emin = min(Te128.energies)
        lowenergy_E = []
        lowenergy_err = []
        lowenergy_gsf = []
        for i, energy in enumerate(Te128_n.energies):
            if energy < Te128_Emin:
                lowenergy_E.append(energy)
                lowenergy_err.append(Te128_n.yerr[i])
                lowenergy_gsf.append(Te128_n.y[i])
            else:
                break
        
        Te128.energies = np.delete(np.concatenate((np.array(lowenergy_E),Te128.energies)), [0,1,2], 0) #delete first three rows
        Te128.y = np.delete(np.concatenate((np.array(lowenergy_gsf), Te128.y)), [0,1,2], 0) #delete first three rows
        Te128.yerr = np.delete(np.concatenate((np.array(lowenergy_err),Te128.yerr)), [0,1,2], 0) #delete first three rows
        Te128.yerrplot = Te128.yerr
        return Te128

def plotgraph(L1n, L2n, chi2, Sb127_nld_ocl, Sb127_nld_fermi, Sb127_nld_lvl, new_rho, chi2_lim_e, nldalpha = 1,fermistart = 37):
    #plot nld
    fig,ax = plt.subplots()
    ax.plot(Sb127_nld_ocl[:,0],Sb127_nld_ocl[:,1],'b.',label='Sb127_nld', alpha = nldalpha)
    ax.plot(Sb127_nld_fermi[fermistart:,0],Sb127_nld_fermi[fermistart:,1],'k--',alpha=0.5,label='CT model')
    ax.plot(Sb127_nld_lvl[:,0],Sb127_nld_lvl[:,1],'k-',label='Known levels')
    ax.errorbar(Sb127_nld_ocl[:,0],Sb127_nld_ocl[:,1],yerr=Sb127_nld_ocl[:,2],color='b', capsize=2, ls='none')
    ax.plot(Sn,new_rho,'go',label='rho at Sn')
    ax.errorbar(Sn,new_rho,yerr=rho_Sn_err, color='g',capsize = 2, ls='none', alpha = nldalpha)
    ax.set_yscale('log')
    if annotate and not plotcloud:
        for index, point in enumerate(Sb127_nld_ocl):
            ax.annotate(str(index), xy = (point[0],point[1]))#, xytext=(row['Sn']+40, row['rho_EB']))
    ax.vlines(x=chi2_lim_e, ymin=0, ymax=new_rho*5, colors='purple', ls='--', alpha = 0.5, label='Fitting interval')
    ax.vlines(x=[Sb127_nld_ocl[L1n,0], Sb127_nld_ocl[L2n,0]], ymin=0, ymax=new_rho*5, colors='red', ls='-.', alpha = 0.5, label='L1 and L2')
    ax.set_ylim(1e-1,new_rho*5)
    ax.set_xlim(-0.5,Sn+0.5)
    ax.set_title('NLD for L1=%d, L2=%d. Chi2=%.2f'%(L1n,L2n,chi2))
    ax.grid()
    ax.legend()
    fig.show()

def plotgraph_cloud(L1n, L2n, Sb127_nld_ocl, ax, Sb127_nld_fermi = None, nldalpha = 1, fermistart = 37, errorbars = False):
    #plot nld
    ax.plot(Sb127_nld_ocl[:,0],Sb127_nld_ocl[:,1],'b.', alpha = nldalpha)
    try:
        ax.plot(Sb127_nld_fermi[fermistart:,0],Sb127_nld_fermi[fermistart:,1],'k--',alpha=nldalpha)
    except:
        pass    
    if errorbars:
        ax.errorbar(Sb127_nld_ocl[:,0],Sb127_nld_ocl[:,1],yerr=Sb127_nld_ocl[:,2],color='b', capsize=2, ls='none', alpha = nldalpha)
        
def import_Anorm_alpha(path):
    string = np.genfromtxt(path)
    Anorm, alpha = [x for x in string if np.isnan(x) == False]
    return Anorm, alpha
    #A =  1.9479 and alpha = 1.5250
    
def import_Bnorm(path):
    return np.genfromtxt(path, skip_header = 3)
    
def rho2D(rho,target_spin,spin_cutoff):
    #calculate D from rho at Bn (Larsen et al. 2011, eq. (20))
    factor = 2*spin_cutoff**2/((target_spin + 1)*np.exp(-(target_spin + 1)**2/(2*spin_cutoff**2)) + target_spin*np.exp(-target_spin**2/(2*spin_cutoff**2)))
    D = factor/rho
    return D
    
    
    