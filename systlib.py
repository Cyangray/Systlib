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
import math
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/francesco/Documents/Talys-calculations')
from readlib import *
import lmfit
from sklearn.linear_model import LinearRegression
from scipy.interpolate import interp1d


#convertion factor?
k_C = 1.439996 #e^2*fm*MeV
hc = 197.3269804 #MeV*fm
dBdE2f = 16*np.pi/27*hc**(-3)*k_C
sqrt2pi = np.sqrt(2*np.pi)
const = 1/(3*np.pi**2*hc**2*10) #mb
k_B = 8.617333262145e1 #Boltzmann konstant, keV*GK^-1

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

def upbend(E, T, a_up, C, sigma):
    return C*np.exp(a_up*E)

def upbend_draw(E,C,a_up):
    return C*np.exp(a_up*E)

def PDR_to_MeVbarn(C,sigma,E0):
    #assumes PDR described by Gaussian
    #analytical result of 1/const * int_0^\inf dE E C 1/(sqrt(2*pi)) exp(-(E-E0)^2/(2*sigma^2))
    return (C*(sigma/np.sqrt(2) + E0))/const #mb*MeV

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
        
        elif function==upbend:
            aname = 'a_up' + str(i+1)
            Cname = 'C' + str(i+1)
            EX.append(par[aname])
            Gamma.append(par[Cname])
            sigma.append(par['T'])
            
            
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

def chisquared(theo,exp,err,DoF=1,method = 'linear',reduced=True):
    if len(theo) == len(exp) and len(exp) == len(err):
        chi2=0
        if method == 'linear':
            for i in range(len(theo)):
                chi2+=((theo[i] - exp[i])/err[i])**2
            if reduced:
                return chi2/(len(theo)-DoF)
            else:
                return chi2
        elif method == 'log':
            for i in range(len(theo)):
                #chi2+=((np.log(theo[i]) - exp[i])/np.log(err[i]))**2
                expmax = exp[i] + err[i]/2
                expmin = exp[i] - err[i]/2
                chi2+=(np.log(theo[i]/exp[i])/np.log(expmax/expmin))**2
                #((np.log(theo[i]) - exp[i])/np.log(err[i]))**2
            if reduced:
                return chi2/(len(theo)-DoF)
            else:
                return chi2
        else:
            print('Method not known')
    else:
        print('Lenths of arrays not matching')

def ToLatex(nucname):
    ''' function translating a string of the form 'NNNXX' indicating the
    name of a nucleus, into something like '$^{NNN}$XX' for LaTeX rendering'''
    nums = ''
    letters = ''
    for char in nucname:
        if char.isnumeric():
            nums += char
        else:
            letters += char
    newstring = '$^{' + nums + '}$' + letters
    return newstring    

def import_ocl(path,a0,a1,channels = 61, fermi=False):
    '''import data generated with the oslo method software, and convert it to something legible'''
    raw_matrix = np.loadtxt(path)
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

def flat_distr_chi2_fade(upperlim, lowerlim, sigma, value):
    #Calculates the chi2 score for a distribution that is max and flat between
    #lowerlim and upperlim, and fading as a normal distribution with standard
    #deviation like sigma, outside these values
    if value < lowerlim:
        return ((lowerlim - value)/sigma)**2
    elif value <= upperlim:
        return 0
    elif value > upperlim:
        return ((upperlim - value)/sigma)**2
    
    
'''
def gaussian(mean,sigma,val):
    #simple gaussian function
    return np.exp(-0.5*((val-mean)/sigma)**2)/(sigma*sqrt2pi)

def normal_dist_weigh(mean, sigma, value):
    #gives the ratio between a gaussian evaluated at the mean, and at a specific value given the standard deviation
    #i.e. the function gives 1 at value=mean, and increases as an inverse gaussian away from the mean
    return np.exp(0.5*((value - mean)/sigma)**2)

def flat_distr_fade(upperlim, lowerlim, sigma, value):
    #Flat distribution, only that instead of suddently jumping from 0 to 1, it fades
    #up to the max value as a gaussian function with given sigma
    if value < lowerlim:
        return normal_dist_weigh(lowerlim,sigma,value)
    elif value <= upperlim:
        return 1.#normal_dist_weigh(value,sigma,value)
    elif value > upperlim:
        return normal_dist_weigh(upperlim,sigma,value)

    

def calc_data_matrix(lst, rng, bins, dens):
    """
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
    """
    
    energies = lst[4].energies #all energy vectors are alike. Pick the 5th, but any would do
    energy_array = np.array([])
    point_array = np.array([])
    inv_chi2_array = np.array([])
    avg_values_array = np.array([])
    upper_unc = np.array([])
    lower_unc = np.array([])
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

def calc_errors(lst, rng, limit_points, median = True):
    """
    Function taking a list of either nld or gsf objects, and returning a matrix
    containing the energy, the median value of the list of values associated to that energy,
    and the value of the gsf or the nld at a certain probability threshold in limit_points
    e.g.: limit_points is [0.05,0.50,0.95], the output matrix value_matrix contain
    the energy in the first column, the median in the second, the lower limit for
    the 90% confidence interval in the third, the median (again) in the fourth,
    and the upper limit for the 90% confidence interval in the fifth.
    if median = False, the wheighted average will be calculated and put in the
    first column, and not the median

    Parameters
    ----------
    lst : list, ndarray
        a list containing either nld or gsf objects
    rng : list, tuple
        the range where the histogram of the probability distribution will be
        calculated
    limit_points : list, ndarray
        the list of thresholds we want to locate in the distribution
    median (optional) : Boolean
        if we want to calculate the median value, and not the mean

    Returns
    -------
    value_matrix : ndarray
        matrix with the energy as first column,
        wheighted average or median as second column
        lower error as third column
        upper error as fourth column
    
    """
    bins = 500
    energies = lst[4].energies #all energy vectors are alike. Pick the 5th, but any would do
    value_matrix = np.zeros((energies.size,(limit_points.size+2)))
    
    #for every energy, make a histogram with numpy.
    for i, energy in enumerate(energies):
        y_stripe = []
        inv_chi2_stripe = []
        for graph in lst:
            if not np.isnan(graph.y[i]):
                y_stripe.append(graph.y[i])
                inv_chi2_stripe.append(1/graph.chi2)
        y_stripe = np.log10(y_stripe)
        
        hist, bin_edges = np.histogram(y_stripe, weights = inv_chi2_stripe, bins = bins, range = rng, density = True)
        binsize = bin_edges[1]-bin_edges[0]
        stripe_binned = [bin_edges[i] + binsize/2 for i in range(hist.size)]
        
        #check when the cumulative sum of the hists gets bigger than the values
        #in limit_points, and store the bin value in the uncertainties array
        
        uncertainties = np.full_like(limit_points,np.nan)#np.zeros_like(limit_points)
        norm_hist = hist*binsize
        cs = np.cumsum(norm_hist)
        if not np.all(np.isnan(cs)):
            for j, quant in enumerate(limit_points):
                if quant < 1.00000000e+00:
                    a = np.where(cs>quant,cs,np.full_like(cs,np.nan))
                    d = np.nanargmin( a)
                else:
                    try:
                        a = np.where(cs>=quant,stripe_binned,np.full_like(stripe_binned,np.nan))
                        d = np.nanargmin(a)
                    except:
                        a = np.where(cs>=quant*0.999,stripe_binned,np.full_like(stripe_binned,np.nan))
                        d = np.nanargmin(a)
                    #d = np.nanargmin(a)
                uncertainties[j] = 10**stripe_binned[d]
            
        """
        running_sum = 0
        j=0
        for stripe_binned_point,hist_point in zip(stripe_binned,hist):
            running_sum += hist_point*binsize
            if running_sum > limit_points[j]:
                uncertainties[j] = 10**(stripe_binned_point-binsize/2)
                j += 1
            if j == uncertainties.size:
                break
        """
        if not median:
            #calculate the weighted average
            avg = np.average(y_stripe, weights = inv_chi2_stripe)
        else:
            median_index = np.where(limit_points == 0.5)[0][0]
            avg = uncertainties[median_index]
            
        
        row = np.hstack([energy,avg,uncertainties])
        value_matrix[i,:] = row 
        
    value_matrix = value_matrix[~np.isnan(value_matrix).any(axis=1)]
    value_matrix = value_matrix[~np.any(value_matrix == 0, axis=1)]
    return value_matrix

def calc_errors_astrorate(lst, rng, limit_points, value = 'ncrate', median = True):
    """
    Function taking a ncrate object, and returning a matrix
    containing the temperature, the median value of the list of values associated to that temperature,
    and the value of the ncrate at a certain probability threshold in limit_points
    e.g.: limit_points is [0.05,0.50,0.95], the output matrix value_matrix contain
    the temperature in GK in the first column, the median in the second, the lower limit for
    the 90% confidence interval in the third, the median (again) in the fourth,
    and the upper limit for the 90% confidence interval in the fifth.
    if median = False, the wheighted average will be calculated and put in the
    first column, and not the median

    Parameters
    ----------
    lst : list, ndarray
        a list containing ncrate objects
    rng : list, tuple
        the range where the histogram of the probability distribution will be
        calculated
    limit_points : list, ndarray
        the list of thresholds we want to locate in the distribution
    median (optional) : Boolean
        if we want to calculate the median value, and not the mean

    Returns
    -------
    value_matrix : ndarray
        matrix with the temperature as first column,
        wheighted average or median as second column
        lower error as third column
        upper error as fourth column
    """
    
    bins = 500
    Temperatures = lst[4].T #all temperature vectors are alike. Pick the 5th, but any would do
    value_matrix = np.zeros((Temperatures.size,(limit_points.size+2)))
    
    #for every energy, make a histogram with numpy.
    for i, T in enumerate(Temperatures):
        y_stripe = []
        inv_chi2_stripe = []
        
        for graph in lst:
            if value == 'ncrate':
                y = graph.ncrate
            else:
                y = graph.MACS
            if not np.isnan(y[i]):
                y_stripe.append(y[i])
                inv_chi2_stripe.append(1/graph.chi2)
        y_stripe = np.log10(y_stripe)
        
        hist, bin_edges = np.histogram(y_stripe, weights = inv_chi2_stripe, bins = bins, range = rng, density = True)
        binsize = bin_edges[1]-bin_edges[0]
        stripe_binned = [bin_edges[i] + binsize/2 for i in range(hist.size)]
        
        #check when the cumulative sum of the hists gets bigger than the values
        #in limit_points, and store the bin value in the uncertainties array
        
        uncertainties = np.full_like(limit_points,np.nan)#np.zeros_like(limit_points)
        norm_hist = hist*binsize
        cs = np.cumsum(norm_hist)
        if not np.all(np.isnan(cs)):
            for j, quant in enumerate(limit_points):
                if quant < 1.00000000e+00:
                    a = np.where(cs>quant,cs,np.full_like(cs,np.nan))
                    d = np.nanargmin( a)
                else:
                    try:
                        a = np.where(cs>=quant,stripe_binned,np.full_like(stripe_binned,np.nan))
                        d = np.nanargmin(a)
                    except:
                        a = np.where(cs>=quant*0.999,stripe_binned,np.full_like(stripe_binned,np.nan))
                        d = np.nanargmin(a)
                    #d = np.nanargmin(a)
                uncertainties[j] = 10**stripe_binned[d]
                
        if not median:
            #calculate the weighted average
            avg = np.average(y_stripe, weights = inv_chi2_stripe)
        else:
            median_index = np.where(limit_points == 0.5)[0][0]
            avg = uncertainties[median_index]
            
        row = np.hstack([T,avg,uncertainties])
        value_matrix[i,:] = row 
        
    value_matrix = value_matrix[~np.isnan(value_matrix).any(axis=1)]
    value_matrix = value_matrix[~np.any(value_matrix == 0, axis=1)]
    return value_matrix

def calc_errors_talys(nlds, limit_points):
    """
    Function taking a list of nld objects, and calculating the median of all the
    values in talys_nld_cnt.txt, and their associated values at confidences in limit_points
    in a similar way of calc_errors. It returns a list of len(limit_points) matrices,
    each corresponding to the nld talys extrapolation at the respective confidence
    in limit_points

    Parameters
    ----------
    nlds : list, ndarray
        a list containing nld objects
    rng : list, tuple
        the range where the histogram of the probability distribution will be
        calculated
    limit_points : list, ndarray
        the list of thresholds we want to locate in the distribution

    Returns
    -------
    talys_matrices : list
        list of matrices, each in the form of talys_nld_cnt.txt, expressing the
        nld at each limit_point (so the one in the middle, corresponding to 
        limit_point = 0.5, will be the matrix corresponding to the median nld
        calculated in calc_errors)
    
    """
    bins = 500
    energies = nlds[4].talys[:,0] #all energy vectors are alike. Pick the 5th, but any would do
    line_length = nlds[4].talys[0,:].size #length of a row
    talys_matrices = [np.zeros_like(nlds[4].talys) for i in limit_points]
    for i, energy in enumerate(energies):
        print(energy)
        for talys_matrix in talys_matrices:
            talys_matrix[i,0] = energy
        
        for j in range(1,line_length):
            y_stripe = []
            inv_chi2_stripe = []
            for curr_nld in nlds:
                if not np.isnan(curr_nld.talys[i,j]):
                    y_stripe.append(curr_nld.talys[i,j])
                    inv_chi2_stripe.append(1/curr_nld.chi2)
            #y_stripe = np.log10(y_stripe)
            hist, bin_edges = np.histogram(y_stripe, weights = inv_chi2_stripe, bins = bins, density = True)
            binsize = bin_edges[1]-bin_edges[0]
            stripe_binned = [bin_edges[i] + binsize/2 for i in range(hist.size)]
            
            #check when the cumulative sum of the hists gets bigger than the values
            #in limit_points, and store the bin value in the uncertainties array
            uncertainties = np.zeros_like(limit_points)
            running_sum = 0
            k=0
            for stripe_binned_point,hist_point in zip(stripe_binned,hist):
                running_sum += hist_point*binsize
                if running_sum >= limit_points[k]:
                    uncertainties[k] = stripe_binned_point #10**stripe_binned_point
                    k += 1
                if k == uncertainties.size:
                    break
            
            #median_index = np.where(limit_points == 0.5)[0][0]
            #avg = uncertainties[median_index]
            for talys_matrix, unc in zip(talys_matrices, uncertainties):
                talys_matrix[i,j] = unc

    return talys_matrices

def print_talys_matrices(talys_matrices, limit_points):
    for talys_matrix, limits in zip(talys_matrices,limit_points):
        filename = 'nld_talys_' + '{:.2f}'.format(limits) + '.txt'
        fmt = "%7.2f %6.3f %9.2E %8.2E %8.2E" + 30*" %8.2E"
        header = "U[MeV]  T[MeV]  NCUMUL   RHOOBS   RHOTOT     J=0      J=1      J=2      J=3      J=4      J=5      J=6      J=7      J=8      J=9     J=10     J=11     J=12     J=13     J=14     J=15     J=16     J=17     J=18     J=19     J=20     J=21     J=22     J=23     J=24     J=25     J=26     J=27     J=28     J=29"
        np.savetxt(filename, talys_matrix, fmt=fmt, header=header)

'''

class gsf:
    def __init__(self, path, label='', energycol = 0, xscol = 1, errcol = None, is_sigma = True, a0 = 0, a1 = 0, is_ocl = False):
        if is_ocl:
            self.rawmat = import_ocl(path,a0,a1)
            errcol = 2
        else:
            self.rawmat = np.loadtxt(path)
        self.path = path
        self.energies = self.x = self.rawmat[:,energycol]
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
            
    def clean_nans(self):
        clean_E = []
        clean_y = []
        clean_yerr = []
        for E, y, yerr in zip(self.energies, self.y, self.yerr):
            if not math.isnan(y):
                clean_E.append(E)
                clean_y.append(y)
                clean_yerr.append(yerr)
        self.energies = np.array(clean_E)
        self.y = np.array(clean_y)
        self.yerr = np.array(clean_yerr)
        self.yerrplot = self.yerr
    
    def delete_point(self, position):
        self.y = np.delete(self.y, position)
        self.energies = np.delete(self.energies, position)
        self.yerr = np.delete(self.yerr, position)
        self.yerrplot = np.delete(self.yerrplot, position)
        
        '''
        if self.error:
            lists = [self.y, self.yerr, self.energies, self.yerrplot]
        else:
            lists = [self.y, self.energies]
        for el in lists:
            el = np.delete(el, position)
        '''
        
        
        
    def plot(self, ax = 0, alpha = 1, ploterrors = True, color = '', style = '.', label = None):
        if label is None:
            label = self.label
            
        if ax:
            if self.error and ploterrors:
                ax.errorbar(self.energies, self.y, fmt = color+style, yerr = self.yerrplot, alpha = alpha, label = label)
            else:
                ax.plot(self.energies, self.y, color+style, alpha = alpha, label = label)
                
                
        else:
            if self.error and ploterrors:
                plt.errorbar(self.energies, self.y, fmt = color+style, yerr = self.yerrplot, alpha = alpha, label = label)
            else:
                plt.plot(self.energies, self.y, color+style, alpha = alpha, label = label)
        
class nld:
    def __init__(self, path=None, label='', energycol = 0, ldcol = 1, errcol = 100, a0 = 0, a1 = 0, is_ocl = False):
        if path:
            if is_ocl:
                self.rawmat = import_ocl(path,a0,a1)
                errcol = 2
            else:
                self.rawmat = np.loadtxt(path)
            self.path = path
            self.energies = self.x = self.rawmat[:,energycol]
            self.label = label
            self.y = self.rawmat[:,ldcol]
            try:
                self.yerr = self.rawmat[:,errcol]
                self.error = True
            except:
                self.error = False
    
    def clean_nans(self):
        clean_E = []
        clean_y = []
        clean_yerr = []
        for E, y, yerr in zip(self.energies, self.y, self.yerr):
            if not y == np.NaN:
                clean_E.append(E)
                clean_y.append(y)
                clean_yerr.append(yerr)
        self.energies = np.array(clean_E)
        self.y = np.array(clean_y)
        self.yerr = np.array(clean_yerr)
    
    def plot(self, ax = 0, alpha = 1, ploterrors = True, color = '', style = 'o', label = None):
        if label is None:
            label = self.label
            
        if ax:
            if self.error and ploterrors:
                ax.errorbar(self.energies, self.y, yerr = self.yerr, alpha = alpha, label = label)
            else:
                ax.plot(self.energies, self.y, color+style, alpha = alpha, label = label)
        else:
            if self.error and ploterrors:
                plt.errorbar(self.energies, self.y, yerr = self.yerr, alpha = alpha, label = label)
            else:
                plt.plot(self.energies, self.y, color+style, alpha = alpha, label = label)

class astrorate:
    def __init__(self, path=None, label=''):
        if path:
            self.path = path
            self.ncrate_mat = readastro_path(path)
            self.T = self.x = self.ncrate_mat[:,0]
            self.ncrate = self.y = self.ncrate_mat[:,1]
            self.MACS = self.ncrate_mat[:,2]
        self.label = label
    
    def plot(self, ycol = 'ncrate', ax = 0, alpha = 1, color = '', style = 'o'):
        if ycol=='ncrate':
            y = self.ncrate
        elif ycol == 'MACS':
            y = self.MACS
        if ax:
            ax.plot(self.T, y, color+style, alpha = alpha, label = self.label)
        else:
            plt.plot(self.T, y, color+style, alpha = alpha, label = self.label)


def load_known_gsf(A,Z,lab=''):
    '''
    function to load the gsf of a nucleus
    TODO: get path as input and override the hard-coded one

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
        return gsf('data/nuclear/Tin/Darmstadt/' + str(int(A)) + 'Sn_Total_GSF_Darmstadt.dat', label = nucleus, energycol = 0, xscol = 1, errcol = 2, is_sigma = False)
    elif (Z2Name(Z) == 'Sn') and (A in (120,124)) and Oslo:
        return gsf('data/nuclear/Tin/Oslo/' + str(int(A)) + 'Sn_GSF.txt', label = 'Sn' + str(int(A)) + '_o', energycol = 0, xscol = 1, errcol = [3,2], is_sigma = False)
    elif (Z2Name(Z) == 'Te') and (A == 128):
        Te128_n = gsf('data/nuclear/128Te-gn3.txt', label = 'Te128_n', energycol = 3, xscol = 0, errcol = 1)
        Te128_2n = gsf('data/nuclear/128Te-gn3-2n.txt', label = 'Te128_2n', energycol = 3, xscol = 0, errcol = 1)
        Te128 = gsf('data/nuclear/128Te-gn3-2n.txt', label='Te128', energycol = 3, xscol = 0, errcol = 1)
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

'''
def plotgraph_cloud(L1n, L2n, Sb127_nld_ocl, ax, Sb127_nld_fermi = None, nldalpha = 1, fermistart = 37, errorbars = False):
    #plot nld
    ax.plot(Sb127_nld_ocl[:,0],Sb127_nld_ocl[:,1],'b.', alpha = nldalpha)
    try:
        ax.plot(Sb127_nld_fermi[fermistart:,0],Sb127_nld_fermi[fermistart:,1],'k--',alpha=nldalpha)
    except:
        pass    
    if errorbars:
        ax.errorbar(Sb127_nld_ocl[:,0],Sb127_nld_ocl[:,1],yerr=Sb127_nld_ocl[:,2],color='b', capsize=2, ls='none', alpha = nldalpha)
'''
        
def import_Anorm_alpha(path):
    string = np.genfromtxt(path)
    Anorm, alpha = [x for x in string if np.isnan(x) == False]
    return Anorm, alpha
    #A =  1.9479 and alpha = 1.5250
    
def import_Bnorm(path):
    #return np.genfromtxt(path, skip_header = 3)[0,0]
    arr = np.loadtxt(path, skiprows = 3)
    return arr.item()
    
def rho2D(rho,target_spin,spin_cutoff):
    #calculate D from rho at Bn (Larsen et al. 2011, eq. (20))
    factor = 2*spin_cutoff**2/((target_spin + 1)*np.exp(-(target_spin + 1)**2/(2*spin_cutoff**2)) + target_spin*np.exp(-target_spin**2/(2*spin_cutoff**2)))
    D = factor/rho
    return D
    
def log_interp1d(xx, yy, **kwargs):
    """ Interpolate a 1-D function.logarithmically """
    logy = np.log(yy)
    lin_interp = interp1d(xx, logy, kind='linear', **kwargs)
    log_interp = lambda zz: np.exp(lin_interp(zz))
    return log_interp

'''
def gen_GDR(min_energy = 0.1):
    #import GLO parameters for Te and Sn
    Te128_par = np.loadtxt('data/nuclear/128Te_params_root')
    Sn126_par = np.loadtxt('data/nuclear/126Sn_params')
    
    #Calculate GLO parameters of 127Sb as the average of Te and Sn
    Sb127_par = (Te128_par + Sn126_par)/2
    
    energies = np.arange(min_energy, 30.1, 0.1)
    GDR = GLO(energies, Sb127_par[3], Sb127_par[0], Sb127_par[1], Sb127_par[2])
    return energies, GDR
'''

def make_E1_M1_files(gsf_folder_path, Sn, A, Z, a0, a1, M1_frac = 0.1, filename = 'strength.nrm', target_folder = None, high_energy_interp=None):
    '''
    Function that takes the path of a Oslo Method generated gsf with energy in the first column and 
    the gsf in the second, and writes two tables for both E1 and M1 ready to be taken
    as input by TALYS.
    '''
    
    # read/write gsf files
    if gsf_folder_path != '':
        if gsf_folder_path[-1] != '/':
            gsf_folder_path = gsf_folder_path + '/'
        gsf_path = gsf_folder_path + filename
    else:
        gsf_path = filename
        
    gsf = import_ocl(gsf_path, a0, a1)
    gsf = gsf[~np.isnan(gsf).any(axis=1)]
    gsf = np.delete(gsf, [-1,-2], 0) #delete last two points
    gsf = gsf[:,:-1]
    if high_energy_interp is not None:
        gsf = np.vstack((gsf,high_energy_interp))
    
    if target_folder is not None:
        if target_folder != '':
            if target_folder[-1] != '/':
                target_folder = target_folder + '/'
        gsf_folder_path = target_folder
    
    fn_gsf_outE1 = gsf_folder_path + "gsfE1.dat"
    fn_gsf_outM1 = gsf_folder_path + "gsfM1.dat"
    
    # The file is/should be writen in [MeV] [MeV^-3] [MeV^-3]
    if gsf[0, 0] == 0:
        gsf = gsf[1:, :]
    Egsf = gsf[:, 0]
    gsfE1 = gsf[:, 1]*(1-M1_frac)
    gsfM1 = gsf[:, 1]*M1_frac

    # REMEMBER that the TALYS functions are given in mb/MeV (Goriely's tables)
    # so we must convert it (simple factor)
    factor_from_mb = 8.6737E-08   # const. factor in mb^(-1) MeV^(-2)
    
    fE1 = log_interp1d(Egsf, gsfE1, fill_value="extrapolate")
    fM1 = log_interp1d(Egsf, gsfM1, fill_value="extrapolate")
    
    Egsf_out = np.arange(0.1, 30.1, 0.1)
    
    header = f" Z=  {Z} A=  {A}\n" + "  U[MeV]  fE1[mb/MeV]"
    # gsfE1 /= factor_from_mb
    np.savetxt(fn_gsf_outE1, np.c_[Egsf_out, fE1(Egsf_out)/factor_from_mb],
               fmt="%9.3f%12.3E", header=header)
    # gsfM1 /= factor_from_mb
    np.savetxt(fn_gsf_outM1, np.c_[Egsf_out, fM1(Egsf_out)/factor_from_mb],
               fmt="%9.3f%12.3E", header=header)
    return Egsf_out, fE1(Egsf_out)/factor_from_mb
    

def make_TALYS_tab_file(talys_nld_path, ocl_nld_path, A, Z):
    '''
    Function that incorporates the talys_nld_cnt.txt produced by counting, into
    the big Zz.tab file from TALYS. The code overwrites the Zz.tab file, so be careful
    '''
    newfile_content = ''
    if Z < 10:
        Zstring = '  ' + str(Z)
    elif Z < 100:
        Zstring = ' ' + str(Z)
    else:
        Zstring = str(Z)
        
    if A < 10:
        Astring = '  ' + str(A)
    elif A < 100:
        Astring = ' ' + str(A)
    else:
        Astring = str(A)
        
    isotope_strip = 'Z='+ Zstring +' A=' + Astring
    isotopeline = 100000
    with open(talys_nld_path, 'r') as read_obj:
        with open(ocl_nld_path, 'r') as ocl_nld_f:
            ocl_nld = ocl_nld_f.readlines()
            for n, line in enumerate(read_obj):
                stripped_line = line
                new_line = stripped_line
                if isotope_strip in line:
                    isotopeline = n
                    
                if n >= isotopeline + 3:
                    if (n - (isotopeline + 3)) < len(ocl_nld):
                        new_line = ocl_nld[n - (isotopeline + 3)]
            
                newfile_content += new_line
    
    #print(newfile_content)
    
    with open(talys_nld_path, 'w') as write_obj:
        write_obj.write(newfile_content)
    
        
def find_chis(vals,chis):
    #function taking as input all chi2-scores associated to the values for a single energy or temperature
    #it finds where the function crosses the chi2+1 line
    whole_mat = np.vstack((vals,chis)).T
    chimin = np.min(chis)
    #lower_mat = whole_mat[np.argwhere(chis<=(chimin+1)),:]
    lower_mat = whole_mat[chis<=(chimin+1)]
    upper_mat = whole_mat[(chis>(chimin+1)) & (chis<(chimin+3))]
    
    min1 = lower_mat[lower_mat[:,0]==np.min(lower_mat[:,0])][0]
    min2 = upper_mat[upper_mat[:,0]==np.min(upper_mat[:,0])][0]
    max1 = lower_mat[lower_mat[:,0]==np.max(lower_mat[:,0])][0]
    max2 = upper_mat[upper_mat[:,0]==np.max(upper_mat[:,0])][0]
    
    # y(x) = A + Bx
    Bmin = (min2[1]-min1[1])/(min2[0]-min1[0])
    Amin = min2[1]-Bmin*min2[0]
    Bmax = (max2[1]-max1[1])/(max2[0]-max1[0])
    Amax = max2[1]-Bmax*max2[0]
    
    # evaluate at Y = chimin + 1
    Y = chimin + 1
    Xmin = (Y-Amin)/Bmin
    Xmax = (Y-Amax)/Bmax
    
    return[Xmin,Xmax]

def calc_errors_chis(lst):
    
    xx = lst[4].x #all energy or temperature vectors are alike. Pick the 5th, but any would do
    val_matrix = np.zeros((xx.size,6))
    for i, x in enumerate(xx):
        chis = []
        vals = []
        row = np.zeros(6)
        counterflag = 0
        for graph in lst:
            if not np.isnan(graph.y[i]):
                chis.append(graph.chi2)
                vals.append(graph.y[i])
                #it may be that all y[i] are nans. Do a check so that the code only saves data if there actually are values to analyze
                counterflag += 1
        if counterflag > 10:
            best_fit = vals[np.argwhere(chis==np.min(chis)).item()]
            errmin, errmax = find_chis(vals,chis)
            row[:] = [x, best_fit, -best_fit+2*errmin, errmin, errmax, 2*errmax-best_fit]
        else:
            row[:] = [x, np.nan, np.nan, np.nan, np.nan, np.nan]
        val_matrix[i,:] = row[:]
    return val_matrix

def calc_errors_chis_MACS(lst):
    
    xx = lst[4].x #all energy or temperature vectors are alike. Pick the 5th, but any would do
    val_matrix = np.zeros((xx.size,6))
    for i, x in enumerate(xx):
        chis = []
        vals = []
        row = np.zeros(6)
        counterflag = 0
        for graph in lst:
            if not np.isnan(graph.MACS[i]):
                chis.append(graph.chi2)
                vals.append(graph.MACS[i])
                #it may be that all y[i] are nans. Do a check so that the code only saves data if there actually are values to analyze
                counterflag += 1
        if counterflag > 10:
            best_fit = vals[np.argwhere(chis==np.min(chis)).item()]
            errmin, errmax = find_chis(vals,chis)
            row[:] = [x, best_fit, -best_fit+2*errmin, errmin, errmax, 2*errmax-best_fit]
        else:
            row[:] = [x, np.nan, np.nan, np.nan, np.nan, np.nan]
        val_matrix[i,:] = row[:]
    return val_matrix
        