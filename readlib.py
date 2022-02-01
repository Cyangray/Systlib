#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:15:13 2020

@author: francesco, updated 1st February 2022

shortened version of readlib.py made for the analysis of 127Sb data
"""

import numpy as np

#constants
k_B = 8.617333262145e1 #Boltzmann konstant, keV*GK^-1
c_vacuum = 299792458 #speed of light in vacuum, m/s
N_A = 6.02214076e23 #Avogadro's number

#create useful dictionaries
Zdict = {'H' :1,
        'He' :2,
        'Li' :3,
        'Be' :4,
        'B' :5,
        'C' :6,
        'N' :7,
        'O' :8,
        'F' :9,
        'Ne' :10,
        'Na' :11,
        'Mg' :12,
        'Al' :13,
        'Si' :14,
        'P' :15,
        'S' :16,
        'Cl' :17,
        'Ar' :18,
        'K' :19,
        'Ca' :20,
        'Sc' :21,
        'Ti' :22,
        'V' :23,
        'Cr' :24,
        'Mn' :25,
        'Fe' :26,
        'Co' :27,
        'Ni' :28,
        'Cu' :29,
        'Zn' :30,
        'Ga' :31,
        'Ge' :32,
        'As' :33,
        'Se' :34,
        'Br' :35,
        'Kr' :36,
        'Rb' :37,
        'Sr' :38,
        'Y' :39,
        'Zr' :40,
        'Nb' :41,
        'Mo' :42,
        'Tc' :43,
        'Ru' :44,
        'Rh' :45,
        'Pd' :46,
        'Ag' :47,
        'Cd' :48,
        'In' :49,
        'Sn' :50,
        'Sb' :51,
        'Te' :52,
        'I' :53,
        'Xe' :54,
        'Cs' :55,
        'Ba' :56,
        'La' :57,
        'Ce' :58,
        'Pr' :59,
        'Nd' :60,
        'Pm' :61,
        'Sm' :62,
        'Eu' :63,
        'Gd' :64,
        'Tb' :65,
        'Dy' :66,
        'Ho' :67,
        'Er' :68,
        'Tm' :69,
        'Yb' :70,
        'Lu' :71,
        'Hf' :72,
        'Ta' :73,
        'W' :74,
        'Re' :75,
        'Os' :76,
        'Ir' :77,
        'Pt' :78,
        'Au' :79,
        'Hg' :80,
        'Tl' :81,
        'Pb' :82,
        'Bi' :83,
        'Po' :84,
        'At' :85,
        'Rn' :86,
        'Fr' :87,
        'Ra' :88,
        'Ac' :89,
        'Th' :90,
        'Pa' :91,
        'U' :92,
        'Np' :93,
        'Pu' :94,
        'Am' :95,
        'Cm' :96,
        'Bk' :97,
        'Cf' :98,
        'Es' :99,
        'Fm' :100,
        'Md' :101,
        'No' :102,
        'Lr' :103,
        'Rf' :104,
        'Db' :105,
        'Sg' :106,
        'Bh' :107,
        'Hs' :108,
        'Mt' :109,
        'Ds' :110,
        'Rg' :111,
        'Cn' :112,
        'Nh' :113,
        'Fl' :114,
        'Mc' :115,
        'Lv' :116,
        'Ts' :117,
        'Og' :118
        }

isodict = {v: k for k, v in Zdict.items()}

ompdict = {1 : 'jlmomp-y',
           2 : 'localomp-n'}

#function library
def Z2Name(Z):
    '''convert Z to corresponding element name, if Z is a number. Otherwise it returns the input string'''
    if isinstance(Z, float):
        Z = int(Z)
    if isinstance(Z, int):
        return isodict[Z]
    else:
        return Z

def Name2Z(Xx):
    '''convert element name to corresponding Z if input is string. Otherwise it returns the input int'''
    if isinstance(Xx, str):
        return Zdict[Xx]
    else:
        return Xx

def search_string_in_file(file_name, string_to_search):
    with open(file_name, 'r') as read_obj:
        for n, line in enumerate(read_obj):
            if string_to_search in line:
                return n

def readstrength(nucleus, A, ldmodel, massmodel, strength, omp):
    #read the strength function from the gsf.txt mashup file
    filepath = 'data/TALYS/' + ZandAtoName(A,nucleus) + 'gsf.txt'
    rowsskip = 83* ( ((ldmodel-1)*3*8*2) + ((massmodel-1)*8*2) + ((strength-1)*2) + ((omp-1)) )
    return np.loadtxt(filepath, skiprows = (rowsskip + 2), max_rows = 81)

def readastro(nucleus, A, ldmodel, massmodel, strength, omp):
    #read the astrorate from the ncrates.txt mashup file
    filepath = 'data/TALYS/' + ZandAtoName(A,nucleus) + 'ncrates.txt' 
    rowsskip = 33* ( ((ldmodel-1)*3*8*2) + ((massmodel-1)*8*2) + ((strength-1)*2) + ((omp-1)) )
    ncrates = np.loadtxt(filepath, skiprows = (rowsskip + 3), max_rows = 30)
    return np.delete(ncrates, [0,1,2,3], 0) #delete first four rows

def readreaclib(nucleus, A, reaclib_file = 'reaclib'):
    Xx = Z2Name(nucleus)
    nucleus1 = Xx.lower() + str(A)
    nucleus1 = (5 - len(nucleus1))*' ' + nucleus1
    nucleus2 = Xx.lower() + str(A + 1)
    nucleus2 = (5 - len(nucleus2))*' ' + nucleus2
    reaction = 'n' + nucleus1 + nucleus2
    try:
        rowsskip = search_string_in_file(reaclib_file, reaction) + 1
        row1 = np.genfromtxt(reaclib_file, skip_header = rowsskip, max_rows = 1, autostrip = True, delimiter = 13)
        row1 = row1[~np.isnan(row1)]
        row2 = np.genfromtxt(reaclib_file, skip_header = rowsskip + 1, max_rows = 1, autostrip = True, delimiter = 13)
        row2 = row2[~np.isnan(row2)]
        return np.concatenate((row1, row2))
    except:
        print("reaction \"" + reaction + "\" not in " + reaclib_file + "!")
        bad_a = np.zeros(7)
        bad_a[:] = np.NaN
        return bad_a
    
def nonsmoker(a, T9):
    return np.exp(a[0] + a[1]*T9**-1 + a[2]*T9**(-1/3) + a[3]*T9**(1/3) + a[4]*T9 + a[5]*T9**(5/3) + a[6]*np.log(T9))

def rate2MACS(rate,target_mass_in_au,T):
    '''
    Parameters
    ----------
    rate : ncrate, in cm^3 mol^-1 s^-1
    target_mass_in_au : self explanatory
    T : Temperature in GK

    Returns
    -------
    MACS: in mb

    '''
    
    m1 = target_mass_in_au
    mn = 1.008664915
    red_mass = (m1*mn)/(m1 + mn) * 931494.10242 #keV/c^2
    v_T = np.sqrt(2*k_B*T/red_mass)*c_vacuum
    return rate/v_T*1e25/(N_A) #mb

def MACS2rate(MACS,target_mass_in_au,T):
    '''
    Parameters
    ----------
    MACS : MACS, in mb
    target_mass_in_au : self explanatory
    T : Temperature in GK

    Returns
    -------
    rate: in cm^3 s^-1 mol^-1
    '''
    
    m1 = target_mass_in_au
    mn = 1.008664915
    red_mass = (m1*mn)/(m1 + mn) * 931494.10242 #keV/c^2
    v_T = np.sqrt(2*k_B*T/red_mass)*c_vacuum
    return MACS*v_T*1e-25*(N_A) #mb

def xs2rate(energies,xss,target_mass_in_au,Ts=None):
    '''
    Parameters
    ----------
    energies : array
        Energy array in keV
    xss : array
        cross section array in mb (must correspond to E)
    target_mass_in_au : float
        self explanatory

    Returns
    -------
    array matrix with temperature and n,gamma rate in cm^3/mole*s

    '''
    mb2cm2 = 1e-27 #b to cm^2
    m1 = target_mass_in_au
    mn = 1.008664915
    red_mass = (m1*mn)/(m1 + mn) * 931494.10242/(c_vacuum*100)**2 #keV s^2/cm^2
    if Ts is None:
        T9extrange = np.array([0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    else:
        T9extrange = Ts
    rate = []
    for T in T9extrange:
        integrand = [energies[i]*xss[i]*np.exp(-energies[i]/(k_B*T)) for i in range(len(energies))]
        integral = np.trapz(integrand,x=energies)
        curr_rate = (8/(np.pi*red_mass))**(1/2)*(k_B*T)**(-3/2)*integral*N_A*mb2cm2
        rate.append(curr_rate)
    return rate


def xs2MACS(energies,xss,Ts=None):
    '''
    Parameters
    ----------
    energies : array
        Energy array in keV
    xss : array
        cross section array in mb (must correspond to E)
    target_mass_in_au : float
        self explanatory

    Returns
    -------
    array matrix with temperature and MACS in mb

    '''
    
    if Ts is None:
        T9extrange = np.array([0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    else:
        T9extrange = Ts
    MACS = []
    for T in T9extrange:
        integrand = [energies[i]*xss[i]*np.exp(-energies[i]/(k_B*T)) for i in range(len(energies))]
        integral = np.trapz(integrand,x=energies)
        curr_MACS = 2/(np.sqrt(np.pi)*(k_B*T)**2)*integral
        MACS.append(curr_MACS)
    return MACS 

def Name2ZandA(inputstr):
    '''Translates an input string in the form Xx123 to [Z, A]'''
    
    Xx = ''
    A = ''
    for character in inputstr:
        if character.isalpha():
            Xx += character
        else:
            A += character
    return [Name2Z(Xx.title()), int(A)]

def ZandAtoName(A,Z):
    '''translates A and Z into a XxxNn isotope name string'''
    Astr = str(A)
    Zstr = Z2Name(Z)
    return Astr + Zstr
    
def translate_nuclides(path):
    '''reads the "nuclides" file from SkyNet to a matrix (list of lists) in the
    form [[Z, A], [Z, A], ...], keeping the line order of "nuclides". This can
    then be associated with the isotope abundances output from SkyNet.'''
    
    text = np.genfromtxt(path, dtype = str)
    outmatrix = []
    for line in text:
        if line == 'n':
            Z = 0
            A = 1
        elif line == 'p':
            Z = 1
            A = 1
        elif line == 'd':
            Z = 1
            A = 2
        elif line == 't':
            Z = 1
            A = 3
        else:
            Z, A = Name2ZandA(line)
        outmatrix.append([Z,A])
    return outmatrix