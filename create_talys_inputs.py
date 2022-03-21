#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 15:26:05 2021

@author: francesco, updated 21st March 2022

Script to loop through the folders made with run_counting.py, and writes for
each gsf and nld a e1strength, m1strength and Sb.tab file to be used as input by talys,
possibly in Saga. In order to keep the number of nld-gsf pairs low, only the pairs
with a chi2 score between chimin and chimin+2 will be translated into talys-input files.
(To calculate all, just put chimin2 to a very big number)
In order to extrapolate for the GDR, the 127Sb gsf is attached
to the 128Te GDR
"""

import os
import numpy as np
from systlib import make_E1_M1_files, make_TALYS_tab_file, load_known_gsf, make_E1_M1_files_simple

a0 = -0.7599
a1 = 0.1562
Sn = 8.380
A = 127
Z = 51
blist = np.linspace(0.40,1.4,50)
base_rho = 375720
L1max = 23
Gglist = np.linspace(67.5,142.5,13)
Te128 = load_known_gsf(128,52)
master_folder = '127Sb-saga-original/'
root_folder = '/home/francesco/Documents/127Sb-experiment/Systlib/'
data_folder = root_folder + 'data/'
dataset_folder = root_folder + 'Make_dataset/'
indexes_to_delete = np.argwhere(Te128.energies<10) #delete values less than 10 MeV
Te128energies = np.delete(Te128.energies, indexes_to_delete, 0) #delete first three rows
Te128y = np.delete(Te128.y, indexes_to_delete, 0) #delete first three rows
Te128mat = np.vstack((Te128energies,Te128y)).T


nlds = np.load('data/generated/nlds.npy', allow_pickle = True)
gsfs = np.load('data/generated/gsfs.npy', allow_pickle = True)

Sn127_M1_par = np.loadtxt('data/generated/127Sn_M1_params')


#load best fits
best_fits = np.load('data/generated/best_fits.npy', allow_pickle = True)
best_gsf = best_fits[1]
best_gsf.clean_nans()
best_gsf.delete_point(-1)
best_nld = best_fits[0]
best_nld.clean_nans()
best_nld.delete_point([-1,-2,-3])
best_ncrate = best_fits[2]

#find chimin (common for nld, gsfs)
chimin = best_gsf.chi2

#calculate chimin+2
chimin2 = chimin + 2

#make new gsf list with only chi<chimin+2 (count elements!)
gsfs_filtered = []
for gsf in gsfs:
    if gsf.chi2<=chimin2:
        gsfs_filtered.append(gsf)

#begin the main loop
try:
    os.mkdir(dataset_folder + '/' + master_folder)
except:
    pass
for count, gsf in enumerate(gsfs_filtered):
    L1 = gsf.L1
    L2 = gsf.L2
    rho = gsf.rho
    Gg = gsf.Gg
    
    #retrieve the nld from 127Sb-database using L1 L2 and rho, build Sb.tab, put in folder
    bstr_int = "{:.2f}".format(gsf.b)
    bstr = bstr_int.translate(bstr_int.maketrans('','', '.'))
    new_rho_str = '{:.6f}'.format(rho)
    new_dir_rho = bstr + '-' + str(int(rho)) + '/'
    new_dir_L1_L2 = 'L1-'+L1+'_L2-'+L2 + '/'
    new_dir_Gg = str(int(Gg)) + '/'
    nld_dir = new_dir_rho + new_dir_L1_L2
    Gg_dir = nld_dir + new_dir_Gg
    os.makedirs(dataset_folder + master_folder + nld_dir, exist_ok = True)
    os.system('cp ' + data_folder + 'TALYS/Sb.tab '+ dataset_folder + master_folder + nld_dir + 'Sb.tab')
    make_TALYS_tab_file(dataset_folder + master_folder + nld_dir + 'Sb.tab', dataset_folder + '127Sb-database/' + nld_dir + 'talys_nld_cnt.txt', A, Z)
    
    #create Gg folder, put E1, M1 and talys input file there - sånn beholder du samme mappestruktur man kan mate i talys
    os.makedirs(dataset_folder + master_folder + Gg_dir)
    E_tal, gsf_tal = make_E1_M1_files(dataset_folder + '127Sb-database/' + Gg_dir, Sn, A, Z, a0, a1, M1 = [Sn127_M1_par[0,0],Sn127_M1_par[1,0],Sn127_M1_par[2,0]], target_folder = dataset_folder + master_folder + Gg_dir, high_energy_interp = Te128mat)
    
    #copy input file to working folder
    os.system('cp ' + dataset_folder + '127Sb ' + dataset_folder + master_folder + Gg_dir + 'input.txt')
    
    #os.system('tar -czf ' + dataset_folder + '/' + master_folder + '/' + Gg_dir + '/' + 'Gg_input.tar.gz .')
    os.chdir(dataset_folder + master_folder + Gg_dir)
    os.system('tar -czf ' + dataset_folder + 'Gg_input.tar.gz ./')
    os.system('rm ./*')
    os.system('mv ' +  dataset_folder + 'Gg_input.tar.gz ./')
    os.chdir(root_folder)
    
    #Don't tar first, count how many gsf pairs I have
    if count%100 == 0:
        print('.')

#zip rho directories
for rho_dir in os.listdir(dataset_folder + master_folder):
    for L1L2 in os.listdir(dataset_folder + master_folder + rho_dir):
        zipping_dir = dataset_folder + master_folder + rho_dir + '/' + L1L2
        os.chdir(zipping_dir)
        os.system('tar -czf ' + dataset_folder + 'L1_L2.tar.gz ./')
        os.system('rm -r ./*')
        os.system('mv ' + dataset_folder + 'L1_L2.tar.gz ./')
        
os.chdir(root_folder)