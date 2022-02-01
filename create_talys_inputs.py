#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 15:26:05 2021

@author: francesco, updated 1st February 2022

Script to loop through the folders made with run_counting.py, and writes for
each gsf and nld a e1strength, m1strength and Sb.tab file to be used as input by talys,
possibly in Saga. In order to extrapolate for the GDR, the 127Sb gsf is attached
to the 128Te GDR
"""

import os
import numpy as np
from systlib import rho2D, make_E1_M1_files, make_TALYS_tab_file, load_known_gsf

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
master_folder = '127Sb-saga-original'
data_folder = '/home/francesco/Documents/127Sb-experiment/Systlib/data'
dataset_folder = '/home/francesco/Documents/127Sb-experiment/Systlib/Make_dataset'
indexes_to_delete = np.argwhere(Te128.energies<10) #delete values less than 10 MeV
Te128energies = np.delete(Te128.energies, indexes_to_delete, 0) #delete first three rows
Te128y = np.delete(Te128.y, indexes_to_delete, 0) #delete first three rows
Te128mat = np.vstack((Te128energies,Te128y)).T


try:
    os.mkdir(dataset_folder + '/' + master_folder)
except:
    pass
os.chdir(dataset_folder + '/' + master_folder)

for b in blist:
    bstr_int = "{:.2f}".format(b)
    bstr = bstr_int.translate(bstr_int.maketrans('','', '.'))
    new_rho = b * base_rho
    new_D = rho2D(new_rho,8,6.445)
    new_rho_str = '{:.6f}'.format(new_rho)
    new_D_str = '{:.6f}'.format(new_D)
    new_dir_rho = bstr + '-' + str(int(new_rho))
    os.mkdir(new_dir_rho)
    os.chdir(new_dir_rho)
        
    for L1n in range(1,L1max):
        L1 = str(L1n)
        if L1n == 1:
            L2_skip = 2
        else:
            L2_skip = 1
        
        for L2n in range(L1n + L2_skip,L1max):
            L2 = str(L2n)
            new_dir_L1_L2 = 'L1-'+L1+'_L2-'+L2
            os.mkdir(new_dir_L1_L2)
            os.chdir(new_dir_L1_L2)
            
            #make .tab files
            #run(['cp','/home/francesco/Documents/127Sb-experiment/countingexp/Sb.tab','Sb.tab'])
            os.system('cp ' + data_folder + '/TALYS/Sb.tab Sb.tab') #>/dev/null
            make_TALYS_tab_file('Sb.tab', dataset_folder + '/127Sb-database/' + new_dir_rho + '/' + 'L1-'+L1+'_L2-'+L2+ '/talys_nld_cnt.txt', A, Z)
            
            for Gg in Gglist:
                Ggstr = str(int(Gg))
                os.mkdir(Ggstr)
                os.chdir(Ggstr)
                
                #make E1, M1 files
                E1_M1_path = dataset_folder + '/127Sb-database/' + new_dir_rho + '/' + 'L1-'+L1+'_L2-'+L2+ '/' + Ggstr
                E1_M1_taget_path = dataset_folder + '/' + master_folder + '/' + new_dir_rho + '/' + 'L1-'+L1+'_L2-'+L2+ '/' + Ggstr
                E_tal, gsf_tal = make_E1_M1_files(E1_M1_path, Sn, A, Z, a0, a1, M1_frac = 0.1, target_folder = E1_M1_taget_path, high_energy_interp = Te128mat)
                
                #copy input file to working folder
                os.system('cp ' + dataset_folder + '/127Sb input.txt')
                
                #tar the files
                os.system('tar -czf Gg_input.tar.gz .')
                
                #delete the untarred files
                os.system('mv Gg_input.tar.gz ../')
                os.system('rm *')
                os.system('mv ../Gg_input.tar.gz ./')
                os.chdir('..')
            
            #make a tar file in every L1_L2 folder
            os.system('tar -czf L1_L2.tar.gz .')
            
            #delete untarred files
            os.system('mv L1_L2.tar.gz ../')
            os.system('rm -r *')
            os.system('mv ../L1_L2.tar.gz ./')
            os.chdir('..')
    
    os.chdir('..')
    print('bstr: ' + bstr)
