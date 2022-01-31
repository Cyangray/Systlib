#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Thu Oct 14 15:26:05 2021

@author: francesco, last update 31st january 2022

Script to loop through counting.c and normalization.c from the Oslo method.
It creates a directory structure with different input.cnt with different values
of L1, L2 and rho, and runs counting on them. The results are written in each its
folder. Then for each of these, normalization.c is run with different values of
<Gamma_gamma>
Preparations: 
1) in the same folder as this script, put copies of counting.dat, 
rhosp.rsg, sigsp.rsg, rhopaw.rsg, sigpaw.rsg files in the same folder
2) run counting once as you normally would. Copy the values in input.cnt in the 
variable "newinput_cnt", changing the values for L1, L2 and rho with the corresponding 
code variables. Use the example for 127Sb here as an example
3) similarly for normalization: run it once as you normally would, and change the
values in newinput_nrm with the one you got, keeping the name of the python variables
here.

NB: it needs a slightly modified (and recompiled) counting.c and normalization, 
where all lines where it asks input (sscanf, I think) must be commented out, so 
that it only imports data from input.cnt and input.nrm.
to make the program run smoother, comment out as well all printed messages in
the c file.
'''

from subprocess import call
import os
import numpy as np

#initiate lists to loop through
#blist: this is a list of factor you multiply base_rho with, and here for example
#it will create 50 values between 0.4*base_rho and 1.4*base_rho
blist = np.linspace(0.40,1.4,50)
base_rho = 375720
#Gglist: list of different <Gamma_gamma>s. here 13 values between 67.5 and 142.5 meV
Gglist = np.linspace(67.5,142.5,13)
#maximum value for lower bound to be used in counting.
L1max = 23

#Variable names:
#master folder name:
master_folder = '127Sb-database'
#full adress of modified counting and normalization codes
counting_code_path = "/home/francesco/oslo-method-software/prog/counting"
normalization_code_path = "/home/francesco/oslo-method-software/prog/normalization"






def rho2D(rho, target_spin, spin_cutoff):
    '''
    calculate D from rho at Bn (Larsen et al. 2011, eq. (20))
    Takes as input rho as 1/MeV, and gives output D as eV
    target_spin is the spin of the target nucleus
    spin_cutoff is self-explanatory - calculate with robin (?)
    '''
    factor = 2*spin_cutoff**2/((target_spin + 1)*np.exp(-(target_spin + 1)**2/(2*spin_cutoff**2)) + target_spin*np.exp(-target_spin**2/(2*spin_cutoff**2)))
    D = factor/(rho)
    return D*1e6

try:
    os.mkdir(master_folder)
except:
    pass
os.chdir(master_folder)
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
        for L2n in range(L1n + L2_skip, L1max):
            L2 = str(L2n)
            new_dir_L1_L2 = 'L1-'+L1+'_L2-'+L2
            os.mkdir(new_dir_L1_L2)
            os.chdir(new_dir_L1_L2)
            os.system('cp ../../../counting.dat counting.dat')
            os.system('cp ../../../rhosp.rsg rhosp.rsg')
            os.system('cp ../../../sigsp.rsg sigsp.rsg')
            os.system('cp ../../../rhopaw.rsg rhopaw.rsg')
            os.system('cp ../../../sigpaw.rsg sigpaw.rsg')
            
            newinput_cnt = ' 127.000000 1.847000 8.380000 ' + new_rho_str +' 90000.000000 \n ' + L1 + ' ' + L2 + ' 43 49 \n 12 17 54 59 \n 1 12.353000 -0.452000 \n 2 0.800000 -1.723957 \n 1 \n 0 -2.689502 2.672232 \n 0 -0.099645 11.078947 \n 1.000000 \n 0.500000 1.000242 -1.000000 1.000498 \n 300.000000 '
            with open('input.cnt', 'w') as write_obj:
                write_obj.write(newinput_cnt)
            call([counting_code_path])

            
            for Gg in Gglist:
                Ggstr = str(int(Gg))
                Gg_input_str = '{:.6f}'.format(Gg)
                os.mkdir(Ggstr)
                os.chdir(Ggstr)
                
                os.system('cp ../rhosp.rsg rhosp.rsg')
                os.system('cp ../rhotmopaw.cnt rhotmopaw.cnt')
                os.system('cp ../sigextpaw.cnt sigextpaw.cnt')
                os.system('cp ../spincut.cnt spincut.cnt')
                os.system('cp ../sigpaw.cnt sigpaw.cnt')
                os.system('cp ../input.cnt input.cnt')
                
                newinput_nrm = ' 0 8.380000 8.000000 \n ' + new_D_str + ' ' + Gg_input_str + ' \n 105.000000 150.000000 '
                with open('input.nrm', 'w') as write_obj:
                    write_obj.write(newinput_nrm)
                
                call([normalization_code_path]);
                os.system('rm rhosp.rsg')
                os.system('rm rhotmopaw.cnt')
                os.system('rm sigextpaw.cnt')
                os.system('rm spincut.cnt')
                os.system('rm sigpaw.cnt')
                os.system('rm input.cnt')
                os.chdir('..')
                
                
            os.system('rm counting.dat')
            os.system('rm rhosp.rsg')
            os.system('rm sigsp.rsg')
            os.system('rm rhopaw.rsg')
            os.system('rm sigpaw.rsg')
            os.chdir('..')
    print('bstr: ' + bstr)
    os.chdir('..')

     
