#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 09:49:41 2021

@author: francesco
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import sys
sys.path.insert(1, '/home/francesco/Documents/Talys-calculations')
from readlib import *
from systlib import *
import seaborn as sns


'''
Imports NLDs created by run_counting.py in the countingexp/ directory, and plots
them all together in a graph. The alpha (transparency) depends on how well each particular
NLD fits to the known levels (by a chi squared test, either linear or logarithmic)
in a given energy range (chi2_lim).
It plots heatmaps of how well the different combinations fit to the chosen fitting interval.
It can also plot the correlation matrices between the three normalization parameters.
It saves the chitable file, to be used in plot_nld_gsf.py
'''

#parameters. play with these
L1max = 23
L2max = 23
plot_graph = False
plotlist=[
    #[1,3],
    #[11,15],
    #[19,22],
    ]
#region of the nld where to evaluate the chi squared test
chi2_lim2 = [14,18] #Fitting interval. Limits are included. [14,18] is a nice region to fit to
chi2_lim = [6,10]
method = 'log'
annotate = False
plot_heatmap = False
plot_correlation_matrix = False

#Don't modify unless you know what you're doing
nrhos = 50
blist = np.linspace(0.53,1.2804,nrhos)
Gglist = np.linspace(80,130,5)
base_rho = 375720
rho_Sn_err = 90000
fermistart=37 #from where to plot the CT model fit
Sn = 8.380
a0 =  -0.7599
a1 =   0.1562

#Initializing main loop. 
#loop idea: loop through all L1 and L2 (the lower region to fit to in counting.c),
#rho and Gg. Each point to a specific rhopaw.cnt produced in run_counting.py. Imports
#the NLD and the GSF, calculate the chi squared test to the region in chi2_lim,
#save all the parameters in the nld and gsf objects, and the objects in lists.
gsfs = []
nlds = []
if plotlist == []:
    plotlist = [[70,70]]
Sb127_nld_lvl = import_ocl('../countingexp/80/053-199131/L1-5_L2-12_B-053_Gg-80/rholev.cnt',a0,a1)
chi2_lim_e = [Sb127_nld_lvl[int(chi2_lim[0]),0], Sb127_nld_lvl[int(chi2_lim[1]),0]]
chi2_lim_e2 = [Sb127_nld_lvl[int(chi2_lim2[0]),0], Sb127_nld_lvl[int(chi2_lim2[1]),0]]
if plot_correlation_matrix:
    #prepare for graph
    fig3, ax3 = plt.subplots(nrows = 3, ncols = 1)
    fig3.set_size_inches(10, 20)
    ax3[0].set_title('Correlation Matrix')
    corr_combinations = [[0,1],[0,2],[1,2]]
    corr_xlims = [[1,3],[1,3],[2e8,4e8]]
    corr_ylims = [[2e8,4e8],[1.4,1.8],[1.4,1.8]]
    corr_combinations_labels = [['A','B'],['A',r'$\alpha$'],['B',r'$\alpha$']]
    #fig3.set_size_inches(20, 10)
#beginning the big nested loop
for Gg in Gglist:
    Ggstr = str(int(Gg))
    for b in blist:
        bstr_int = "{:.2f}".format(b)
        bstr = bstr_int.translate(bstr_int.maketrans('','', '.'))
        new_rho = b * base_rho
        new_rho_str = '{:.6f}'.format(new_rho)
        new_dir_rho = bstr + '-' + str(int(new_rho))
        print('%d, %s'%(Gg, bstr_int)) #print something to show the progression
        
        for L1n in range(1,L1max):
            for L2n in range(L1n+1,L2max):
                if not (L1n == 1 and L2n == 2):
                    #import nlds from the 127Sb experiment handeled with counting, and the fitted CT models.
                    L1 = str(L1n)
                    L2 = str(L2n)
                    nld_dir = 'L1-'+L1+'_L2-'+L2+'_B-' + bstr + '_Gg-' + Ggstr
                    
                    curr_nld = nld('../countingexp/' + Ggstr +'/' + new_dir_rho + '/' + nld_dir + '/rhopaw.cnt',a0 = a0, a1 = a1, is_ocl = True)
                    curr_gsf = gsf('../countingexp/' + Ggstr +'/' + new_dir_rho + '/' + nld_dir + '/strength.nrm', a0 = a0, a1 = a1, is_sigma = False, is_ocl = True)
                    
                    Sb127_nld_ocl = import_ocl('../countingexp/' + Ggstr +'/' + new_dir_rho + '/' + nld_dir + '/rhopaw.cnt',a0,a1)
                    Sb127_nld_fermi = import_ocl('../countingexp/' + Ggstr +'/' + new_dir_rho + '/' + nld_dir +  '/fermigas.cnt',a0,a1,fermi=True)
                    Anorm, alpha = import_Anorm_alpha('../countingexp/' + Ggstr +'/' + new_dir_rho + '/' + nld_dir +  '/alpha.txt')
                    Bnorm = import_Bnorm('../countingexp/' + Ggstr +'/' + new_dir_rho + '/' + nld_dir +  '/input.nrm')
    
                    #calculate the reduced chi2
                    lvl_values = np.concatenate([Sb127_nld_lvl[chi2_lim[0]:(chi2_lim[1]+1),1], Sb127_nld_lvl[chi2_lim2[0]:(chi2_lim2[1]+1),1]])
                    ocl_values = np.concatenate([curr_nld.y[chi2_lim[0]:(chi2_lim[1]+1)], curr_nld.y[chi2_lim2[0]:(chi2_lim2[1]+1)]])
                    ocl_errs = np.concatenate([curr_nld.yerr[chi2_lim[0]:(chi2_lim[1]+1)], curr_nld.yerr[chi2_lim2[0]:(chi2_lim2[1]+1)]])
                    
                    #ocl_values = np.concatenate([Sb127_nld_ocl[chi2_lim[0]:(chi2_lim[1]+1),1], Sb127_nld_ocl[chi2_lim2[0]:(chi2_lim2[1]+1),1]])
                    #ocl_errs = np.concatenate([Sb127_nld_ocl[chi2_lim[0]:(chi2_lim[1]+1),2], Sb127_nld_ocl[chi2_lim2[0]:(chi2_lim2[1]+1),2]])
                    chi2 = chisquared(lvl_values, ocl_values, ocl_errs, DoF=0, method = method)
                    
                    #store values in objects
                    curr_nld.L1 = curr_gsf.L1 = L1
                    curr_nld.L2 = curr_gsf.L2 = L2
                    curr_nld.Anorm = curr_gsf.Anorm = Anorm
                    curr_nld.Bnorm = curr_gsf.Bnorm = Bnorm
                    curr_nld.alpha = curr_gsf.alpha = alpha
                    curr_nld.Gg = curr_gsf.Gg = Gg
                    curr_nld.rho = curr_gsf.rho = new_rho
                    curr_nld.b = curr_gsf.b = b
                    curr_nld.chi2 = curr_gsf.chi2 = chi2
                    
                    #add nld and gsf to list
                    nlds.append(curr_nld)
                    gsfs.append(curr_gsf)
                    #plot simple nlds
                    for L1L2 in plotlist:
                        if L1n == L1L2[0] and L2n == L1L2[1]:
                            plotgraph(L1n,L2n,chi2,Sb127_nld_ocl, Sb127_nld_fermi, Sb127_nld_lvl,new_rho,chi2_lim_e,fermistart = fermistart)
                    
                    
        #plot heatmap
        if plot_heatmap:
            hmap = np.empty((L1max,L2max))
            hmap[:] = np.nan
            for row in chitable:
                if method=='log':
                    #hmap[int(row[0]),int(row[1])] = row[2]
                    hmap[int(row[0]),int(row[1])] = np.log(row[2])
                elif method=='linear':
                    hmap[int(row[0]),int(row[1])] = row[2]
                    #hmap[int(row[0]),int(row[1])] = np.log(row[2])
            fig,ax = plt.subplots()
            ax = sns.heatmap(hmap, linewidth=0.5, cmap = 'magma')
            ax.set_xlabel('L2')
            ax.set_ylabel('L1')
            ax.set_title(r'$\chi^2$, %.2f-%.2f MeV'%(chi2_lim_e[0],chi2_lim_e[1]) + method + r' fit, L1=%d, L2=%d, $\rho$=%d' %((chi2_lim[0]),(chi2_lim[1]),new_rho))
            fig.show()
            
            
#chitable = np.vstack((L1s,L2s,chis,Anorms,Bnorms,alphas)).T 

#plot correlation matrix between Anorm, Bnorm and alpha
if plot_correlation_matrix:
    norm_factor= 0
    anti_chis = [1/nld.chi2 for nld in nlds]
    norm_factor = max(anti_chis)
        
    if nrhos == 5:
        colors = ['b','r','g','y','c']
    else:
        colors = ['k' for x in range(nlds)]
    
    Anorms = [nld.Anorm for nld in nlds]
    Bnorms = [nld.Bnorm for nld in nlds]
    alphas = [nld.alpha for nld in nlds]
    params = np.vstack((Anorms, Bnorms, alphas)).T
    for i, corr_combination in enumerate(corr_combinations):
        for j, nld in enumerate(nlds):
            alpha_list = [(1/nld.chi2)/norm_factor for nld in nlds]
            rgba_colors = np.zeros((len(nlds),4))
            r,g,b,a = to_rgba(colors[j])
            rgba_colors[:,0] = r
            rgba_colors[:,1] = g
            rgba_colors[:,2] = b
            # the fourth column needs to be your alphas
            rgba_colors[:, 3] = alpha_list
            ax3[i].scatter(params[j, corr_combination[0]], params[j, corr_combination[1]] ,color=rgba_colors)
            
            
            #ax3[i].scatter(chitable[:,corr_combination[0]],chitable[:,corr_combination[1]],color=rgba_colors)#, label=r'$\rho$=%d'%rhos[j])
            #ax3[i].scatter(chitable[0,corr_combination[0]],chitable[0,corr_combination[1]],color=colors[j], alpha=1, label=r'$\rho$=%d'%rhos[j])
        
        ax3[i].set_xlabel(corr_combinations_labels[i][0])
        ax3[i].set_ylabel(corr_combinations_labels[i][1])
        ax3[i].set_xlim(corr_xlims[i])
        ax3[i].set_ylim(corr_ylims[i])
        ax3[i].grid()
        #ax3[i].legend()
    #fig3.tight_layout()
    fig3.show()
    

#np.save('chitables.npy', chitables)
np.save('nlds.npy', nlds)
np.save('gsfs.npy', gsfs)












