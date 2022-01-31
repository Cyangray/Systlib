#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 21:35:23 2022

@author: francesco

Produce a file in the strength.nrm format to be red by the root fitting algorithm
"""

from systlib import *
from sklearn.linear_model import LinearRegression
import numpy as np

Te128 = load_known_gsf(128, 52)
regr = LinearRegression()
indices = np.array(list(range(0,52))).reshape((-1, 1))
energies = Te128.energies
regr.fit(indices, energies)

a0 = regr.intercept_
a1 = regr.coef_[0]
print('a0 = %s, a1 = %s'%(a0,a1))
strength_arr = []
for i in range(102):
    if i < 52:
        strength_arr.append(Te128.y[i])
    else:
        strength_arr.append(Te128.yerr[i-52])
        
np.savetxt('strength128Te.nrm',strength_arr, fmt = '  %.7e')