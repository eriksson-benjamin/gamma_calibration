# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 16:36:49 2022

@author: bener807
"""

import sys
sys.path.insert(0, 'C:/python/TOFu/functions/')
import tofu_functions as dfs
import numpy as np


def merge_data():
    adq14 = dfs.get_dictionaries('ADQ14')
    counts = np.zeros([411, 20])
    for i, detector in enumerate(adq14):
        path = f'../data/pulse_areas/{detector}.txt'
        p = np.loadtxt(path, dtype='float')
        counts[:, i] = p[:, 1]
    np.savetxt('adq14.txt', counts)

    adq412 = dfs.get_dictionaries('ADQ412')
    counts = np.zeros([142, 17])
    for i, detector in enumerate(adq412):
        path = f'../data/pulse_areas/{detector}.txt'
        p = np.loadtxt(path, dtype='float')
        counts[:, i] = p[:, 1]
    np.savetxt('adq412.txt', counts)

merge_data()