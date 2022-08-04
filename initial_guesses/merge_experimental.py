# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 16:36:49 2022

@author: bener807
"""

import sys
sys.path.insert(0, 'C:/python/TOFu/functions/')
import tofu_functions as dfs
import numpy as np


def merge_data(directory):
    """Merge pulse areas to single file to paste into compton_fit.xlsx."""
    adq14 = dfs.get_dictionaries('ADQ14')
    counts = np.zeros([411, 20])
    for i, detector in enumerate(adq14):
        path = f'../data/pulse_areas/{directory}/{detector}.txt'
        p = np.loadtxt(path, dtype='float')
        counts[:, i] = p[:, 1]
    np.savetxt('adq14.txt', counts)

    adq412 = dfs.get_dictionaries('ADQ412')
    counts = np.zeros([142, 17])
    for i, detector in enumerate(adq412):
        path = f'../data/pulse_areas/{directory}/{detector}.txt'
        p = np.loadtxt(path, dtype='float')
        counts[:, i] = p[:, 1]
    np.savetxt('adq412.txt', counts)


if __name__ == '__main__':
    merge_data('20-11-2020')
