# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 10:37:13 2022

@author: bener807
"""

import pandas as pd
import numpy as np
import sys
sys.path.insert(0, 'C:/python/TOFu/functions/')
import tofu_functions as dfs

df = pd.read_excel('S2.xlsx', 'starting_guesses_mcnp', usecols='C:AH')

keys = list(df.keys())
detectors = dfs.get_dictionaries('S2')
for key, detector in zip(keys, detectors):
    vals = list(df[key][77:86])
    print(vals)
    np.savetxt(f'{detector}.txt', vals)