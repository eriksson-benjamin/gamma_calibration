# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 13:47:56 2022

@author: bener807
"""

import sys
sys.path.insert(0, 'C:/python/TOFu/functions/')
import tofu_functions as dfs
import numpy as np


def write_fit_parameters(detectors, directory, file_name):
    """Write the fit parameters (offset and multiplier) to text file."""
    with open(file_name, 'w') as handle:
        handle.write('# Offset 	Multiplier\n')
    for detector in detectors:
        path = f'output/fit_parameters/{directory}/{detector}.txt'
        p = np.loadtxt(path, dtype='str')

        # Grab offset/multiplier, divide by normalization constant
        offset = float(p[3][1]) / float(p[-1][1])
        multiplier = float(p[4][1]) / float(p[-1][1])

        # Write to file
        with open(file_name, 'a') as handle:
            handle.write(f'{offset} {multiplier}\n')
    return p


if __name__ == '__main__':
    # Write S1 parameter file
    detectors = dfs.get_dictionaries('S1')
    p = write_fit_parameters(detectors, '20-11-2020',
                             'energy_calibration_S1.txt')

    # Write S2 parameter file
    detectors = dfs.get_dictionaries('S2')
    p = write_fit_parameters(detectors, '20-11-2020',
                             'energy_calibration_S2.txt')
