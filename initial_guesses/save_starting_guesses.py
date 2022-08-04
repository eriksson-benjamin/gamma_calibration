# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:11:35 2022

@author: bener807
"""

import pandas as pd
import os
import numpy as np


def read_excel(column, row):
    """Return excel cell value."""
    data = pd.read_excel('compton_fit.xlsx', 'Fits', index_col=None,
                         usecols=column, header=row - 1, nrows=0)
    return data.columns.values[0]


def fit_parameters():
    """Return fit parameters from excel sheet."""
    alpha = read_excel('E', 10)
    beta = read_excel('E', 11)
    gamma = read_excel('E', 12)
    offset = read_excel('E', 15)
    multiplier = read_excel('E', 16)
    primary = read_excel('E', 19)
    secondary = read_excel('E', 20)
    hi = read_excel('L', 10)
    lo = read_excel('L', 11)

    return [alpha, beta, gamma, offset, multiplier, primary, secondary, lo, hi]


def write_to_file(parameters):
    """Write starting guesses to file."""
    data_frame = pd.read_excel('compton_fit.xlsx', 'Fits')
    detector_name = data_frame['Choose measured spectrum'][24]
    file_name = f'parameters/{detector_name}.txt'

    if os.path.exists(file_name):
        u_inp = input(f'Overwrite {file_name} (y/n): ')
        if u_inp != 'y':
            print('File was not overwritten.')
            return 0

    # Write fit parameters to file
    par_names = ['alpha', 'beta', 'gamma', 'offset', 'multiplier', 'I_511',
                 'I_1275', 'lim_lo', 'lim_hi']
    with open(file_name, 'w') as handle:
        for pn, par in zip(par_names, parameters):
            handle.write(f'{pn.ljust(10)} {par}\n')
    print(f'{file_name} saved.')


if __name__ == '__main__':
    parameters = fit_parameters()
    write_to_file(parameters)
