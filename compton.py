# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 13:36:24 2022

@author: bener807
"""


import sys
import scipy.optimize as optimize
from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, 'C:/python/useful_definitions/')
import useful_defs as udfs
import scipy
udfs.set_nes_plot_style()
sys.path.insert(0, 'C:/python/TOFu/functions/')
import tofu_functions as dfs
import os


'''
Takes a simulated and a measured Na-22 spectrum and applies alpha, beta, gamma
broadening to the simulation as well as normalisation of the intensities
of 511 keV and 1275 keV gammas. Offset and multiplier is applied to the x-axis
of the measured spectrum. Least squares fit is performed by varying all seven
parameters (alpha, beta, gamma, I_511, I_1275, offset, multiplier).
'''


def read_mcnp(detector):
    """Return simulated 511 keV and 1275 keV gamma spectra."""
    path = 'MCNP/output'
    if detector[:2] == 'S1':
        sim_511 = np.loadtxt(f'{path}/S1/binned_{detector[3:]}_511.txt')
        sim_1275 = np.loadtxt(f'{path}/S1/binned_{detector[3:]}_1275.txt')

    elif detector[:2] == 'S2':
        sim_511 = np.loadtxt(f'{path}/S2/binned_511.txt')
        sim_1275 = np.loadtxt(f'{path}/S2/binned_1275.txt')
    else:
        raise NameError(f'Unknown detector: {detector}')

    return sim_511, sim_1275


def read_experimental(detector, directory):
    """Return the experimental Na-22 data."""
    path_exp = f'data/pulse_areas/Na22/{directory}'
    path_bgr = f'data/pulse_areas/background/{directory}'

    exp = np.loadtxt(f'{path_exp}/{detector}.txt')
    bgr = np.loadtxt(f'{path_bgr}/{detector}.txt')
    return exp, bgr


def starting_guesses(detector, directory):
    """Return starting guesses for given detector."""
    p = np.loadtxt(f'initial_guesses/parameters/{directory}/{detector}.txt',
                   dtype='float', usecols=1)
    parameters = p[:-2]
    limits = p[-2:]

    return parameters, limits


def resolution_broadening(sim_511, sim_1275, alpha, beta, gamma,
                          I_511, I_1275):
    """Calculate the resolution broadening matrix."""
    # Shifted/unshifted matrix energy axis ()
    if np.all(sim_511[:, 0] == sim_1275[:, 0]):
        # Energy axes
        erg_axis_1 = sim_511[:, 0]
        erg_axis_2 = np.append([-100], erg_axis_1[:-1])

        # Counts
        cts_511 = sim_511[:, 1]
        cts_1275 = sim_1275[:, 1]
    else:
        raise ValueError('Simulated energy axes are non-identical.')

    # Empty matrix
    matrix = np.zeros(2 * [len(erg_axis_1)])
    sim_total = np.zeros(len(erg_axis_1))

    for i, (erg_1, erg_2) in enumerate(zip(erg_axis_1, erg_axis_2)):
        # Calculate resolution for given alpha/beta/gamma
        resolution = np.sqrt(alpha**2 + beta**2 / erg_1 + gamma**2 / erg_1**2)

        # Calculate resolution broadening matrix
        sigma = resolution * erg_axis_1 / (2 * np.sqrt(2 * np.log(2)))
        cdf_1 = norm.cdf(x=erg_1, loc=erg_axis_1, scale=sigma)
        cdf_2 = norm.cdf(x=erg_2, loc=erg_axis_1, scale=sigma)
        matrix[:, i] = (I_511 * cts_511 + I_1275 * cts_1275) * (cdf_1 - cdf_2)
        sim_total[i] = np.sum(matrix[:, i]) * 1E+6

    return sim_total, erg_axis_1


def warp_axes(parameters, sim_511, sim_1275, exp, bgr, limits):
    """
    Broaden simulation and warp axis of experimental data.

    Notes
    -----
    - Broaden simulated spectrum and apply offset/multiplier to experimental
      spectrum.
    - Axis is normalized to ensure offset/multiplier are reasonably
      large for minimization algorithm.
    - To apply offset/multiplier to data after gamma calibration we need to
      divide by the normalization (i.e. bin_norm multiplied by bin width)
    """
    # Get free parameters
    alpha = parameters[0]
    beta = parameters[1]
    gamma = parameters[2]
    offset = parameters[3]
    multiplier = parameters[4]
    I_511 = parameters[5]
    I_1275 = parameters[6]

    # Calculate broadened simulated spectrum
    sim_total, sim_axis = resolution_broadening(sim_511, sim_1275, alpha, beta,
                                                gamma, I_511, I_1275)

    # Apply offset/multiplier to experimental binning, normalize binning
    bin_nums = np.arange(0, len(exp))
    bin_norm = 1000
    exp_axis = (offset + multiplier * bin_nums) / bin_norm
    return (exp_axis, exp[:, 1]), (exp_axis, bgr[:, 1]), (sim_axis, sim_total)


def fit_function(parameters, sim_511, sim_1275, exp, bgr, limits):
    """Minimize chi^2."""
    exp, bgr, sim = warp_axes(parameters, sim_511, sim_1275, exp, bgr, limits)
    exp_axis = exp[0]
    exp_counts = exp[1]
    sim_axis = sim[0]
    sim_counts = sim[1]
    bgr_counts = bgr[1]

    # Cut out fit range for comparison
    lims = (np.searchsorted(exp_axis, limits[0]),
            np.searchsorted(exp_axis, limits[1]))
    exp_ax = exp_axis[lims[0]:lims[1]]
    exp_cts = exp_counts[lims[0]:lims[1]]
    bgr_cts = bgr_counts[lims[0]:lims[1]]

    # Interpolate between simulated/experimental data
    sim_cts = udfs.interpolate_new_axis(exp_ax, sim_axis, sim_counts)

    # Add background spectrum
    sim_cts += bgr_cts

    # Calculate chi2
    chi2 = np.sum((exp_cts - sim_cts)**2 / sim_cts) / (len(exp_cts) - 7)

    print(f'chi2 = {chi2}')

    return chi2


def plot_comparison(parameters, sim_511, sim_1275, exp, bgr, limits, detector):
    """
    Plot MCNP spectrum against experimental Na-22 spectrum.

    Notes
    -----
    1. MCNP spectrum is broadened.
    2. Measured background component is added to MCNP spectrum.
    """
    # Warp experimental axis
    exp, bgr, sim = warp_axes(parameters, sim_511, sim_1275, exp, bgr, limits)

    # Interpolate simulated axis onto experimental
    sim_cts = udfs.interpolate_new_axis(exp[0], sim[0], sim[1])

    # Plot
    plt.figure()

    # Experimental data
    plt.plot(exp[0], exp[1], 'k.', label='Na-22', markersize=1)
    plt.errorbar(exp[0], exp[1], yerr=np.sqrt(exp[1]), linestyle='None',
                 color='k')

    # Simulated data
    plt.plot(sim[0], sim[1], 'C0--', label='MCNP')

    # Simulated data + background
    plt.plot(bgr[0], sim_cts + bgr[1], 'r-', label='MCNP+background')

    # Background component
    plt.plot(bgr[0], bgr[1], 'C1-.', label='background')

    # Configure plots
    plt.title(detector, loc='left')
    plt.xlabel('$E_{ee}$ $(MeV_{ee})$')
    plt.ylabel('Intensity (a.u.)')
    plt.yscale('log')
    plt.legend()


def write_parameters(detector, directory, parameters):
    """Write fit parameters to file."""
    par_names = ['alpha', 'beta', 'gamma', 'offset', 'multiplier', 'I_511',
                 'I_1275', 'lim_lo', 'lim_hi', 'norm_factor']
    file_name = f'output/fit_parameters/{directory}/{detector}.txt'

    # Check if file already exists
    if os.path.exists(file_name):
        ans = input(f'Overwrite {file_name}? (y/n): ')
        if ans != 'y':
            print('File was not overwritten.')
            return 0

    # Save parameters to file
    with open(file_name, 'w') as handle:
        for pn, par in zip(par_names, parameters):
            handle.write(f'{pn.ljust(11)} {par}\n')


def main(detector, directory):
    """Perform fit and write parameters to file."""
    # Read simulated/experimental data and initial parameter guesses
    sim_511, sim_1275 = read_mcnp(detector)
    exp, bgr = read_experimental(detector, directory)
    parameters, limits = starting_guesses(detector, directory)

    # Plot initial guess
    while True:
        plot_comparison(parameters, sim_511, sim_1275,
                        exp, bgr, limits, detector)
        plt.savefig('initial_guesses/init.png')
        # plt.close('all')
        ans = input('Continue with fit? (y/n) ')
        if ans == 'y':
            break
        else:
            parameters, limits = starting_guesses(detector, directory)

    # Minimize test statistic with given bounds
    bnds = ((0, None), (0, None), (0, None), (None, None), (0, None),
            (0, None), (0, None))
    popt = scipy.optimize.minimize(fit_function, parameters, bounds=bnds,
                                   args=(sim_511, sim_1275, exp, bgr, limits))

    # Plot resulting fit
    plot_comparison(popt.x, sim_511, sim_1275, exp, bgr, limits, detector)

    # Print parameters
    parameters = np.append(popt.x, limits)

    # Normalization factor is 1/(bin_width*bin_norm)
    bin_norm = 1000.
    norm_factor = (np.diff(exp[:, 0])[0] * bin_norm)
    parameters = np.append(parameters, norm_factor)

    par_names = ['alpha', 'beta', 'gamma', 'offset', 'multiplier', 'I_511',
                 'I_1275', 'lim_lo', 'lim_hi', 'norm_factor']

    for pn, par in zip(par_names, parameters):
        print(f'{pn.ljust(11)} {par}\n')

    # Save parameters to file
    write_parameters(detector, directory, parameters)


if __name__ == '__main__':
    # Main analysis
    detector = sys.argv[1]
    directory = sys.argv[2]
    main(detector, directory)
