# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 14:07:09 2022

@author: bener807
"""

import compton as cpt
import matplotlib.pyplot as plt
import numpy as np


def plot_gamma_spectra(detector, directory):
    """
    Plot the experimental/simulated Na-22 spectra for technical paper.

    Notes
    -----
    Simulated gamma spectrum is not broadened here.
    """
    # Read experimental/simulated data
    sim_511, sim_1275 = cpt.read_mcnp(detector)
    exp = cpt.read_experimental(detector, directory)

    # Plot
    # ----
    fig, (ax1, ax2) = plt.subplots(2, 1)
    fig.set_figheight(8)
    fig.set_figwidth(4)
    plt.subplots_adjust(hspace=0.3)

    ax1.plot(exp[:, 0] / 1E6, exp[:, 1], 'k.', markersize=1.5)
    ax1.errorbar(exp[:, 0] / 1E6, exp[:, 1], yerr=np.sqrt(exp[:, 1]),
                 color='k', linestyle='None')
    ax2.plot(sim_511[:, 0], sim_511[:, 1] + sim_1275[:, 1], 'r--')

    # Scale
    ax1.set_yscale('log')
    ax2.set_yscale('log')

    # Labels
    ax1.set_xlabel('Integrated charge (a.u.)')
    ax1.set_ylabel('Counts')
    ax2.set_xlabel('$E_{ee}$ $(MeV_{ee})$')
    ax2.set_ylabel('Counts (a.u.)')

    # Limits
    ax1.set_xlim(-0.05, 0.42)
    ax1.set_ylim(1, 2E4)
    ax2.set_xlim(0.1, 1.4)
    ax2.set_ylim(1E-6, 1E-2)

    # Arrows
    ax1.annotate("", xy=(0.078, 5500), xytext=(0.14, 5500),
                 arrowprops=dict(arrowstyle="->"))
    ax1.annotate("", xy=(0.238, 400), xytext=(0.3, 400),
                 arrowprops=dict(arrowstyle="->"))
    ax2.annotate("", xy=(0.369, 2.8E-3), xytext=(0.52, 2.8E-3),
                 arrowprops=dict(arrowstyle="->"))
    ax2.annotate("", xy=(1.117, 6.3E-4), xytext=(1.268, 6.3E-4),
                 arrowprops=dict(arrowstyle="->"))

    ax1.text(0.9, 0.92, '(a)', transform=ax1.transAxes)
    ax2.text(0.9, 0.92, '(b)', transform=ax2.transAxes)


def read_parameters(detector, directory):
    """Return starting guesses for given detector."""
    p = np.loadtxt(f'output/fit_parameters/{directory}/{detector}.txt',
                   dtype='float', usecols=1)
    parameters = p[:-2]
    limits = p[-2:]

    return parameters, limits


def plot_fitted_spectra(detectors, directory):
    """Plot the experimental Na-22 with the broadened fit applied to it."""
    fig, axes = plt.subplots(2, 2, sharex=False, sharey=False)
    fig.set_figheight(6)
    fig.set_figwidth(6)
    letters = ['(a)', '(b)', '(c)', '(d)']
    x_lims = [(0.14, 1.32), (0.145, 1.4), (0.2, 1.49), (0.25, 1.4)]
    y_lims = [(60, 4.1E4), (50, 1.5E4), (5, 4E3), (60, 2E4)]
    for i, ax in enumerate(axes.flatten()):
        # Read experimental/simulated data and fit parameters
        sim_511, sim_1275 = cpt.read_mcnp(detectors[i])
        exp = cpt.read_experimental(detectors[i], directory)
        parameters, limits = read_parameters(detectors[i], directory)

        # Apply fit parameters
        exp, sim = cpt.warp_axes(parameters, sim_511, sim_1275, exp, limits)

        # Plot
        # ----
        ax.plot(exp[0], exp[1], 'k.', label='Na-22', markersize=1.5)
        ax.errorbar(exp[0], exp[1], yerr=np.sqrt(exp[1]), linestyle='None',
                    color='k')
        ax.plot(sim[0], sim[1], 'r--', label='MCNP')
        ax.set_yscale('log')

        ax.text(0.87, 0.87, letters[i], transform=ax.transAxes)
        ax.set_xlim(x_lims[i])
        ax.set_ylim(y_lims[i])
        ax.set_title(detectors[i].replace('_', '-'), loc='left')
    plt.subplots_adjust(hspace=0.18, wspace=0.18)
    fig.text(0.5, 0.1, '$E_{ee}$ ($MeV_{ee}$)', transform=fig.transFigure)
    fig.text(0.08, 0.52, 'counts', rotation='vertical')


if __name__ == '__main__':
    plot_gamma_spectra('S2_01', '20-11-2020')

    detectors = ['S1_05', 'S2_01', 'S2_02', 'S2_03']
    plot_fitted_spectra(detectors, '20-11-2020')
