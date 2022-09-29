# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 13:32:56 2022

@author: bener807
"""

"""
Calculates the pulse areas of the experimental data acquired using the Na-22
source for all S1 and S2 detectors. 

Notes
-----
Integration limits of 10-30 ns are used.
Output is saved in data/pulse_areas.
"""

import sys
sys.path.insert(0, 'C:/python/useful_definitions/')
import useful_defs as udfs
sys.path.insert(0, 'C:/python/TOFu/functions/')
import tofu_functions as dfs
import numpy as np
import matplotlib.pyplot as plt


def import_data(detector, rec_len, path):
    """Return the raw pulse waveform data."""
    board, channel = dfs.get_board_name(detector)
    path = f'{path}/M11D-B{board}_DT{channel}'
    p = udfs.import_pulses(full_path=path, record_length=rec_len)
    return p


def calculate_integral(pulses, bias, rec_len, detector):
    """Return the integral under each pulse waveform."""
    # Perform baseline reduction
    dat = dfs.baseline_reduction(pulses)

    # Remove bad pulses
    dat, _ = dfs.cleanup(dat, dx=1, detector_name=detector, bias_level=bias)

    # Perform sinc interpolation
    dat_sinc = dfs.sinc_interpolation(dat, np.arange(0, rec_len),
                                      np.arange(0, rec_len, 0.1))

    # Perform integration between 10-30 ns
    dat_sinc = dat_sinc[:, 100:300]

    # Integrate
    integrals = -np.trapz(dat_sinc, dx=0.1, axis=1)

    return integrals


def bin_integrals(integrals, bin_edges):
    """Return the binned integrals."""
    counts, bins = np.histogram(integrals, bins=bin_edges)
    bin_centers = bin_edges[1:] - np.diff(bin_edges) / 2

    return counts, bin_centers


def write_to_file(events, bin_centers, detector):
    """Save histogram data to file."""
    to_write = np.array([bin_centers, events]).T
    np.savetxt(f'{detector}.txt', to_write,
               header='Bin centers (areas) Counts')


def plot_histogram(bin_centers, counts, fig_name):
    """Plot the histogram."""
    plt.figure(fig_name)
    plt.plot(bin_centers, counts)
    plt.xlabel('Area')
    plt.ylabel('Counts')


def main(detector, bias_level, rec_len, bin_edges, path):
    """Calculate area under pulses for each detector."""
    # Import pulses
    pulses = import_data(detector, rec_len, path)

    # Calculate area under pulse
    integrals = calculate_integral(pulses, bias_level, rec_len, detector)

    # Bin and write data to file
    counts, bin_centers = bin_integrals(integrals, bin_edges)
    write_to_file(counts, bin_centers, detector)

    return counts, bin_centers


if __name__ == '__main__':
    # Set directory
    directory = '20-11-2020'

    # Select Na22 or background spectrum
    path = f'data/raw_data/Na22/{directory}'
    # path = f'data/raw_data/background/{directory}'

    # Run for ADQ412
    adq412 = dfs.get_dictionaries('ADQ412')
    for detector in adq412:
        counts, bin_centers = main(detector, 1600, 56,
                                   np.arange(0, 5E4, 40), path)
        plot_histogram(bin_centers, counts, 'ADQ412')

    # Run for ADQ14
    adq14 = dfs.get_dictionaries('ADQ14')
    for detector in adq14:
        if detector[:2] == 'S1':
            bias_level = 27000
            bin_width = 300
        else:
            bias_level = 27000
            bin_width = 800

        counts, bin_centers = main(detector, bias_level, 64,
                                   np.arange(0, 1.2345E6, bin_width), path)
        plot_histogram(bin_centers, counts, 'ADQ14')
