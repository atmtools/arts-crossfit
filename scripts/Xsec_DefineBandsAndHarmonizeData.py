#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:53:05 2020

@author: Manfred Brath

This script splits/defines the band in which the cross section data for each species is split.
Furthermore, for each species and each band the spectral resolution is set to a common resolution.
This is the highest spectral resolution within the specific band.

If desired, the script produces for each species an overview plot of the raw data and the definde bands.

Usage:
Set the corresponding paths within the next cell. Then run the script and follow the instructions within the
console/terminal.

"""

import glob
import json
import os
from copy import deepcopy
from gzip import GzipFile

import numpy as np
from scipy.interpolate import interp1d

import Xsec_aux_functions as xaf

# %% constants, paths

# band configuration filename
config_name = 'band_config_list.json'

# folder of harmonized data
harmonized_folder = '../data/harmonized_data/'

# folder of overview of rawdata and defined bands
plot_folder = '../plots/xsec_rawdata/'

script_path = os.getcwd()

# find the HitranXsec data as json file
filelist = glob.glob(script_path + '/../data/HitranXsecJson/*json.gz')
filelist.sort()

counter = -1
for files in filelist:
    counter = counter + 1
    dummy = os.path.basename(files) + ' ==>> ' + str(counter)
    print(dummy)

print('')
print('----------------------------------------------------------------------')
print('now go through the data')
print('')
# %% ask what do you want?

print('read existing config(0) or make a new one(1)?')
new_config_flag = input('==>>')
new_config_flag = bool(int(new_config_flag))

if new_config_flag == 1:
    print('display raw data for band definition? yes(1) no(0) ')
    overview_flag = input('==>>')
    overview_flag = bool(int(overview_flag))

# %% define/load configdata

if new_config_flag == 1:

    # check if file already exist
    old_flag = os.path.isfile(config_name)

    if old_flag:
        with open(config_name, "r") as read_file:
            config_list_old = json.load(read_file)

    config_list = [[]] * len(filelist)

    cnt = -1
    for filename in filelist:

        cnt = cnt + 1

        data_name = os.path.basename(filename)
        species = data_name.split('.')[0]

        with GzipFile(filename) as f:
            xsec_data = json.loads(f.read().decode("utf-8"))

        number_of_sets = len(xsec_data)

        # print filename
        dummy = os.path.basename(filename) + ' #' + str(cnt)
        print(dummy)

        config_data = [species]

        # get spectral information of data
        wvn_intervalls = np.zeros((number_of_sets, 2))
        for j in range(number_of_sets):
            number_of_freqs = len(xsec_data[j][0]['xsec'])
            number_of_obs = len(xsec_data[j])
            wvn_min = xsec_data[j][0]['wmin']
            wvn_max = xsec_data[j][0]['wmax']
            dw = (wvn_max - wvn_min) / number_of_freqs

            wvn_intervalls[j, 0] = wvn_min
            wvn_intervalls[j, 1] = wvn_max

            print(f"set: {j:.0f}")
            print(f"Ns = {number_of_obs:.0f}")
            print(f"Nf = {number_of_freqs:.0f}")
            print(f"wvn_min = {wvn_min:.2f}")
            print(f"wvn_max = {wvn_max:.2f}")
            print(f"dw      = {dw:.7f}")
            print('')

        # Check overlapp of bands
        overlapp = 0
        if number_of_sets > 1:

            for u in range(number_of_sets):
                for v in range(number_of_sets):

                    if u != v:
                        ol = xaf.getOverlap(wvn_intervalls[u, :], wvn_intervalls[v, :])

                        overlapp = overlapp + int(ol > 0)

        if overlapp == 0:
            for j in range(number_of_sets):
                config_data.append(wvn_intervalls[j, :].tolist())

            print('As we have no overlapp, we take the values directly from the data.')
            config_list[cnt] = config_data

            if overview_flag == True:

                # Plot rawdata with new bands
                fig1, ax1 = xaf.plot_raw_data(xsec_data, species)

                for j in range(0, len(ax1)):
                    ax1[j] = xaf.make_band_patches(ax1[j], config_list[cnt][1:],
                                                   [1e-24, 1e-15], edgecolor='None')

        else:

            if overview_flag == True:
                # #Show overview plots
                fig1, ax1 = xaf.plot_raw_data(xsec_data, species)

            input1 = 0
            if old_flag:

                print('Existing band definition')
                for band_idx in range(1, len(config_list_old[cnt])):
                    print(f"band: {band_idx:.0f}")
                    print(config_list_old[cnt][band_idx])

                if overview_flag == True:
                    for j in range(0, len(ax1)):
                        ax1[j] = xaf.make_band_patches(ax1[j], config_list_old[cnt][1:],
                                                       [1e-24, 1e-15], edgecolor='None')

                    xaf.plt.show()

                print('to use old values press enter or to define new ones type number of bands')
                N_bands = input('==>>')

                if len(N_bands) == 0:
                    config_list[cnt] = config_list_old[cnt]
            else:
                print('Number of bands, leave empty for skipping')
                N_bands = input('==>>')

            if len(N_bands) > 0:

                try:
                    N_bands = int(N_bands)

                    if N_bands > number_of_sets:
                        raise ('really so many bands??')
                except:
                    print('I think the number of bands was not correct!')
                    print('Try again!')
                    N_bands = input('==>>')
                    N_bands = int(N_bands)

                for j in range(N_bands):

                    print('')
                    print('input the edges of band ' + str(j))
                    print('wvn_min, wvn_max [cm^-1]')
                    band_specs_in = input('==>>')
                    band_specs_in = band_specs_in.split(',')

                    try:
                        if len(band_specs_in) != 2:
                            raise ('interval definition not correct')

                    except:
                        print('Interval definition not correct')
                        print('Try again!')
                        band_specs_in = input('==>>')
                        band_specs_in = band_specs_in.split(',')

                    dump = [float(idx) for idx in band_specs_in]
                    config_data.append(dump)

                config_list[cnt] = config_data

                if overview_flag == True:

                    # Plot rawdata with new bands
                    fig1, ax1 = xaf.plot_raw_data(xsec_data, species)

                    for j in range(0, len(ax1)):
                        ax1[j] = xaf.make_band_patches(ax1[j], config_list[cnt][1:],
                                                       [1e-24, 1e-15], edgecolor='None')

                    xaf.plt.show()

        print('')
        print('-------------------------------------------------------------------')
        if overview_flag == True:

            if not os.path.exists(plot_folder):
                os.makedirs(plot_folder)

            plotname1 = plot_folder + species + '.png'
            fig1.savefig(plotname1, dpi=300, bbox_inches='tight')
            xaf.plt.close(fig=fig1)

    with open(config_name, 'w') as fout:
        json.dump(config_list, fout)


else:

    with open(config_name, "r") as read_file:
        config_list = json.load(read_file)

# %% harmonize data in frequency

print('   ')
print('------------------------------------------------------------------')
print('start harmonizing the data...')
print('   ')

for f_i in range(len(filelist)):

    filename = filelist[f_i]
    data_name = os.path.basename(filename)
    species = data_name.split('.')[0]

    print('   ')
    print('------------------------------------------------------------------')
    print('Species: ' + species)
    print('file: ' + filename)
    print('   ')

    with GzipFile(filename) as f:
        xsec_data = json.loads(f.read().decode("utf-8"))

    config = [config_i for config_i in config_list if config_i[0] == species][0]

    # loop over defined bands
    wvn_max_previous = np.inf
    for i in range(len(config) - 1):

        # Wave number interval of defined band
        wvn_min = config[i + 1][0]
        wvn_max = config[i + 1][1]
        wvn_int_def = np.array([wvn_min, wvn_max])

        print('Defined band #' + str(i))
        print('wvn min = ' + str(wvn_min))
        print('wvn max = ' + str(wvn_max))
        print('  ')

        # get the highest spectral resolution of the data for the defined band
        dw = []
        Nfs = []
        # loop over bands (data)
        for j in range(len(xsec_data)):

            # loop over spectra within band (data)
            for k in range(len(xsec_data[j])):

                wmax_data = xsec_data[j][k]['wmax']
                wmin_data = xsec_data[j][k]['wmin']
                Nf = len(xsec_data[j][k]['xsec'])

                # Wave number interval of data
                wvn_int_data = np.array([wmin_data, wmax_data])

                ol = xaf.getOverlap(wvn_int_def, wvn_int_data)

                if ol == 0:
                    continue
                else:
                    dwtemp = (wmax_data - wmin_data) / Nf

                    dw.append(dwtemp)
                    Nfs.append(Nf)

        dw = np.array(dw)
        Nfs = np.array(Nfs)

        # Highest spectral resolution
        delta_w = np.min(dw)

        wvn = np.linspace(wvn_min, wvn_max, round((wvn_max - wvn_min) / delta_w))

        # check that the edges are not double
        if wvn_min == wvn_max_previous:
            wvn = wvn[1:]

        wvn_max_previous = wvn_max

        # allocate
        xsec_band_i = np.zeros((len(dw), len(wvn)))
        data_band_i = []

        # now we harmonize the data
        # loop over bands (data)
        cnt_data = -1

        for j in range(len(xsec_data)):

            # loop over spectra within band (data)
            for k in range(len(xsec_data[j])):

                wmax_data = xsec_data[j][k]['wmax']
                wmin_data = xsec_data[j][k]['wmin']
                Nf = len(xsec_data[j][k]['xsec'])
                xsec_data_jk = np.array(xsec_data[j][k]['xsec'])  # in [cm^2]

                # Wave number interval of data
                wvn_int_data = np.array([wmin_data, wmax_data])

                ol = xaf.getOverlap(wvn_int_def, wvn_int_data)

                if ol == 0:
                    continue
                else:
                    cnt_data = cnt_data + 1

                    if cnt_data + 1 > len(dw):
                        print('verdammt')
                        raise (' damn it!')

                    # wavenumber of data
                    wvn_data = np.linspace(wmin_data, wmax_data, Nf)

                    # check for negative values in data
                    logic_neg_data = xsec_data_jk < 0

                    if np.sum(logic_neg_data) > 0:
                        print(str(np.sum(logic_neg_data)) + ' Negative x-sections encountered in data \n')

                    # set up interpolation object
                    fun = interp1d(wvn_data, xsec_data_jk, kind='linear', bounds_error=False)

                    # the actual interpolation
                    xsec_temp = fun(wvn)

                    # just for checking
                    xsec_int_chk = np.trapz(xsec_temp, wvn)

                    # integrate cross sections in units of cm²/cm
                    xsec_data_interval = deepcopy(xsec_data_jk)
                    logic_interval = np.logical_and(wvn_min <= wvn_data, wvn_data < wvn_max)
                    xsec_data_interval[~logic_interval] = 0.
                    xsec_int_temp = np.trapz(xsec_data_interval, wvn_data)

                    # difference between both integrations
                    dxsec_int_temp = xsec_int_chk / xsec_int_temp - 1

                    # check for negative values after interpolations
                    logic_neg = xsec_temp < 0
                    if np.sum(logic_neg) > 0:
                        print(str(np.sum(logic_neg)) + ' Negative x-sections encountered after interpolation')

                    xsec_temp[logic_neg] = 0.
                    xsec_band_i[cnt_data, :] = xsec_temp

                    temp = {}
                    temp['xsec'] = list(xsec_temp / 10000)  # store them in [m^2]
                    temp['wmin'] = wvn[0]
                    temp['wmax'] = wvn[-1]
                    temp['fmin'] = wvn[0] * xaf.c0 * 100
                    temp['fmax'] = wvn[-1] * xaf.c0 * 100
                    temp['temperature'] = np.float64(xsec_data[j][k]['temperature'])
                    temp['xscfile'] = xsec_data[j][k]['xscfile']
                    temp['species'] = xsec_data[j][k]['species']
                    temp['IntXsec_cm2_per_cm'] = xsec_int_temp
                    temp['DeltaIntXsec_relative'] = dxsec_int_temp

                    if xsec_data[j][k]['pressure'] == 0.0:
                        temp['pressure'] = np.float64(101325.)
                    else:
                        temp['pressure'] = np.float64(xsec_data[j][k]['pressure'])

                    data_band_i.append(temp)

        # export harmonized data
        harmonized_data_name = species + '.band' + str(i) + '.xsec'

        with GzipFile(harmonized_folder + harmonized_data_name + '.json.gz', "w") as f:
            f.write(json.dumps(data_band_i).encode("utf-8"))
