#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:53:05 2020

@author: Manfred Brath

This script calculates the absorption cross section fit coefficients from the harmonized
Hitran absorption cross sections and stores the output for each species as ARTS-XsecRecord.
If desired, the script produces 9 figures per species and band as overview of the fitting
process.

Usage:
Set the options at the end of this file and run the script.

"""

import glob
import json
import os
from gzip import GzipFile

import numpy as np
import pyarts
from typhon import constants

import Xsec_aux_functions as xaf


# %%

def process_xsec_coefficients(species, harmonized_folder, coeff_folder, main_plot_folder, store_coeffs=True,
                              plotting=True):
    script_path = os.getcwd()

    # get the files
    filelist = glob.glob(script_path + '/' + harmonized_folder + species + '.*.json.gz')
    filelist.sort()

    print('let us see, what we got!')
    print('  ')
    counter = -1
    for files in filelist:
        counter = counter + 1
        dummy = os.path.basename(files) + ' # ' + str(counter)
        print(dummy)

    print('')
    print('----------------------------------------------------------------------')
    print('now go through the data')
    print('')

    # %% load harmonized

    Xsec_processed_data_array = [[]] * len(filelist)
    min_pressures = np.zeros((len(filelist)))
    max_pressures = np.zeros((len(filelist)))
    min_temperatures = np.zeros((len(filelist)))
    max_temperatures = np.zeros((len(filelist)))

    for f_i in range(len(filelist)):

        filename = filelist[f_i]
        data_name = os.path.basename(filename)
        species = data_name.split('.')[0]
        band_no = int(data_name.split('.')[1][4:])

        # remove minus sign of species name, because ARTS species tags are defined without them
        dump = species.split('-')
        species_arts = ''
        species_arts = species_arts.join(dump)

        print('   ')
        print('------------------------------------------------------------------')
        print('Species: ' + species)
        print('Band #' + str(band_no))
        print('file: ' + filename)
        print('   ')

        with GzipFile(filename) as f:
            data = json.loads(f.read().decode("utf-8"))

        wvn = np.linspace(data[0]['wmin'], data[0]['wmax'], len(data[0]['xsec']))

        # allocate
        fit_goodness = [[]] * len(wvn)
        fit_coeffs = np.zeros((6, len(wvn)))
        N_data = np.zeros(len(wvn))
        MinP = np.zeros(len(wvn)) * np.nan
        MaxP = np.zeros(len(wvn)) * np.nan
        MinT = np.zeros(len(wvn)) * np.nan
        MaxT = np.zeros(len(wvn)) * np.nan

        T = np.zeros(len(data))
        P = np.zeros(len(data))

        for i in range(len(data)):
            T[i] = data[i]['temperature']
            P[i] = data[i]['pressure']

        for idx in range(len(wvn)):

            if idx % 100 == 0:
                print(str(idx))

            # allocate
            Xsec = np.zeros(len(data))

            for i in range(len(data)):
                Xsec[i] = data[i]['xsec'][idx]

            fit_result = xaf.fit_xsec_data(T, P, Xsec, min_deltaLogP=0.2, min_deltaT=20.)

            fit_coeffs[:, idx] = fit_result['coefficients']
            N_data[idx] = fit_result['NumberOfPoints']
            MinP[idx] = fit_result['MinP']
            MaxP[idx] = fit_result['MaxP']
            MinT[idx] = fit_result['MinT']
            MaxT[idx] = fit_result['MaxT']
            fit_goodness[idx] = fit_result['residuum']

        # %% prepare output

        y_var = 'y=log10(p/p_0]), p_0=1Pa'

        # frequency in Hz
        freq = wvn * constants.c * 100

        s_data_temp = pyarts.classes.GriddedField2()
        s_data_temp.gridnames = ['frequency grid [Hz]', 'fit coefficients [m]']
        s_data_temp.grids = [pyarts.classes.Vector(freq),
                             pyarts.classes.ArrayOfString(['p00', 'p10', 'p01', 'p20', 'p11', 'p02'])]
        s_data_temp.data = fit_coeffs.transpose()
        s_data_temp.name = (species + '-band_' + str(band_no))
        s_data_temp.dataname = ('species: ' + species_arts + ' \n'
                                + 'fit of hitran cross section data \n'
                                + 'version XXX \n'
                                + 'fit model: \n'
                                + '    xsec=(' + fit_result['formula'] + ')**2 \n'
                                + '    x=T/T_0, T_0=1K \n'
                                + '    ' + y_var + ' \n'
                                + '    xsec in m**2')

        Xsec_processed_data_array[f_i] = s_data_temp
        min_pressures[f_i] = np.min(MinP)
        max_pressures[f_i] = np.max(MaxP)
        min_temperatures[f_i] = np.min(MinT)
        max_temperatures[f_i] = np.max(MaxT)

        # %% same validation

        if plotting:

            RMSE = np.zeros(np.shape(T)) * np.nan
            bias = np.zeros(np.shape(T)) * np.nan
            StDev = np.zeros(np.shape(T)) * np.nan
            XsecMean = np.zeros(np.shape(T)) * np.nan
            XsecSum = np.zeros(np.shape(T)) * np.nan
            XsecMeanFit = np.zeros(np.shape(T)) * np.nan
            XsecSumFit = np.zeros(np.shape(T)) * np.nan
            XsecInt = np.zeros(np.shape(T)) * np.nan
            DeltaXsecInt = np.zeros(np.shape(T)) * np.nan
            XsecIntFit = np.zeros(np.shape(T)) * np.nan
            XsecIntFitFull = np.zeros(np.shape(T)) * np.nan
            R2 = np.zeros(np.size(fit_coeffs, axis=1))
            XsecTest = np.zeros((np.size(T), np.size(fit_coeffs, axis=1))) * np.nan
            XsecTestFit = np.zeros((np.size(T), np.size(fit_coeffs, axis=1))) * np.nan
            logic_notnan = np.zeros((np.size(T), np.size(fit_coeffs, axis=1)), dtype=bool)

            for i in range(len(RMSE)):
                # fitted spectrum
                XsecTestFit[i, :] = xaf.calculate_xsec(T[i], P[i], fit_coeffs)

                # observed
                XsecTest[i, :] = data[i]['xsec']

                # Integrated Xsec from fit in cm^2/cm
                logic_i = np.isnan(XsecTest[i, :])
                temp = XsecTestFit[i, :] * 1
                temp[logic_i] = 0.
                XsecIntFit[i] = np.trapz(temp, wvn) * 1e4  # Integrate, only where data exists
                XsecIntFitFull[i] = np.trapz(XsecTestFit[i, :], wvn) * 1e4  # Integrate full spectrum

                # Integrated Xsec from rawdata in cm^2/cm
                XsecInt[i] = data[i]['IntXsec_cm2_per_cm']

                # Error of the rawdata due to interpolation
                DeltaXsecInt[i] = data[i]['DeltaIntXsec_relative']

                logic_notnan[i, :] = ~np.isnan(XsecTest[i, :])

                RMSE[i] = np.sqrt(np.nanmean((XsecTestFit[i, :] - XsecTest[i, :]) ** 2))
                bias[i] = np.nanmean(XsecTestFit[i, :] - XsecTest[i, :])
                StDev[i] = np.nanstd(XsecTest[i, :])
                XsecMean[i] = np.nanmean(XsecTest[i, :])
                XsecSum[i] = np.nansum(XsecTest[i, :])
                XsecMeanFit[i] = np.nanmean(XsecTestFit[i, :])
                XsecSumFit[i] = np.nansum(XsecTestFit[i, :])

            # Calculate Goodnes of fit
            for i in range(np.size(fit_coeffs, axis=1)):
                R2[i] = xaf.calc_Rsquare(XsecTest[:, i], XsecTestFit[:, i], np.sum(fit_coeffs[:, i] > 0))

            # %% some preparation for  plotting

            # make random selection of data to show
            nop = 9
            index = np.random.permutation(len(T))
            index = index[0:(np.min([len(data), nop]))]

            print('Start with plotting!')

            # %% plot xsec for some points

            print('Plotting fig0 ')

            fig0, axs0 = xaf.default_figure(3, 3)

            for i in range(len(index)):

                row = int(i / 3)
                col = i - (row * 3)

                ax = axs0[row, col]

                p = P[i]
                t = T[i]

                # fitted spectrum
                # xsec_fit = xaf.calculate_xsec(t, p, fit_coeffs, density_flag=density_flag)
                xsec_fit = xaf.calculate_xsec_fullmodel(t, p, fit_coeffs,
                                                        minT=min_temperatures[f_i],
                                                        maxT=max_temperatures[f_i],
                                                        minP=min_pressures[f_i],
                                                        maxP=max_pressures[f_i])

                # observed
                xsec = data[i]['xsec']

                plot_title = f"{p:.2f}$\,$Pa $-$ {t:.0f}$\,$K"

                if col == 0:
                    ylabel = '$a_{xsec}$ m$^2$'
                else:
                    ylabel = None

                if row == 2:
                    xlabel = 'Wavenumber [cm$^{-1}$]'
                else:
                    xlabel = None

                if row == 0 and col == 0:
                    legend = True
                else:
                    legend = False

                axs0[row, col] = xaf.plot_xsec(wvn, xsec, xsec_fit, ax, xlim=None, xlabel=xlabel, ylabel=ylabel,
                                               plot_title=plot_title, legend=legend)

            # %% plot differences for some points

            print('Plotting fig1 ')

            fig1, axs1 = xaf.default_figure(3, 3)

            for i in range(len(index)):

                row = int(i / 3)
                col = i - (row * 3)

                ax = axs1[row, col]

                p = P[i]
                t = T[i]

                # fitted spectrum
                xsec_fit = xaf.calculate_xsec(t, p, fit_coeffs)

                # observed
                xsec = data[i]['xsec']

                # difference
                dxsec = xsec - xsec_fit

                plot_title = f"{p:.2f}$\,$Pa $-$ {t:.0f}$\,$K"

                if col == 0:
                    ylabel = '$\Delta a_{xsec}$ m$^2$'
                else:
                    ylabel = None

                if row == 2:
                    xlabel = 'Wavenumber [cm$^{-1}$]'
                else:
                    xlabel = None

                axs1[row, col] = xaf.plot_xsec(wvn, dxsec, [], ax, xlim=None, xlabel=xlabel, ylabel=ylabel,
                                               plot_title=plot_title)

            # %% plot overview

            print('Plotting fig2')

            fig2, axs2 = xaf.default_figure(3, 3, sharey='all', sharex='all',
                                            width_in_cm=29.7, height_in_cm=29.7)

            z11 = RMSE / XsecMean * 100
            fig2, axs2[0, 0] = xaf.scatter_plot(T, P, z11, fig2, axs2[0, 0], clim=[0, 40],
                                                plot_title='$\\frac{RMSE}{<a_{xsec,obs}>}$',
                                                cbar_label='[$\%$]')
            axs2[0, 0].invert_yaxis()

            z12 = bias / XsecMean * 100
            fig2, axs2[0, 1] = xaf.scatter_plot(T, P, z12, fig2, axs2[0, 1], clim=[-5, 5],
                                                plot_title='$\\frac{<\\Delta a_{xsec}>}{<a_{xsec,obs}>}$',
                                                cbar_label='[$\\%$]', cmap='difference')

            z21 = StDev
            fig2, axs2[1, 0] = xaf.scatter_plot(T, P, z21, fig2, axs2[1, 0],
                                                plot_title='std($a_{xsec,obs}$)',
                                                cbar_label='[m$^2$]')

            z22 = XsecInt
            z22_max = np.max([XsecInt.max(), XsecIntFit.max()])
            z22_min = np.min([XsecInt.min(), XsecIntFit.min()])
            fig2, axs2[1, 1] = xaf.scatter_plot(T, P, z22, fig2, axs2[1, 1], clim=[z22_min, z22_max],
                                                plot_title='$\int a_{xsec,obs,raw}$',
                                                cbar_label='[cm$^2$ cm$^{-1}$]')

            z31 = XsecIntFitFull
            fig2, axs2[2, 0] = xaf.scatter_plot(T, P, z31, fig2, axs2[2, 0], clim=[z22_min, z22_max],
                                                plot_title='$\int a_{xsec,fit}$',
                                                cbar_label='[cm$^2$ cm$^{-1}$]')

            z32 = (XsecIntFit / XsecInt - 1) * 100
            z32_limit = np.max([abs(z32.min()), abs(z32.max())])
            fig2, axs2[2, 1] = xaf.scatter_plot(T, P, z32, fig2, axs2[2, 1], clim=[-z32_limit, z32_limit],
                                                plot_title='$(\int a_{xsec,fit}/\int a_{xsec,obs,raw})-1$',
                                                cbar_label='[$\\%$]', cmap='difference')

            z33 = DeltaXsecInt * 100
            z33_limit = np.max([abs(z33.min()), abs(z33.max())])
            fig2, axs2[2, 2] = xaf.scatter_plot(T, P, z33, fig2, axs2[2, 2], clim=[-z33_limit, z33_limit],
                                                plot_title='$(\int a_{xsec,obs,inter}/\int a_{xsec,obs,raw})-1$',
                                                cbar_label='[$\\%$]', cmap='difference')

            # %% plot coefficients

            print('Plotting fig3')

            fig3, axs3 = xaf.default_figure(3, 2, sharey='none', sharex='all')

            coeff_names = fit_result['coeff_names']
            formula = fit_result['formula']

            for i in range(np.size(fit_coeffs, axis=0)):

                row = int(i / 2)
                col = i - (row * 2)

                ax = axs3[row, col]

                plot_title = coeff_names[i]

                if col == 0:
                    ylabel = '[m]'
                else:
                    ylabel = None

                if row == 2:
                    xlabel = 'Wavenumber [cm$^{-1}$]'
                else:
                    xlabel = None

                axs3[row, col] = xaf.plot_xsec(wvn, fit_coeffs[i, :], [], ax, xlim=None, xlabel=xlabel, ylabel=ylabel,
                                               plot_title=plot_title)

            y_var = '$y=\\log_{10} \\frac{p}{p_0}$'

            sup_text = species + ': $a_{xsec}=(' + formula + ')^2$; $x=T/T_0$, ' + y_var

            fig3.suptitle(sup_text)

            # %% edge cases

            print('Plotting fig4')

            T_i = [np.min(MinT) * 0.9, np.max(MaxT) / 0.9]
            P_i = [10., np.min(MinP), np.max(MaxP)]

            fig4, axs4 = xaf.default_figure(len(T_i) * len(P_i), 1, width_in_cm=20.9, height_in_cm=29.7)
            fig4.subplots_adjust(hspace=0.5)

            for i in range(len(T_i)):
                for j in range(len(P_i)):

                    row = i * len(P_i) + j

                    ax = axs4[row]

                    p = P_i[j]
                    t = T_i[i]

                    # fitted spectrum
                    xsec_fit = xaf.calculate_xsec(t, p, fit_coeffs)

                    plot_title = f"{p:.2f}$\,$Pa $-$ {t:.0f}$\,$K"

                    ylabel = '$a_{xsec}$ [m$^2$]'

                    if row == len(T_i) * len(P_i) - 1:
                        xlabel = 'Wavenumber [cm$^{-1}$]'
                    else:
                        xlabel = None

                    axs4[row] = xaf.plot_xsec(wvn, xsec_fit, [], ax, xlim=None, xlabel=xlabel, ylabel=ylabel,
                                              plot_title=plot_title)

            # %% edge cases T-derivative

            print('Plotting fig4.1')

            T_i = [np.min(MinT), np.max(MaxT)]
            P_i = [10., np.min(MinP), np.max(MaxP)]

            fig41, axs41 = xaf.default_figure(len(T_i) * len(P_i), 1, width_in_cm=20.9,
                                              height_in_cm=29.7, sharey=False)
            fig41.subplots_adjust(hspace=0.5)

            for i in range(len(T_i)):
                for j in range(len(P_i)):

                    row = i * len(P_i) + j

                    # axs=axs41[row]

                    p = P_i[j]
                    t = T_i[i]
                    deltaT = t * 1e-4

                    # T derivative of fitted spectrum
                    DxsecDT, _ = xaf.xsec_derivative(t, p, fit_coeffs)

                    # check derivative
                    DxsecDTnum = (xaf.calculate_xsec(t + deltaT, p, fit_coeffs)
                                  - xaf.calculate_xsec(t, p, fit_coeffs)) / deltaT

                    plot_title = f"{p:.2f}$\,$Pa $-$ {t:.0f}$\,$K"

                    ylabel = '$\partial_T a_{xsec}$ [m$^2\,$K$^{-1}$ ]'

                    if row == len(T_i) * len(P_i) - 1:
                        xlabel = 'Wavenumber [cm$^{-1}$]'
                    else:
                        xlabel = None

                    axs41[row] = xaf.plot_xsec(wvn, DxsecDT, DxsecDTnum, axs41[row], xlim=None, xlabel=xlabel,
                                               ylabel=ylabel,
                                               plot_title=plot_title)

            # %% edge cases P-derivative

            print('Plotting fig4.2')

            T_i = [np.min(MinT), np.max(MaxT)]
            P_i = [10., np.min(MinP), np.max(MaxP)]

            fig42, axs42 = xaf.default_figure(len(T_i) * len(P_i), 1, width_in_cm=20.9,
                                              height_in_cm=29.7, sharey=False)
            fig42.subplots_adjust(hspace=0.5)

            for i in range(len(T_i)):
                for j in range(len(P_i)):

                    row = i * len(P_i) + j

                    p = P_i[j]
                    t = T_i[i]
                    deltaP = p ** 1e-4

                    # P derivative of fitted spectrum
                    _, DxsecDP = xaf.xsec_derivative(t, p, fit_coeffs)

                    # check derivative
                    DxsecDPnum = (xaf.calculate_xsec(t, p + deltaP, fit_coeffs)
                                  - xaf.calculate_xsec(t, p, fit_coeffs)) / deltaP

                    plot_title = f"{p:.2f}$\,$Pa $-$ {t:.0f}$\,$K"

                    ylabel = '$\partial_P a_{xsec}$ [m$^2\,$Pa$^{-1}$ ]'

                    if row == len(T_i) * len(P_i) - 1:
                        xlabel = 'Wavenumber [cm$^{-1}$]'
                    else:
                        xlabel = None

                    axs41[row] = xaf.plot_xsec(wvn, DxsecDP, DxsecDPnum, axs42[row], xlim=None, xlabel=xlabel,
                                               ylabel=ylabel,
                                               plot_title=plot_title)

            # %% show full space

            print('Plotting fig5')

            Ptest = np.logspace(3, np.log10(101325), 100)
            Ttest = np.linspace(180, 300, 121)

            Xmean = np.zeros((len(Ptest), len(Ttest))) * np.nan

            for i in range(len(Ptest)):
                for j in range(len(Ttest)):
                    # %fitted spectrum
                    xsec_fit = xaf.calculate_xsec(Ttest[j], Ptest[i], fit_coeffs)

                    # mean over frequency of fitted spectrum
                    Xmean[i, j] = np.mean(xsec_fit)

            Xmean[Ptest < np.min(MinP), :] = np.nan
            Xmean[Ptest > np.max(MaxP), :] = np.nan

            Xmean_p = np.nanmean(Xmean, axis=0)
            Xmean_r = Xmean * np.nan

            # change in pressure relative to mean over pressure
            for j in range(len(Ttest)):
                Xmean_r[:, j] = Xmean[:, j] / Xmean_p[j] - 1.

            fig5, axs5 = xaf.default_figure(1, 2)

            xlabel = 'Temperature [K]'
            ylabel = 'Pressure [hPa]'
            cbar_label = '[m$^2$]'
            title = '$\\overline{a_{xsec,f}}$ (fit)'

            fig5, axs5[0], pcm, cbar = xaf.pcolor_plot(Ttest, Ptest / 100, Xmean, fig5, axs5[0],
                                                       np.nanmin(Xmean), np.nanmax(Xmean),
                                                       xlabel=xlabel, ylabel=ylabel,
                                                       cmap='temperature', title=title,
                                                       cbar_label=cbar_label)

            axs5[0].invert_yaxis()

            title = '$\\frac{\\overline{a_{xsec,f}} (p,T) } {< \\overline{a_{xsec,f}} >_p} -1 $'
            cbar_label = '[$\\%$]'
            fig5, axs5[1], pcm, cbar = xaf.pcolor_plot(Ttest, Ptest / 100, Xmean_r * 100, fig5, axs5[1],
                                                       -5., 5.,
                                                       xlabel=xlabel,
                                                       cmap='difference', title=title,
                                                       cbar_label=cbar_label)

            # %% some fit meta information

            print('Plotting fig6')

            fig6, axs6 = xaf.default_figure(6, 1, width_in_cm=20.9, height_in_cm=29.7, sharey='none')
            fig6.subplots_adjust(hspace=0.5)

            # Number of points per fit
            plot_title = 'Number of points per fit'
            ylabel = '$N_{obs}$ [$\\,$]'

            axs6[0] = xaf.plot_xsec(wvn, N_data, [], axs6[0], xlim=None, xlabel=None, ylabel=ylabel,
                                    plot_title=plot_title)

            # minimum temperature of fit
            plot_title = 'Minimum temperature of fit'
            ylabel = '$T_{min}$ [K]'
            axs6[1] = xaf.plot_xsec(wvn, MinT, [], axs6[1], xlim=None, xlabel=None, ylabel=ylabel,
                                    plot_title=plot_title)

            # maximum temperature of fit
            plot_title = 'Maximum temperature of fit'
            ylabel = '$T_{min}$ [K]'
            axs6[2] = xaf.plot_xsec(wvn, MaxT, [], axs6[2], xlim=None, xlabel=None, ylabel=ylabel,
                                    plot_title=plot_title)

            # minimum pressure of fit
            plot_title = 'Minimum pressure of fit'
            ylabel = '$p_{min}$ [hPa]'
            axs6[3] = xaf.plot_xsec(wvn, MinP / 100, [], axs6[3], xlim=None, xlabel=None, ylabel=ylabel,
                                    plot_title=plot_title)
            # axs6[3].set_yscale('log')

            # maximum temperature of fit
            plot_title = 'Maximum pressure of fit'
            ylabel = '$p_{min}$ [hPa]'
            xlabel = 'Wavenumber [cm$^{-1}$]'
            axs6[4] = xaf.plot_xsec(wvn, MaxP / 100, [], axs6[4], xlim=None, xlabel=xlabel, ylabel=ylabel,
                                    plot_title=plot_title)
            # axs6[4].set_yscale('log')

            plot_title = 'Goodness of fit'
            ylabel = '$R^2$'
            xlabel = 'Wavenumber [cm$^{-1}$]'
            axs6[5] = xaf.plot_xsec(wvn, R2, [], axs6[5], xlim=None, xlabel=xlabel, ylabel=ylabel,
                                    plot_title=plot_title)

            # %% save figures

            print('Saving figures')

            plotfolder = os.path.join(script_path, main_plot_folder, species_arts)

            if not os.path.exists(plotfolder):
                os.makedirs(plotfolder)

            band_str = '.band' + str(band_no)

            plotname0 = os.path.join(plotfolder, species_arts + band_str + '.fit_and_obs_spectra.pdf')
            fig0.savefig(plotname0)

            plotname1 = os.path.join(plotfolder, species_arts + band_str + '.delta_xsec_spectra.pdf')
            fig1.savefig(plotname1)

            plotname2 = os.path.join(plotfolder, species_arts + band_str + '.stat_overview.pdf')
            fig2.savefig(plotname2)

            plotname3 = os.path.join(plotfolder, species_arts + band_str + '.fit_coeffs.pdf')
            fig3.savefig(plotname3)

            plotname4 = os.path.join(plotfolder, species_arts + band_str + '.example_edge_case.pdf')
            fig4.savefig(plotname4)

            plotname41 = os.path.join(plotfolder, species_arts + band_str + '.example_edge_case.Tderivative.pdf')
            fig41.savefig(plotname41)

            plotname42 = os.path.join(plotfolder, species_arts + band_str + '.example_edge_case.Pderivative.pdf')
            fig42.savefig(plotname42)

            plotname5 = os.path.join(plotfolder, species_arts + band_str + '.change_in_p.pdf')
            fig5.savefig(plotname5)

            plotname6 = os.path.join(plotfolder, species_arts + band_str + '.fit_limits.pdf')
            fig6.savefig(plotname6)

            xaf.plt.close('all')

    # %%

    if store_coeffs == True:
        xsec_record = pyarts.classes.XsecRecord()
        xsec_record.spec=species_arts
        xsec_record.fitminpressures=min_pressures
        xsec_record.fitmaxpressures=max_pressures
        xsec_record.fitmintemperatures=min_temperatures
        xsec_record.fitmaxtemperatures=max_temperatures
        xsec_record.fitcoeffs=Xsec_processed_data_array
        xsec_record.version=2


        coeff_folder = os.path.join(script_path, coeff_folder)
        fid = '.xml'

        if not os.path.exists(coeff_folder):
            os.makedirs(coeff_folder)

        coeff_file_name = os.path.join(coeff_folder, species_arts + fid)

        pyarts.xml.save(xsec_record, coeff_file_name, precision='.14e', format='binary')

        print('Saving coefficients')


# %%

if __name__ == '__main__':

    # show plots?
    plotting = True

    # store coefficients?
    store_coeffs = True

    # folder of harmonized data
    harmonized_folder = '../data/harmonized_data/'

    # main plot folder
    main_plot_folder = '../plots/'

    # coefficients folder
    coeff_folder = '../coefficients/'

    script_path = os.getcwd()

    filelist = glob.glob(script_path + '/' + harmonized_folder + '*.json.gz')
    filelist.sort()

    # Get species
    all_species = []
    for file in filelist:

        filename = os.path.basename(file)
        species = filename.split('.')[0]

        if species not in all_species:
            all_species.append(species)


    for species in all_species:
        process_xsec_coefficients(species, harmonized_folder, coeff_folder, main_plot_folder, store_coeffs=store_coeffs,
                                  plotting=plotting)
