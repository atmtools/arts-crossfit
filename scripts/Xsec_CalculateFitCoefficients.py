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
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from gzip import GzipFile

import numpy as np
import pyarts
from typhon import constants

import Xsec_aux_functions as xaf


# %%

def process_xsec_coefficients(species, harmonized_folder, coeff_folder, main_plot_folder, store_coeffs=True,
                              plotting=True):
    # get the files
    filelist = glob.glob(harmonized_folder + species + '.*.json.gz')
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

    # %% Do the Fitting

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
        fit_coeffs = np.zeros((5, len(wvn)))
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

            fit_result = xaf.fit_xsec_data(T, P, Xsec, min_deltaP=80000, min_deltaT=40.)

            fit_coeffs[:, idx] = fit_result['coefficients']
            N_data[idx] = fit_result['NumberOfPoints']
            MinP[idx] = fit_result['MinP']
            MaxP[idx] = fit_result['MaxP']
            MinT[idx] = fit_result['MinT']
            MaxT[idx] = fit_result['MaxT']
            fit_goodness[idx] = fit_result['sum of residuals']

        # %% prepare output

        # frequency in Hz
        freq = wvn * constants.c * 100

        s_data_temp = pyarts.classes.GriddedField2()
        s_data_temp.gridnames = ['frequency grid [Hz]', 'fit coefficients [m]']
        s_data_temp.grids = [pyarts.classes.Vector(freq),
                             pyarts.classes.ArrayOfString(['p00', 'p10', 'p01', 'p20', 'p02'])]
        s_data_temp.data = fit_coeffs.transpose()
        s_data_temp.name = (species + '-band_' + str(band_no))

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
            L_raw = np.zeros(np.shape(T)) * np.nan
            dw_raw = np.zeros(np.shape(T)) * np.nan


            for i in range(len(RMSE)):
                # fitted spectrum
                XsecTestFit[i, :] = xaf.calculate_xsec_fullmodel(T[i], P[i], fit_coeffs)

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

                # bandwidth of the hitran spectrum
                L_raw[i]=(data[i]['wmax_rawdata']-data[i]['wmin_rawdata'])/(wvn.max()-wvn.min())

                #wavenumber resolution of original hitran spectrum
                dw_raw[i]=data[i]['DeltaWvnOfRawdata']

                deltaXsec=(XsecTestFit[i, :] - XsecTest[i, :])
                logic_delta=np.isnan(deltaXsec)
                N_not_nan=np.sum(~logic_delta)
                deltaXsec[logic_delta]=0.
                dw=np.nanmean(np.diff(wvn))
                RMSE[i] = np.sqrt(np.trapz( deltaXsec**2,wvn)/dw/N_not_nan)*1e4
                bias[i] = np.trapz(deltaXsec,wvn)/dw/N_not_nan*1e4
                StDev[i] = np.nanstd(XsecTest[i, :])
                XsecMean[i] = np.nanmean(XsecTest[i, :])
                XsecSum[i] = np.nansum(XsecTest[i, :])
                XsecMeanFit[i] = np.nanmean(XsecTestFit[i, :])
                XsecSumFit[i] = np.nansum(XsecTestFit[i, :])

            # Calculate Goodnes of fit
            for i in range(np.size(fit_coeffs, axis=1)):
                R2[i] = xaf.calc_Rsquare(XsecTest[:, i], XsecTestFit[:, i], np.sum(fit_coeffs[:, i] > 0))

            # %% some preparation for  plotting

            # # make random selection of data to show
            # nop = 9
            # index = np.random.permutation(len(T))
            # index = index[0:(np.min([len(data), nop]))]
            idx=np.argsort(T)

            #number of plots
            nop=9;
            if len(P)<nop:
                index=idx
            else:
                selection=np.round(np.linspace(0,len(idx)-1,nop))
                selection=selection.astype(int)
                index=idx[selection]


            print('Start with plotting!')

            # %% plot xsec for some points

            print('Plotting fig0 ')

            fig0, axs0 = xaf.default_figure(3, 3)

            for i in range(len(index)):

                row = int(i / 3)
                col = i - (row * 3)

                ax = axs0[row, col]

                p = P[index[i]]
                t = T[index[i]]

                # fitted spectrum
                # xsec_fit = xaf.calculate_xsec(t, p, fit_coeffs, density_flag=density_flag)
                xsec_fit = xaf.calculate_xsec_fullmodel(t, p, fit_coeffs)
                xsec_fit=xsec_fit*1e4

                # observed
                xsec = np.array(data[index[i]]['xsec'])*1e4

                plot_title = f"{p:.2f}$\,$Pa $-$ {t:.0f}$\,$K"

                if col == 0:
                    ylabel = '$a_{xsec}$ [cm$^2$]'
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

                p = P[index[i]]
                t = T[index[i]]

                # fitted spectrum
                xsec_fit = xaf.calculate_xsec(t, p, fit_coeffs)

                # observed
                xsec = data[index[i]]['xsec']

                # difference
                dxsec = (xsec - xsec_fit)*1e4

                plot_title = f"{p:.2f}$\,$Pa $-$ {t:.0f}$\,$K"

                if col == 0:
                    ylabel = '$\Delta a_{xsec}$ [cm$^2$]'
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

            z11 = RMSE 
            fig2, axs2[0, 0] = xaf.scatter_plot(T, P, z11, fig2, axs2[0, 0], #clim=[0, 40],
                                                plot_title='$RMSE$',
                                                cbar_label='[cm$^2$]')
            axs2[0, 0].invert_yaxis()

            z12 = bias 
            max_abs_bias=np.max(np.abs(bias))
            fig2, axs2[0, 1] = xaf.scatter_plot(T, P, z12, fig2, axs2[0, 1], clim=[-max_abs_bias, max_abs_bias],
                                                plot_title='$bias$',
                                                cbar_label='[cm$^2$]',
                                                cmap='difference')

            z13 = StDev*1e4
            fig2, axs2[0, 2] = xaf.scatter_plot(T, P, z13, fig2, axs2[0, 2],
                                                plot_title='std($a_{xsec,obs}$)',
                                                cbar_label='[cm$^2$]')

            z21 = XsecInt
            z21_max = np.max([XsecInt.max(), XsecIntFitFull.max()])
            z21_min = np.min([XsecInt.min(), XsecIntFitFull.min()])
            fig2, axs2[1, 0] = xaf.scatter_plot(T, P, z21, fig2, axs2[1, 0], clim=[z21_min, z21_max],
                                                plot_title='$\int a_{xsec,obs,raw}$',
                                                cbar_label='[cm$^2$ cm$^{-1}$]',
                                                cmap='temperature')

            z22 = XsecIntFit
            fig2, axs2[1, 1] = xaf.scatter_plot(T, P, z22, fig2, axs2[1, 1], clim=[z21_min, z21_max],
                                                plot_title='$\int a_{xsec,fit}$',
                                                cbar_label='[cm$^2$ cm$^{-1}$]',
                                                cmap='temperature')

            z23 = XsecIntFitFull
            fig2, axs2[1, 2] = xaf.scatter_plot(T, P, z23, fig2, axs2[1, 2], clim=[z21_min, z21_max],
                                                plot_title='$\int a_{xsec,fit,full}$',
                                                cbar_label='[cm$^2$ cm$^{-1}$]',
                                                cmap='temperature')

            z31 = L_raw*100
            fig2, axs2[2, 0] = xaf.scatter_plot(T, P, z31, fig2, axs2[2, 0],
                                                clim=[0, 100],
                                                plot_title='Overlap between Hitran and band',
                                                cbar_label='[$\\%$]')

            z32 = (XsecIntFit / XsecInt - 1) * 100
            z32_limit = np.max([abs(z32.min()), abs(z32.max())])
            fig2, axs2[2, 1] = xaf.scatter_plot(T, P, z32, fig2, axs2[2, 1], clim=[-z32_limit, z32_limit],
                                                plot_title='$(\int a_{xsec,fit}/\int a_{xsec,obs,raw})-1$',
                                                cbar_label='[$\\%$]', cmap='difference')

            z33 = dw_raw
            z33_limit = z33.max()
            fig2, axs2[2, 2] = xaf.scatter_plot(T, P, z33, fig2, axs2[2, 2], clim=[0, z33_limit],
                                                plot_title='wavenumber resolution (Hitran)',
                                                cbar_label='[$cm^{-1}$]')

            # %% plot coefficients

            print('Plotting fig3')

            fig3, axs3 = xaf.default_figure(5, 1, sharey='none', sharex='all',
                                            width_in_cm=20.9, height_in_cm=29.7)

            coeff_names = fit_result['coeff_names']
            formula = fit_result['formula']

            for i in range(np.size(fit_coeffs, axis=0)):

                ax = axs3[i]

                plot_title = coeff_names[i]

                ylabel = '[m]'

                if i == np.size(fit_coeffs, axis=0)-1:
                    xlabel = 'Wavenumber [cm$^{-1}$]'
                else:
                    xlabel = None

                axs3[i] = xaf.plot_xsec(wvn, fit_coeffs[i, :], [], ax, xlim=None, xlabel=xlabel, ylabel=ylabel,
                                               plot_title=plot_title)


            sup_text = species + ': $a_{xsec}=' + formula + '$; $x=T/T_0$, $y=p/p_0$'

            fig3.suptitle(sup_text)

            # %% edge cases

            print('Plotting fig4')

            T_i = [np.min(MinT) * 0.9, np.max(MaxT) / 0.9]
            P_i = [10., np.min(MinP), np.max(MaxP)]

            fig4, axs4 = xaf.default_figure(len(T_i) * len(P_i), 1,
                                            width_in_cm=20.9, height_in_cm=29.7)
            fig4.subplots_adjust(hspace=0.5)

            for i in range(len(T_i)):
                for j in range(len(P_i)):

                    row = i * len(P_i) + j

                    ax = axs4[row]

                    p = P_i[j]
                    t = T_i[i]

                    # fitted spectrum
                    xsec_fit = xaf.calculate_xsec_fullmodel(t, p, fit_coeffs)

                    plot_title = f"{p:.2f}$\,$Pa $-$ {t:.0f}$\,$K"

                    ylabel = '$a_{xsec}$ [cm$^2$]'

                    if row == len(T_i) * len(P_i) - 1:
                        xlabel = 'Wavenumber [cm$^{-1}$]'
                    else:
                        xlabel = None

                    axs4[row] = xaf.plot_xsec(wvn, xsec_fit*1e4, [], ax, xlim=None, xlabel=xlabel, ylabel=ylabel,
                                              plot_title=plot_title)

            # %% edge cases T-derivative

            print('Plotting fig4.1')

            T_i = [np.min(MinT), np.max(MaxT)]
            P_i = [np.max([1., np.min(MinP)]), np.max([np.max(MaxP), 101325.])]

            fig41, axs41 = xaf.default_figure(len(T_i) * len(P_i), 1, width_in_cm=20.9,
                                              height_in_cm=29.7, sharey=False)
            fig41.subplots_adjust(hspace=0.5)

            for i in range(len(T_i)):
                for j in range(len(P_i)):

                    row = i * len(P_i) + j

                    # axs=axs41[row]

                    p = P_i[j]
                    t = T_i[i]
                    deltaT = t * 1e-5

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

            fig42, axs42 = xaf.default_figure(len(T_i) * len(P_i), 1, width_in_cm=20.9,
                                              height_in_cm=29.7, sharey=False)
            fig42.subplots_adjust(hspace=0.5)

            for i in range(len(T_i)):
                for j in range(len(P_i)):

                    row = i * len(P_i) + j

                    p = P_i[j]
                    t = T_i[i]
                    deltaP = p * 1e-5

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

                    axs42[row] = xaf.plot_xsec(wvn, DxsecDP, DxsecDPnum, axs42[row], xlim=None, xlabel=xlabel,
                                               ylabel=ylabel,
                                               plot_title=plot_title)

            # %% show full space

            print('Plotting fig5')

            Ptest = np.logspace(1, np.log10(101325), 200)
            Ttest = np.linspace(175, 325, 151)

            Xmean = np.zeros((len(Ptest), len(Ttest))) * np.nan

            for i in range(len(Ptest)):
                for j in range(len(Ttest)):
                    # %fitted spectrum
                    xsec_fit = xaf.calculate_xsec_fullmodel(Ttest[j], Ptest[i],
                                                            fit_coeffs)

                    # Integral over frequency of fitted spectrum
                    Xmean[i, j] = np.trapz(xsec_fit, wvn) * 1e4

            fig5, axs5 = xaf.default_figure(1, 2, sharey=False)

            xlabel = 'Temperature [K]'
            ylabel = 'Pressure [hPa]'
            # cbar_label = '[m$^2$]'
            cbar_label = '[cm$^2$ cm$^{-1}$]'
            title = r'$\int a_{xsec} $d$\nu$ (fit)'

            fig5, axs5[0], pcm, cbar = xaf.pcolor_plot(Ttest, Ptest / 100, Xmean, fig5, axs5[0],
                                                    np.percentile(Xmean, 1), np.percentile(Xmean, 99),
                                                    xlabel=xlabel, ylabel=ylabel,
                                                    cmap='temperature', title=title,
                                                    cbar_label=cbar_label)

            axs5[0].scatter(T, P / 100, 50, XsecInt, cmap='temperature', #edgecolors='w',
                         vmin=cbar.vmin, vmax=cbar.vmax, zorder=1e11)

            axs5[0].invert_yaxis()


            #plot plot xsec_int as fcn of T for surface and 100 hPa
            P_surf,idx_surf=xaf.find_nearest(Ptest, 101325.)
            P_100hPa,idx_100hPa=xaf.find_nearest(Ptest, 10000.)

            #logical to select xsec from raw data
            logic_surf=np.abs(P-101325)<=101325*0.10
            logic_100hPa=np.abs(P-5000)<=5000


            axs5[1] = xaf.plot_xsec(Ttest, Xmean[idx_surf,:],Xmean[idx_100hPa,:],
                                    axs5[1],
                                    xlabel=xlabel,
                                    ylabel=cbar_label,
                                    plot_title=title,
                                    legend=True,
                                    labels=['1013.25hPa','50hPa'],
                                    linewidth=1.5)

            colors=xaf.cmap_matlab_lines()
            axs5[1].plot(T[logic_surf],XsecInt[logic_surf],marker='o',
                         linestyle='None',color=colors[0,:],label='Hitran,surf')
            axs5[1].plot(T[logic_100hPa],XsecInt[logic_100hPa],marker='+',
                         linestyle='None',color=colors[1,:],label='Hitran,50hPa$\pm$50hPa')

            axs5[1].yaxis.tick_right()
            axs5[1].yaxis.set_label_position("right")
            axs5[1].legend()


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
            axs6[4] = xaf.plot_xsec(wvn, MaxP / 100, [], axs6[4], xlim=None, xlabel=xlabel, ylabel=ylabel,
                                    plot_title=plot_title)
            # axs6[4].set_yscale('log')

            plot_title = 'Goodness of fit'
            ylabel = '$R^2$'
            xlabel = 'Wavenumber [cm$^{-1}$]'
            axs6[5] = xaf.plot_xsec(wvn, R2, [], axs6[5], xlim=None, xlabel=xlabel, ylabel=ylabel,
                                    plot_title=plot_title)
            
            # %% show p-T dependency for some wavenumbers

            print('Plotting fig7')
            
            #get frequency indices at nop positions uniformly distributed in amplitude
            xsec_test = xaf.calculate_xsec_fullmodel(293.15, 101325.,fit_coeffs)
            index_sort=np.argsort(xsec_test)            
            index_temp=np.asanyarray(np.linspace(0,len(fit_coeffs[0,:])-1,nop),dtype=int)
            index_wvn=index_sort[index_temp]
            
            
            #get data for these nine points
            Xsec_nu=np.zeros((len(P),len(index_wvn)))
            
            for i in range(len(P)):
                Xsec_nu[i,:]=np.array(data[i]['xsec'])[index_wvn]*1e4
            
            
            Xsec_nu_fit = np.zeros((len(Ptest), len(Ttest), len(index_wvn))) * np.nan
            
            for ii in range(len(Ptest)):
                for jj in range(len(Ttest)):

                    # %fitted spectrum
                    xsec_fit = xaf.calculate_xsec_fullmodel(Ttest[jj], Ptest[ii],
                                                                fit_coeffs[:,index_wvn])

                    # Fitted Xsec in cm^-2
                    Xsec_nu_fit[ii, jj, :] = xsec_fit * 1e4
                    
            fig7, axs7 = xaf.default_figure(3, 3)                                          

            for i in range(len(index_wvn)):
            
                row = int(i / 3)
                col = i - (row * 3)
            
                ax = axs7[row, col]
                
                if col == 0:
                    ylabel = 'Pressure [hPa]'
                else:
                    ylabel = None

                if row == 2:
                    xlabel = 'Temperature [K]'
                else:
                    xlabel = None
                
                cbar_label = '[cm$^2$]'
                title = '$a_{xsec}$' + f'({wvn[index_wvn[i]]:.3f}' + 'cm$^{-1}$)'
    
                fig7, ax, pcm, cbar = xaf.pcolor_plot(Ttest, Ptest / 100, Xsec_nu_fit[:,:,i], fig7, ax,
                                                        np.percentile(Xsec_nu_fit[:,:,i], 1), np.percentile(Xsec_nu_fit[:,:,i], 99),
                                                        xlabel=xlabel, ylabel=ylabel,
                                                        cmap='temperature', title=title,
                                                        cbar_label=cbar_label)
    
                ax.scatter(T, P / 100, 50, Xsec_nu[:,i], cmap='temperature', edgecolors='w',
                             vmin=cbar.vmin, vmax=cbar.vmax, zorder=1e11)
    
                ax.invert_yaxis()

            # %% save figures

            print('Saving figures')

            band_str = 'band_' + str(band_no)

            plotfolder = os.path.join(main_plot_folder, species_arts, band_str)

            if not os.path.exists(plotfolder):
                os.makedirs(plotfolder)

            plotname0 = os.path.join(plotfolder, species_arts + '.fit_and_obs_spectra.pdf')
            fig0.savefig(plotname0)

            plotname1 = os.path.join(plotfolder, species_arts + '.delta_xsec_spectra.pdf')
            fig1.savefig(plotname1)

            plotname2 = os.path.join(plotfolder, species_arts + '.stat_overview.pdf')
            fig2.savefig(plotname2)

            plotname3 = os.path.join(plotfolder, species_arts + '.fit_coeffs.pdf')
            fig3.savefig(plotname3)

            plotname4 = os.path.join(plotfolder, species_arts + '.example_edge_case.pdf')
            fig4.savefig(plotname4)

            plotname41 = os.path.join(plotfolder, species_arts + '.example_edge_case.Tderivative.pdf')
            fig41.savefig(plotname41)

            plotname42 = os.path.join(plotfolder, species_arts + '.example_edge_case.Pderivative.pdf')
            fig42.savefig(plotname42)

            plotname5 = os.path.join(plotfolder, species_arts + '.change_in_p.pdf')
            fig5.savefig(plotname5)

            plotname6 = os.path.join(plotfolder, species_arts + '.fit_limits.pdf')
            fig6.savefig(plotname6)
            
            plotname7 = os.path.join(plotfolder, species_arts + '.Xsec_nu_as_fcn_of_TP.pdf')
            fig7.savefig(plotname7)

            xaf.plt.close('all')

    # %%

    if store_coeffs == True:
        xsec_record = pyarts.classes.XsecRecord()
        xsec_record.spec = species_arts
        xsec_record.fitminpressures = min_pressures
        xsec_record.fitmaxpressures = max_pressures
        xsec_record.fitmintemperatures = min_temperatures
        xsec_record.fitmaxtemperatures = max_temperatures
        xsec_record.fitcoeffs = Xsec_processed_data_array
        xsec_record.version = 2

        fid = '.xml'

        if not os.path.exists(coeff_folder):
            os.makedirs(coeff_folder)

        coeff_file_name = os.path.join(coeff_folder, species_arts + fid)

        pyarts.xml.save(xsec_record, coeff_file_name, precision='.14e', format='binary')

        print('Saving coefficients')


def parse_args():
    """Parse commandline arguments"""
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p",
                        "--plots",
                        action="store_true",
                        help="Generate diagnostic plots.")
    parser.add_argument("-n",
                        "--dry-run",
                        action="store_true",
                        help="Don't store coefficients.")

    return parser.parse_args()


# %%
if __name__ == '__main__':
    args = parse_args()

    # show plots?
    plotting = args.plots

    # store coefficients?
    store_coeffs = not args.dry_run

    script_path = os.path.dirname(os.path.realpath(__file__))

    # folder of harmonized data
    harmonized_folder = os.path.join(script_path, '../data/harmonized_data/')

    # main plot folder
    main_plot_folder = os.path.join(script_path, '../plots/')

    # coefficients folder
    coeff_folder = os.path.join(script_path, '../coefficients/')

    filelist = glob.glob(harmonized_folder + '*.json.gz')
    filelist.sort()

    # Get species
    all_species = []
    for file in filelist:

        filename = os.path.basename(file)
        species = filename.split('.')[0]

        if species not in all_species:
            all_species.append(species)

    # all_species=[all_species[31]]

    for species in all_species:
        process_xsec_coefficients(species, harmonized_folder, coeff_folder, main_plot_folder, store_coeffs=store_coeffs,
                                  plotting=plotting)
