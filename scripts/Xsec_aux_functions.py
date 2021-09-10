#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 16:10:11 2020

@author: Manfred Brath

This file contains the functions that are needed for the harmonization of the Hitran
absorption cross section data and for the calculations of the fit coefficients.
"""
from glob import glob
from os import path

import matplotlib.patches as ptch
import matplotlib.pyplot as plt
import numpy as np
import pyarts
from matplotlib.font_manager import FontProperties
from scipy.interpolate import interp1d
from scipy.linalg import lstsq

# %% constants

# speed of light
c0 = 299792458.0  # [m/s]


# %% fit related functions


def fit_poly21(xdata, ydata, zdata):
    '''
    2d semi quadratic fit:
    z = p00 + p10*x + p01*y + p20*x**2

    Args:
        xdata:  vector
                independent data.
        ydata:  vector
                independent data.
        zdata:  vector
                data, which depends on xdata and ydata.

    Returns:
        poly:   vector
                coefficients of fit, see above for the order.
        res:    float
                summed residuums.
        rnk:    int
                Effective rank of design matrix M.
        s:      ndarray or None
                Singular values of M.

    '''

    M = np.ones((len(xdata), 4))
    M[:, 1] = xdata  # p01
    M[:, 2] = ydata  # p10
    M[:, 3] = xdata ** 2  # p20

    poly, res, rnk, s = lstsq(M, zdata)

    return poly, res, rnk, s


def fit_poly11(xdata, ydata, zdata):
    '''
    2d linear fit:
    z = p00 + p10*x + p01*y

    Args:
        xdata:  vector
                independent data.
        ydata:  vector
                independent data.
        zdata:  vector
                data, which depends on xdata and ydata.

    Returns:
        poly:   vector
                coefficients of fit, see above for the order.
        res:    float
                summed residuums.
        rnk:    int
                Effective rank of design matrix M.
        s:      ndarray or None
                Singular values of M.

    '''

    M = np.ones((len(xdata), 3))
    M[:, 1] = xdata  # p10
    M[:, 2] = ydata  # p01

    poly, res, rnk, s = lstsq(M, zdata)

    return poly, res, rnk, s


def fit_poly2(xdata, zdata):
    '''
    1d quadratic fit:
    z = p0 + p1*x + p2*x**2

    Args:
        xdata:  vector
                independent data.
        zdata:  vector
                data, which depends on xdata.

    Returns:
        poly:   vector
                coefficients of fit, see above for the order.
        res:    float
                summed residuums.
        rnk:    int
                Effective rank of design matrix M.
        s:      ndarray or None
                Singular values of M.

    '''

    # 1d quadratic fit:
    # z = p0 + p1*x + p2*x**2

    M = np.ones((len(xdata), 3))
    M[:, 1] = xdata  # p1
    M[:, 2] = xdata ** 2  # p2

    poly, res, rnk, s = lstsq(M, zdata)

    return poly, res, rnk, s


def fit_poly1(xdata, zdata):
    '''
    1d linear fit:
    z = p0 + p1*x

    Args:
        xdata:  vector
                independent data
        zdata:  vector
                data, which depends on xdata

    Returns:
        poly:   vector
                coefficients of fit, see above for the order
        res:    float
                summed residuums
        rnk:    int
                Effective rank of design matrix M
        s:      ndarray or None
                Singular values of M

    '''

    M = np.ones((len(xdata), 2))
    M[:, 1] = xdata  # p1

    poly, res, rnk, s = lstsq(M, zdata)

    return poly, res, rnk, s


def calc_Rsquare(y, yfit, Ncoeffs):
    '''
    calculates the adjusted R-square statistic
    Args:
        y:  vector
            true value.
        yfit: vector
            fited value.
        Ncoeffs: int
            number of fit coefficients.

    Returns:
        rsquare: float
            adjusted R-square statistic.

    '''

    Delta_y = y - yfit
    Var_y = y - np.nanmean(y)
    SSE = np.nansum(Delta_y ** 2)
    SST = np.nansum(Var_y ** 2)
    n = len(y)

    if SST == 0 or (n - Ncoeffs) == 0:
        SST = np.nan
        n = np.nan

    return 1. - SSE * (n - 1.) / (SST * (n - Ncoeffs))


def calculate_xsec(T, P, coeffs):
    '''
    Low level function to calculate the absorption cross section from the fitted
    coefficients

    Args:
        T: float
            Temperature in K.
        P: float
            Pressure in Pa.
        coeffs: matrix
            fit coefficients.

    Returns:
        Xsec: vector
            Absorption cross section in m**2.


    The fit model
    2d quadratic fit:
    z= p00 + p10*x + p01*y + p20*x**2
    z=Xsec
    x=T
    y=P

    coeffs[0,:]           p00
    coeffs[1,:]           p10
    coeffs[2,:]           p01
    coeffs[3,:]           p20
    '''

    # distinguish if we calculate xsec for a lot of frequencies
    if len(np.shape(coeffs)) > 1:
        poly = np.zeros(4)
        poly[0] = 1
        poly[1] = T
        poly[2] = P
        poly[3] = T ** 2

        # allocate
        Xsec = np.zeros(np.shape(coeffs))

        for i in range(4):
            Xsec[i, :] = coeffs[i, :] * poly[i]

    # or for a lot of states
    else:
        poly = np.zeros((4, len(T)))
        poly[0, :] = 1.
        poly[1, :] = T
        poly[2, :] = P
        poly[3, :] = T ** 2

        # allocate
        Xsec = np.zeros((len(coeffs), len(T)))

        for i in range(4):
            Xsec[i, :] = coeffs[i] * poly[i, :]

    Xsec = np.sum(Xsec, axis=0)

    return Xsec


def calculate_xsec_fullmodel(T, P, coeffs):
    '''
    Function to calculate the absorption cross section from the fitted
    coefficients including check for negative values.

    Args:
        T: float
            Temperature in K.
        P: float
            Pressure in Pa.
        coeffs: matrix
            fit coefficients.

    Returns:
        Xsec: vector
            Absorption cross section in m**2.


    The fit model
    2d quadratic fit:
    z= p00 + p10*x + p01*y + p20*x**2

    z=Xsec
    x=T
    y=P

    coeffs[0,:]           p00
    coeffs[1,:]           p10
    coeffs[2,:]           p01
    coeffs[3,:]           p20
    '''

    # calculate raw xsecs
    xsec = calculate_xsec(T, P, coeffs)

    # Check for negative values and remove them without introducing bias, meaning
    # the integral over the spectrum must not change.
    logic = xsec < 0
    if np.sum(logic) > 0:

        # original sum over spectrum
        sumX_org = np.sum(xsec)

        # remove negative values
        xsec[logic] = 0

        if sumX_org >= 0:
            # estimate ratio between altered and original sum of spectrum
            w = sumX_org / np.sum(xsec)

            # scale altered spectrum
            xsec = xsec * w

    return xsec


def xsec_derivative(T, P, coeffs):
    '''
    Fucntion to calculate the derivative of the absorption cross section
    from the fitted coefficients.

    Args:
        T: float
            Temperature in K.
        P: float
            Pressure in Pa.
        coeffs: matrix
            fit coefficients.

    Returns:
        DxsecDT: vector
            Temperature derivative.
        DxsecDp: vector
            Pressure derivative.

    The fit model
    2d quadratic fit:
    z= p00 + p10*x + p01*y + p20*x**2

    z=Xsec
    x=T
    y=P

    coeffs[0,:]           p00
    coeffs[1,:]           p10
    coeffs[2,:]           p01
    coeffs[3,:]           p20

    '''

    # p00 = coeffs[0, :]
    p10 = coeffs[1, :]
    p01 = coeffs[2, :]
    p20 = coeffs[3, :]

    DxsecDT = p10 + 2 * p20 * T

    DxsecDp = p01

    return DxsecDT, DxsecDp


def fit_xsec_data(T, P, Xsec, min_deltaP=80000, min_deltaT=40., cnt_limit=2, k_outlier=1.5):
    '''
    FUnction to calculate the fit of the xsec at an arbitrary frequency

    Args:
        T: vector
            temperatures.
        P: vector
            pressures same length as `T`.
        Xsec: vector
            cross section same length as `T`.
        min_deltaP: float, optional
            minimum variability of sqrt(`P`) for fit. Defaults to 100
        min_deltaT: float, optional
            minimum variability of `T` for fit. Defaults to 20.
        cnt_limit:  integer, optional
            maximum number of iteration of the fit due to outlier removal.
        k_outlier:  float, optional
            scaling factor for outlier detection

    Returns:
        fit_result: dictionary
            results of the fit.

    The fit model
    2d quadratic fit:
    z= p00 + p10*x + p01*y + p20*x**2

    z=Xsec
    x=T
    y=P

    coeffs[0,:]           p00
    coeffs[1,:]           p10
    coeffs[2,:]           p01
    coeffs[3,:]           p20
    '''

    # FofP = P

    # FofXsec = Xsec

    # check for bad values
    logic_inf = np.isinf(T) | np.isinf(P) | np.isinf(Xsec)
    logic_nan = np.isnan(T) | np.isnan(P) | np.isnan(Xsec)
    logic_bad = logic_inf | logic_nan

    if np.sum(logic_bad) < len(T):

        # remove bad values
        xData = T[~logic_bad]
        yData = P[~logic_bad]
        zData = Xsec[~logic_bad]

        # get number of unique temperatures and pressures
        N_Tunique = np.size(np.unique(xData))
        N_Punique = np.size(np.unique(yData))

        # get some information about the distribution of data
        Ndata = np.sum(~logic_bad)

        Delta_P = max(yData) - min(yData)
        Delta_T = max(xData) - min(xData)

        cnt = 0
        while cnt < cnt_limit:

            # quadratic fit in temperature and linear in pressure
            if (Delta_P >= min_deltaP and Delta_T > 2 * min_deltaT and Ndata > 5
                    and N_Tunique > 4 and N_Punique > 1):

                p, res, rnk, s = fit_poly21(xData, yData, zData)

                coeffs = np.zeros(4)
                coeffs[0] = p[0]
                coeffs[1] = p[1]
                coeffs[2] = p[2]
                coeffs[3] = p[3]

            # linear fit in temperature and pressure
            elif (Delta_P >= min_deltaP and Delta_T > min_deltaT and Ndata > 3
                  and N_Tunique > 1 and N_Punique > 1):

                p, res, rnk, s = fit_poly11(xData, yData, zData)

                coeffs = np.zeros(4)
                coeffs[0] = p[0]
                coeffs[1] = p[1]
                coeffs[2] = p[2]

            # quadratic fit in temperature
            elif Delta_T > 2 * min_deltaT and N_Tunique > 4 and N_Punique == 1:

                p, res, rnk, s = fit_poly2(xData, zData)

                coeffs = np.zeros(4)
                coeffs[0] = p[0]
                coeffs[1] = p[1]
                coeffs[3] = p[2]

            # linear fit in temperature
            elif Delta_T > min_deltaT and N_Tunique > 2:
                p, res, rnk, s = fit_poly1(xData, zData)

                coeffs = np.zeros(4)
                coeffs[0] = p[0]
                coeffs[1] = p[1]

            # linear fit in pressure
            elif Delta_P > min_deltaP and N_Punique > 2:

                p, res, rnk, s = fit_poly1(yData, zData)

                coeffs = np.zeros(4)
                coeffs[0] = p[0]
                coeffs[2] = p[1]

            # no fit, just median value
            else:
                coeffs = np.zeros(4)
                coeffs[0] = np.median(zData)

                res = np.sum((zData - coeffs[0]) ** 2)
                rnk = np.nan
                s = np.nan

            if k_outlier > 0 or np.sum(coeffs == 0) >= 4:
                # Calculate residuals
                zData_fit = calculate_xsec(xData, yData, coeffs)

                residuals = zData_fit - zData

                # Check for outlier
                logic_out = np.logical_and(abs(residuals) > np.std(zData) * k_outlier, np.std(zData) > 0.)

                if np.sum(logic_out) > 0:
                    cnt += 1

                    if cnt < cnt_limit:
                        xData = xData[~logic_out]
                        yData = yData[~logic_out]
                        zData = zData[~logic_out]
                        Ndata = np.sum(~logic_out)
                else:
                    cnt = cnt_limit

        MinP = min(yData)
        MaxP = max(yData)

        MinT = min(xData)
        MaxT = max(xData)

        fit_result = {}
        fit_result['formula'] = 'p00 + p10*x + p01*y + p20*x**2'
        fit_result['coeff_names'] = ['p00', 'p10', 'p01', 'p20']
        fit_result['coefficients'] = coeffs
        fit_result['sum of residuals'] = res
        fit_result['rank'] = rnk
        fit_result['sgl_val'] = s
        fit_result['MinP'] = MinP
        fit_result['MaxP'] = MaxP
        fit_result['MinT'] = MinT
        fit_result['MaxT'] = MaxT
        fit_result['NumberOfPoints'] = Ndata
        # fit_result['R2']=R2

    else:
        fit_result = {}
        fit_result['formula'] = 'p00 + p10*x + p01*y + p20*x**2'
        fit_result['coeff_names'] = ['p00', 'p10', 'p01', 'p20']
        fit_result['coefficients'] = np.zeros(4)
        fit_result['sum of residuals'] = np.nan
        fit_result['rank'] = np.nan
        fit_result['sgl_val'] = np.nan
        fit_result['MinP'] = np.inf
        fit_result['MaxP'] = -np.inf
        fit_result['MinT'] = np.inf
        fit_result['MaxT'] = -np.inf
        fit_result['NumberOfPoints'] = 0

    return fit_result


# %% Apllication functions

def get_coeff_species(coefficients_folder):
    '''
    Convinience function that returns a list with the species inside the coefficients folders.

    Args:
        coefficients_folder (str): Path of the coefficients folder

    Returns:
        all_species (List): List with the Name of each species within the
                            coefficients folder

    '''

    filelist = glob(coefficients_folder + '*.xml.bin')
    filelist.sort()

    # Get species
    all_species = []
    for file in filelist:

        filename = path.basename(file)
        species = filename.split('.')[0]

        if species not in all_species:
            all_species.append(species)

    return all_species


def load_xsec_data(species, coeff_folder):
    '''
    Load the xsec data

    Args:
        species (str): Species name.
        coeff_folder (str): Path to the coefficient folder.

    Returns:
        xsec_data (XsecRecord): Xsec data.

    '''
    xsec_file = path.join(coeff_folder, f'{species}.xml')
    xsec_data = pyarts.xml.load(xsec_file)

    return xsec_data


def calculate_cross_sections(wvn_user, xsec_data, temperature=273., pressure=1013e2):
    '''
    Calculates absorption cross sections for desired wavenumbers.

    Args:
        wvn_user (Vector): Wavenumbers in [cm⁻¹].
        xsec_data (XsecRecord): Xsec data.
        temperature (Float, optional): Temperature in [K]. Default to 273..
        pressure (Float, optional): Pressure in [Pa]. Default to 1013e2.

    Returns:
        xsec_user (Vector): Absorption cross section in [m²].

    '''

    # convert desired wavenumber to frequency in [Hz]
    freq_user = wvn_user * c0 * 100

    # xsec_coeffs=xsec_data.fitcoeffs

    xsec_user = np.zeros(np.shape(wvn_user))

    for m in range(len(xsec_data.fitcoeffs)):
        # frequency of data in [Hz]
        freq_data = xsec_data.fitcoeffs[m].grids[0].data

        # fit coefficients of band m
        coeffs_m = xsec_data.fitcoeffs[m].data.data.transpose()

        # Calculate the cross section on their internal frequency grid
        xsec_temp = calculate_xsec_fullmodel(temperature, pressure, coeffs_m)

        # Interpolate cross sections to user grid
        f_int = interp1d(freq_data, xsec_temp, fill_value=0., bounds_error=False)
        xsec_user_m = f_int(freq_user)

        xsec_user = xsec_user + xsec_user_m

    return xsec_user


# %% aux function

def find_nearest(a, a0):
    '''
    Element in nd array `a` closest to the scalar value `a0`
    (Finds the needle in the haystack)

    Args:
        a (nd array): Haystack.
        a0 (scalar): needle

    Returns:
        found value (scalar): most similar needle
        index of found value

    '''

    idx = np.abs(a - a0).argmin()
    return a.flat[idx], idx


def getOverlap(a, b):
    '''
    Function to calculate the overlap between two ranges given
    by the edges of  each range.

    Args:
        a: vector
            Edges of range 1.
        b:  vector
            Edges of range 2.

    Returns:
        overlap: float
            overlap

    '''

    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def suggest_banddefinition(wvn_intervalls, dws):
    '''
    Function to define the band limits

    Args:
        wvn_intervalls (Matrix): Matrix with the wavenumber intervals of the
                                 set of spectra.
        dws (Vector): Vector with the wavenumber resolution for each set of
                      sprectra.

    Returns:
        band_limits_final (List): List of band limits.

    '''

    number_of_sets = len(dws)

    # sort the interval edges
    sorted_limits = np.sort(wvn_intervalls.flatten())

    band_limits = []
    band_dw = []

    # set the band limits according to the sorted limits, but only if there is
    # overlap with the data

    for i in range(np.size(sorted_limits) - 1):

        # temporary band limit
        band_limit_i = [sorted_limits[i], sorted_limits[i + 1]]

        # Check overlap with data
        ol_i = np.zeros(number_of_sets)

        for j in range(number_of_sets):
            ol_i[j] = getOverlap(wvn_intervalls[j, :], band_limit_i)

        if np.sum(ol_i) > 0:
            # get minimum dw per band (highest resolution)
            dw_i = np.min([dws[j] for j in range(number_of_sets) if ol_i[j] > 0])

            # Check for too small bands. A band needs at least 2 points.
            # This means dw_i must <= band_end-band_start
            if dw_i <= band_limit_i[1] - band_limit_i[0]:
                band_limits.append(band_limit_i)
                band_dw.append(dw_i)

    # now check if two adjacent bands have the same wavenumber resolution. If
    # two adjacent band have the same resoltion, we can join them.
    band_limits_final = []
    marker = np.zeros(len(band_limits))

    cnt = 0
    for i in range(len(band_limits) - 1):

        adjacent = (band_limits[i][1] - band_limits[i + 1][0] == 0.)
        same_dw = (band_dw[i] == band_dw[i + 1])

        crit = adjacent and same_dw

        if not crit:
            cnt += 1

        marker[i + 1] = cnt

    for i in range(cnt + 1):
        temp = [band_limits[j] for j in range(len(band_limits)) if marker[j] == i]

        band_limits_final.append([np.min(temp), np.max(temp)])

    return band_limits_final


# %% plotting routines

def cmap_matlab_lines():
    '''

    Returns:
        cmap: matrix
            Color  map with matlab like line colors.

    '''

    cmap = np.array([[0, 0.44701, 0.74101, 1],
                     [0.85001, 0.32501, 0.09801, 1],
                     [0.92901, 0.69401, 0.12501, 1],
                     [0.49401, 0.18401, 0.55601, 1],
                     [0.46601, 0.67401, 0.18801, 1],
                     [0.30101, 0.74501, 0.93301, 1],
                     [0.63501, 0.07801, 0.18401, 1]])

    return cmap


def default_figure(rows, columns, width_in_cm=29.7, height_in_cm=20.9,
                   sharey='all', sharex='all', dpi=150):
    '''
    simple function to define basic properties of a figure

    Args:
        rows: int
            rows of plots/axis.
        columns: int
            columns of plots/axis.
        width_in_cm: float
            figure width in cm.
        height_in_cm: float
            figure height in cm.
        sharey: str
            marker which y-axis are shared.
        sharex: str
            marker which x-axis are shared.
        dpi: float
            resolution for inline plots.

    Returns:
        fig: matplotlib figure object
            figure object.
        ax: matplotlib axis object or ndarray of axis objects
            matplotlib axis object or ndarray of axis objects.

    '''

    fig, ax = plt.subplots(rows, columns, sharey=sharey, sharex=sharex, dpi=dpi)
    fig.set_size_inches(width_in_cm / 2.54, h=height_in_cm / 2.54)

    return fig, ax


def set_tick_font(ax, font_name):
    '''
    Function to set tick font of x- and y-axis

    Args:
        ax: matplotlib axis object
            axis object.

        font_name: str
            font name.

    Returns:

    '''

    for tick in ax.get_xticklabels():
        tick.set_fontname(font_name)

    for tick in ax.get_yticklabels():
        tick.set_fontname(font_name)


def default_plot_format(ax, font_name=None):
    '''
    simple function to define basic properties of a plot

    Args:
        ax: matplotlib axis object
            axis object

        font_name: str
            font name

    Returns:
        ax: matplotlib axis object
            axis object

        font: font properties object
            font properties

    '''

    font = FontProperties()
    if font_name is not None:
        font.set_name(font_name)

    ax.set_prop_cycle(color=cmap_matlab_lines())

    ax.grid(which='both', linestyle=':', linewidth=0.25)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axes.tick_params(direction='in', which='both')

    return ax, font


def plot_xsec(wvn, Xsec, XsecFit, ax, xlim=None, xlabel=None, ylabel=None,
              plot_title=None, legend=False, font_name=None, labels=['obs', 'fit'],
              linewidth=0.25, fontsize=10, formatter=True):
    '''
    Wrapper to plot up to two cross sections in a plot. If only one cross section
    should be plotted set XsecFit to an empty list.

    Args:
        wvn: vector
            wavenumbers.
        Xsec: vector
            cross sections 1.
        XsecFit:
            cross sections 2.
        ax: matplotlib axis object
            axis object.
        xlim: vector
            x-axis limits.
        xlabel: str
            label of x-axis.
        ylabel: str
            label of y-axis.
        plot_title: str
            plot title.
        legend: boolean
            flag to switch on or off the plot legend. default is False.
        font_name: str
            font name.
        labels: list of strings
            labels 
        linewidth: float or vector/list of length 2
            line width of plotted lines
        fontsize: float
            font size
        formatter: boolean
            flag to switch on or off the internal plot formatter


    Returns:
        ax: matplotlib axis object or ndarray of axis objects

    '''
    if formatter:
        ax, font = default_plot_format(ax, font_name)
    else:
        font = FontProperties()
        if font_name is not None:
            font.set_name(font_name)

    if np.size(linewidth) == 1:
        linewidth = np.ones(2) * linewidth

    ax.plot(wvn, Xsec, label=labels[0], linewidth=linewidth[0])

    if len(XsecFit) > 0:
        ax.plot(wvn, XsecFit, '-.', label=labels[1], linewidth=linewidth[1])

    if xlim != None:
        ax.set_xlim(xlim[0], xlim[1])

    if xlabel != None:
        ax.set_xlabel(xlabel, fontproperties=font)

    if ylabel != None:
        ax.set_ylabel(ylabel, fontproperties=font)

    if plot_title != None:
        ax.set_title(plot_title, fontproperties=font)  # Add a title to the axes.
        ax.title.set_fontsize(fontsize)

    if legend:
        ax.legend()

    return ax


def scatter_plot(T, P, data, fig, ax, clim=None, xlabel='Temperature [K]',
                 ylabel='Pressure [hPa]', plot_title='', cbar_label='',
                 font_name=None, cmap='speed'):
    '''

    Args:
        T: vector
            temperatures.
        P: vector
            pressures same lenghth as `T`.
        data: vector
            data same length as `T`.
        fig: matplotlib figure object
            figure object.
        ax: matplotlib axis object or ndarray of axis objects
            matplotlib axis object or ndarray of axis objects.
        clim: vector or None
            value limits for the coloring.
        xlabel: str
            label of x-axis
        ylabel: str
            label of y-axis
        plot_title: str
            plot title
        cbar_label: str
            label of colorbar
        font_name: str
             font name
        cmap: str
            name of colormap

    Returns:
        fig: matplotlib figure object
            figure object.
        ax: matplotlib axis object or ndarray of axis objects
            matplotlib axis object or ndarray of axis objects.

    '''

    ax, font = default_plot_format(ax, font_name)

    if clim == None:
        clim = [None, None]

    MarkerSize = 50
    sca = ax.scatter(T, P / 100, MarkerSize, data, cmap=cmap, vmin=clim[0], vmax=clim[1])
    # ax.set_yscale('log')

    cbar = fig.colorbar(sca, ax=ax, shrink=1)
    cbar.set_label(cbar_label, fontproperties=font)

    ax.set_xlabel(xlabel, fontproperties=font)
    ax.set_ylabel(ylabel, fontproperties=font)
    ax.set_title(plot_title, fontproperties=font)

    return fig, ax


def pcolor_plot(x, y, Z, fig, ax, minZ, maxZ, font_name=None, xlabel=None, ylabel=None,
                cmap=None, title=None, cbar_label=None):
    '''
    wrapper to plot a 2d field

    Args:
        x: vector
            grid in x direction.
        y: vector
            grid in y direction.
        Z: matrix
            2d field, which will be plotted, with dimensions according to x and y
        fig: matplotlib figure object
            figure object.
        ax: matplotlib axis object or ndarray of axis objects
            matplotlib axis object or ndarray of axis objects.
        minZ: float
            minimum value of colorbar.
        maxZ: float
            minimum value of colorbar.
        font_name: str
             font name.
        xlabel: str
            label of x-axis.
        ylabel: str
            label of y-axis.
        cmap: str
            name of colormap.
        title: str
            plot title
        cbar_label: str
            label of colorbar.

    Returns:
        fig: matplotlib figure object
            figure object.
        ax: matplotlib axis object or ndarray of axis objects
            axis object or ndarray of axis objects.
        pcm: matplotlib pcolormesh object
            pcolormesh object.

        cbar: matplotlib colorbar object
            colorbar object.

    '''

    ax, font = default_plot_format(ax, font_name)

    if cmap == None:
        cmap = plt.get_cmap("Blues")

    # make plot and add colorbar
    pcm = ax.pcolormesh(x, y, Z, shading='nearest', cmap=cmap, vmin=minZ, vmax=maxZ)
    pcm.set_rasterized(True)
    cbar = fig.colorbar(pcm, ax=ax, shrink=1)
    # ax.set_yscale('log')

    # set the Make-Up and writings
    ax.set_title(title, fontproperties=font)
    ax.set_xlabel(xlabel, fontproperties=font)
    ax.set_ylabel(ylabel, fontproperties=font)

    cbar.set_label(cbar_label, fontproperties=font)

    return fig, ax, pcm, cbar


def make_band_patches(ax, bandwidths, verticalwidth, cmap=None, edgecolor='None',
                      alpha=0.25, zorder=-1):
    '''
    function to plot patches in plot to mark the band ranges
    Args:
        ax: matplotlib axis object or ndarray of axis objects
            axis object or ndarray of axis objects.
        bandwidths: list of two component vectors
            lower and upper border of defined bands

        verticalwidth: vector
            upper and lower vertical border for plotting
        cmap: Colormap
            Colormap(default: none)
        edgecolor: colormarker
            colormarker for the edges of the patches
        alpha: float
            alpha blending value, between 0 (transparent) and 1 (opaque).
        zorder: float
            Set the zorder for the artist. Artists with lower zorder values are drawn first.


    Returns:
        ax: matplotlib axis object or ndarray of axis objects
            axis object or ndarray of axis objects.

    '''

    if cmap is None:
        cmap = cmap_matlab_lines()

    for x, i in zip(bandwidths, range(len(bandwidths))):
        idx = i % np.size(cmap, axis=0)

        color = cmap[idx, :]

        patch = ptch.Rectangle((x[0], verticalwidth[0]), x[1] - x[0],
                               verticalwidth[1] - verticalwidth[0], facecolor=color,
                               alpha=alpha, edgecolor=edgecolor, zorder=zorder)

        ax.add_patch(patch)

    return ax


def plot_raw_data(xsec_data, species, font_name=None, max_num=10000):
    '''
    Function to plot overviews of the raw xsec data
    Args:
        xsec_data: list of xsec data
        species: str
            name of the species
        font_name: str
            font name
        max_num: int
            defines how many points of a spectrum is shown. If a
            spectrum has more more points than only a subset is shown.


    Returns:
        fig: matplotlib figure object
            figure object.
        ax: matplotlib axis object or ndarray of axis objects
            matplotlib axis object or ndarray of axis objects.

    '''

    number_of_sets = len(xsec_data)

    fig1, ax1 = default_figure(number_of_sets, 1, width_in_cm=20.9, height_in_cm=29.7)

    if number_of_sets == 1:
        ax1 = [ax1]

    for j in range(number_of_sets):

        ax1[j], font = default_plot_format(ax1[j], font_name)

        dw = []
        N_wvn = []
        for k in range(len(xsec_data[j])):

            wvn = np.linspace(xsec_data[j][k]['wmin'], xsec_data[j][k]['wmax'],
                              len(xsec_data[j][k]['xsec']))

            N_wvn.append(len(wvn))

            dw_k = (xsec_data[j][k]['wmax'] - xsec_data[j][k]['wmin']) / len(xsec_data[j][k]['xsec'])
            dw.append(dw_k)

            XSECS = xsec_data[j][k]['xsec']

            # if xsec are too detailed, make it it coarser.
            # it is just an overview plot, so not all details are needed.
            if len(XSECS) > max_num:
                idx = int(np.round(len(XSECS) / max_num))
                ax1[j].plot(wvn[0::idx], XSECS[0::idx], linewidth=0.1)
            else:
                ax1[j].plot(wvn, XSECS, linewidth=0.1)

        # Get highest resolution and number of samples
        dw = np.min(dw)
        N_wvn = np.max(N_wvn)

        ax1[j].set_yscale('log')
        ax1[j].set_ylim(1e-24, 1e-15)
        ax1[j].grid(which='both', linestyle=':', linewidth=0.25)
        ax1[j].set_ylabel('$a_{xsec} $[cm$^2$]')
        ax1[j].set_title(species + ': set ' + str(j) + '; $N_{obs}=$' + str(len(xsec_data[j])) +
                         '; $N_{sample,max}=$' + f'{N_wvn}' + '; $dwvn_{min}=$' + rf'{dw:.3f}' + 'cm$^{-1}$')
        ax1[j].title.set_fontsize(8)

        if j == number_of_sets:
            ax1[j].set_xlabel('wavenumber [cm$^{-1}$]')

    return fig1, ax1
