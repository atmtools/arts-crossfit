#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:53:05 2020

@author: Manfred Brath

This script calculates the absorption cross section for a desired species.
It is intended as an application example.

Usage:
Run the script

"""

import os

import numpy as np

import Xsec_aux_functions as xaf

# %% paths

# main plot folder
plot_folder = '../plots/Xsecs/'

# coefficients folder
coeff_folder = '../coefficients/'

# %% Calculate and plot cross sections

# Desired species, To get an overview (list) of the
# species use the function get_coeff_species within Xsec_aux_functions.
species = 'CFC11'

# Load cross section data
xsec_data = xaf.load_xsec_data(species, coeff_folder)

# Set temperature and pressure
temperature = 270  # [K]
pressure = 500e2  # [Pa]

# desired wavenumber [cm⁻¹]
wvn = np.linspace(750, 1250, 50001)

# Calculate cross sections
xsec = xaf.calculate_cross_sections(wvn, xsec_data, temperature=temperature, pressure=pressure)

# Plot cross sections
fig0, ax0 = xaf.default_figure(1, 1)
ax0 = xaf.plot_xsec(wvn, xsec, [], ax0,
                    xlabel='Wavenumber [cm$^{-1}$]',
                    ylabel='$a_{xsec}$ m$^2$',
                    plot_title=f"{pressure:.2f}$\,$Pa $-$ {temperature:.0f}$\,$K")

# save figure
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)

plotname0 = os.path.join(plot_folder, f'{species}.absoprtion_spectrum.{pressure:.2f}Pa.{temperature:.0f}K.pdf')
fig0.savefig(plotname0)
