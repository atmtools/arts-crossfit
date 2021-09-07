#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 11:41:02 2021

@author: u242031
"""


import numpy as np
import Xsec_aux_functions as xaf


number_of_sets = 6

wvn_intervalls = np.zeros((number_of_sets, 2))
dws = np.zeros(number_of_sets)
Nfs = np.zeros(number_of_sets)

wvn_intervalls[0, :] = [500.0157000, 6500.0125000]
wvn_intervalls[1, :] = [680.011,750.0080000]
wvn_intervalls[2, :] = [1060.9942000,1165.0002000]
wvn_intervalls[3, :] = [1075.0000000,1165.0100000]
wvn_intervalls[4, :] = [1170.0100000,1380.0000000]
wvn_intervalls[5, :] = [1219.9963000,1285.0007000]


Nfs[0]=99570.
Nfs[1]=4650.
Nfs[2]=34516.
Nfs[3]=2990.
Nfs[4]=6970.
Nfs[5]=21573.

dws=(wvn_intervalls[:,1]-wvn_intervalls[:,0])/(Nfs-1)


#for testing set dw for interval 2 to value of interval 0
dws[1]=dws[0]

def suggest_banddefinition(wvn_intervalls, dws):

    sorted_limits=np.sort(wvn_intervalls.flatten())

    band_limits=[]
    band_dw=[]

    for i in range(np.size(sorted_limits)-1):

        band_limit_i=[sorted_limits[i],sorted_limits[i+1]]

        ol_i=np.zeros(number_of_sets)

        for j in range(number_of_sets):

            ol_i[j]=xaf.getOverlap(wvn_intervalls[j, :] ,band_limit_i)

        if np.sum(ol_i)>0:
            #get minimum dw per band (highest resolution)
            dw_i = np.min([dws[j] for j in range(number_of_sets) if ol_i[j]>0])

            #Check for too small bands. A band needs at least 2 points.
            #This means dw_i must <= band_end-band_start
            if dw_i <= band_limit_i[1] - band_limit_i[0]:

                band_limits.append(band_limit_i)
                band_dw.append(dw_i)


    #now check if two adjacent bands have the same wavenumber resolution
    band_limits_final=[]

    marker=np.zeros(len(band_limits))

    cnt=0
    for i in range(len(band_limits)-1):

        adjacent = (band_limits[i][1]-band_limits[i+1][0] == 0.)
        same_dw = (band_dw[i] == band_dw[i+1])

        crit= adjacent and same_dw

        if not crit:
            cnt+=1

        marker[i+1]=cnt

    for i in range(cnt+1):

        temp=[band_limits[j] for j in range(len(band_limits)) if marker[j] == i ]

        band_limits_final.append([np.min(temp),np.max(temp)])

    return band_limits_final



band_limits_final = suggest_banddefinition(wvn_intervalls, dws)

print(band_limits_final)

