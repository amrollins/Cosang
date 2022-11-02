#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 11:00:27 2021

@author: annelia
"""
#Analysis on sample halo

###
#NOTES
###
#Hardcoded:
#csv inputs
#center of main galaxy position
#virial radius
#virial mass
#
#User inputs:
#radius of observation?
###

import h5py
import numpy as np
import matplotlib.pyplot as plt
import math 

f = h5py.File('test_tag.h5', 'r') #opens the file
f_key = list(f.keys())
print(f_key) ##shows what keys we have to work with

##gets the data and calculates the values we want
met = f['Metallicity']
met_for_calc = np.nan_to_num(np.array(met))
fixed_met = np.log10(met_for_calc * (10 ** 10) /  (0.015))
stell_m = f['StellarMass']
stell_m_for_calc = np.nan_to_num(np.array(stell_m), neginf = 0.0)
avg_met = np.sum(met_for_calc * stell_m_for_calc) / np.sum(stell_m_for_calc)
x = f['X']
x_for_calc = np.nan_to_num(np.array(x))
y =  f['Y']
y_for_calc = np.nan_to_num(np.array(y))
z = f['Z']
z_for_calc = np.nan_to_num(np.array(z))
rad = np.sqrt(((x_for_calc) ** 2) + ((y_for_calc) ** 2) + ((z_for_calc) ** 2))
rad_kpc = rad * 1000

##plot of mettalicity as a function of radius
plt.plot(rad_kpc,fixed_met,'.',markersize=1, c = 'blue')
m, b = np.polyfit(rad_kpc, fixed_met, 1)
y_vals = (m * rad_kpc) + b
plt.plot(rad_kpc, y_vals, c = 'black', linewidth = 100)
plt.title('Metallicity v. Radius')
plt.xlabel('Radius [kpc]')
plt.ylabel('Metallicity [Fe/H]')
plt.axis('auto')
plt.show()
plt.close()

##plot of metallicity weighted by solar mass
heatmap, xedges, yedges = np.histogram2d(rad_kpc, np.nan_to_num(fixed_met, neginf = 0.0), weights=stell_m_for_calc,bins=50)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.clf()
plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('autumn')
plt.colorbar(label='log10(Stellar Mass)')
plt.title('Metallicity v. Radius, Stellar Mass weighted')
#plt.axis('auto')
plt.axes(xlim=(45000, 60000), ylim=(0, 10), autoscale_on=False)
plt.xlabel('Radius [kpc]')
plt.ylabel('Metallicity [Fe/H]')
plt.show()
