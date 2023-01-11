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

gal_coord = [27.5708,29.1913,27.5166] ##hardcoded, make into argument, maybe figure out a process to get this
list_met = []
list_x = []
list_y = []
list_z = []
list_stell_m = []


for j in args.path: #FIXME, make arguments for path to snapshot files
    f = h5py.File(str(j), 'r') #opens the file
    f_key = list(f.keys())
    print(f_key) ##shows what keys we have to work with
    
    ##gets the data and puts it in the right format
    def make_for_calc(ds, neginf = False, center_val = []):
        if neginf == True:
            new = np.nan_to_num(np.array(ds), neginf = 0.0)
        else:
            new = np.nan_to_num(np.array(ds))
    return new
    met = make_for_calc(f['Metallicity'])
    stell_m = make_for_calc(f['StellarMass'], neginf = True)
    x = make_for_calc(f['X'])
    y = make_for_calc(f['Y'])
    z = make_for_calc(f['Z'])

    new_met = []
    new_stell_m = []
    new_x = []
    new_y = []
    new_z = []
    for i in range(len(met)):
        if x[i] != 0 and y[i] != 0:
            new_met.append(met[i])
            new_stell_m.append(stell_m[i])
            new_x.append(x[i] - gal_coord[0] - 1.8383) #correction factors, figure it out
            new_y.append(y[i] - gal_coord[1] - 1.7024)
            new_z.append(z[i] - gal_coord[2] - 4.2)
    
    list_x.append(new_x)
    list_y.append(new_y)
    list_z.append(new_z)
    list_met.append(new_met)
    list_stell_m.append(new_stell_m)

met = np.array(list_met)
stell_m = np.array(list_stell_m)
x = np.array(list_x)
y = np.array(list_y)
z = np.array(list_z)

avg_met = np.sum(met * stell_m) / np.sum(stell_m)
rad = np.sqrt(((x) ** 2) + ((y) ** 2) + ((z) ** 2))
rad_kpc = rad * 1000

##plot of mettalicity as a function of radius
plt.plot(rad_kpc,met,'.',markersize=1, c = 'blue')
m, b = np.polyfit(rad_kpc, met, 1)
y_vals = (m * rad_kpc) + b
plt.plot(rad_kpc, y_vals, c = 'black', linewidth = 100)
plt.axis('auto')
plt.title('Metallicity v. Radius')
plt.xlabel('Radius [kpc]')
plt.ylabel('Metallicity [Fe/H]')
plt.show()
plt.close()

##plot of metallicity weighted by solar mass
heatmap, xedges, yedges = np.histogram2d(rad_kpc, met, weights=stell_m_for_calc,bins=50)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.clf()
plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('autumn')
plt.colorbar(label='log10(Stellar Mass)')
plt.title('Metallicity v. Radius, Stellar Mass weighted')
plt.axis('auto')
plt.xlabel('Radius [kpc]')
plt.ylabel('Metallicity [Fe/H]')
plt.show()
plt.close()

f.close()
