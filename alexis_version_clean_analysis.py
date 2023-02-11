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
import sys
import os
import argparse
import glob

parser = argparse.ArgumentParser(description = "Inputs of CoSANG Analysis")
parser.add_argument('--dir', action='store', dest='dir', default='', help='Path to tag files')
parser.add_argument('--gal_cent', action='store', dest='gal_cent', default=[27.5708,29.1913,27.5166], type=list, help='coordinates of galaxy center, defaults to m12i')

args = parser.parse_args()
dic_args = vars(args)


gal_coord = args.gal_cent ##maybe figure out a process to get this
list_met = []
list_x = []
list_y = []
list_z = []
list_stell_m = []


def make_for_calc(ds, neginf = False, center_val = []):
    if neginf == True:
        new = np.nan_to_num(np.array(ds), neginf = 0.0)
    else:
        new = np.nan_to_num(np.array(ds))
    return new

adress = args.dir+'/*.h5'

for h5name in glob.glob(adress): #reads in all of the files in the directory
    f = h5py.File(str(h5name), 'r') #opens the file
    f_key = list(f.keys())
    print(f_key) ##shows what keys we have to work with
    met = make_for_calc(f['Metallicity'])
    stell_m = make_for_calc(f['StellarMass'], neginf = True)
    x = make_for_calc(f['X'])
    y = make_for_calc(f['Y'])
    z = make_for_calc(f['Z'])

    for i in range(len(met)):
        if x[i] != 0 and y[i] != 0 and stell_m[i] > 0 and met[i]!=0:
            list_met.append(met[i])
            list_stell_m.append(stell_m[i])
            list_x.append(x[i] - gal_coord[0] - 1.8383) #correction factors, figure it out
            list_y.append(y[i] - gal_coord[1] - 1.7024)
            list_z.append(z[i] - gal_coord[2] - 4.2)

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
plt.savefig('rad_v_met.png')
plt.close()

##plot of metallicity weighted by solar mass
heatmap, xedges, yedges = np.histogram2d(rad_kpc, met, weights=stell_m ,bins=50)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.clf()
with np.errstate(divide='ignore'):plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('autumn')
plt.colorbar(label='log10(Stellar Mass)')
plt.title('Metallicity v. Radius, Stellar Mass weighted')
plt.axis('auto')
plt.xlabel('Radius [kpc]')
plt.ylabel('Metallicity [metal fraction]')
plt.savefig('stell_m_weight_rad_v_met.png')
plt.close()
