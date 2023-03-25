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
import argparse
import glob

plt.rcParams['agg.path.chunksize'] = 10000

parser = argparse.ArgumentParser(description = "Inputs of CoSANG Analysis")
parser.add_argument('--dir', action='store', dest='dir', default='', help='Path to tag files')
##parser.add_argument('--gal_cent', action='store', dest='gal_cent', default=[27.5708,29.1913,27.5166], type=list, help='coordinates of galaxy center, defaults to m12i')

args = parser.parse_args()
dic_args = vars(args)


##gal_coord = args.gal_cent ##maybe figure out a process to get this
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
    ##new_x = []
    ##new_y = []
    ##new_z = []

    for i in range(len(met)):
        if x[i] != 0 and y[i] != 0 and stell_m[i] > 0 and met[i]>0.000000196: ##reasonable metalliciity and stell_m limits
            list_met.append(met[i])
            list_stell_m.append(stell_m[i])
            list_x.append(x[i] - 29.355671145527047) #correction factors, figure it out
            list_y.append(y[i] - 31.02505460762871)
            list_z.append(z[i] - 32.4927901257948)

##gal_x = (max(new_x) + min(new_x)) / 2
##gal_y = (max(new_y) + min(new_y)) / 2
##gal_z = (max(new_z) + min(new_z)) / 2
##print("Center coodinate:", gal_x, gal_y, gal_z)
##for j in range(len(list_met)):
  ##  list_x.append(new_x[j] - gal_x)
  ##  list_y.append(new_y[j] - gal_y)
  ##  list_z.append(new_z[j] - gal_z)

met = np.array(list_met)
stell_m = np.array(list_stell_m)
x = np.array(list_x)
y = np.array(list_y)
z = np.array(list_z)

avg_met = np.sum(met * stell_m) / np.sum(stell_m)
solar_frac_met = np.log10(met / 0.0196) ##makes metallicity in solar fraction 
rad = np.sqrt(((x) ** 2) + ((y) ** 2) + ((z) ** 2))
##print('x:', x)
##print('y:', y)
##print('z:', z)
rad_kpc = rad * 1000
x_kpc = x * 1000
y_kpc = y * 1000
z_kpc = z * 1000


##plot of mettalicity as a function of radius
plt.plot(rad_kpc, met,'.',markersize=1, c = 'blue')
m, b = np.polyfit(rad_kpc, met, 1)
y_vals = (m * rad_kpc) + b
plt.plot(rad_kpc, y_vals, c = 'black')
plt.axis('auto')
plt.title('Metallicity v. Radius')
plt.xlabel('Radius [kpc]')
plt.ylabel('Metallicity, [metal fraction]')
plt.savefig('rad_v_met.png')
plt.close()


##plot of metallicity weighted by solar mass
heatmap, xedges, yedges = np.histogram2d(rad_kpc, solar_frac_met, weights=stell_m ,bins=150)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
f = plt.figure()
f.set_figwidth(9)
f.set_figheight(9)
plt.clf()
with np.errstate(divide='ignore'):plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('viridis')
plt.colorbar(label='log10(Stellar Mass)')
plt.title('Metallicity v. Radius, Stellar Mass weighted')
plt.axis('auto')
plt.xlabel('Radius [kpc]')
plt.ylabel('Metallicity [log10(fraction of Sun)]')
plt.savefig('stell_m_weight_rad_v_met.png')

k = 1
bin_width = max(rad_kpc)/20
avg_met_values =[]
rad_line = []
met_err = []
while k * bin_width <= max(rad_kpc):
    met_line = []
    for j in range(len(solar_frac_met)):
        if rad_kpc[j] <= k * bin_width:
            met_line.append(solar_frac_met[j])
    avg_met_values.append(sum(met_line) / len(met_line))
    met_err.append(np.std(met_line))
    rad_line.append( (k-0.5) * bin_width)
    k+=1    

plt.plot(rad_line, avg_met_values, markersize=3, c='black')
plt.errorbar(rad_line, avg_met_values, yerr=met_err, c='black')
plt.savefig('stell_m_rad_met_w_avg_line.png')
plt.close()

high_met_x = []
high_met_y = []
high_met_z = []
high_stell_m = []

mid_met_x =[]
mid_met_y = []
mid_met_z = []
mid_stell_m = []

low_met_y = []
low_met_x = []
low_met_z =[]
low_stell_m = []

one_third_met = (max(met) - min(met)) / 3
two_third_met = 2 * one_third_met
print("making new plots")

for m in range(len(met)):
    if met[m] < one_third_met:
        low_met_x.append(x_kpc[m])
        low_met_y.append(y_kpc[m])
        low_met_z.append(z_kpc[m])
        low_stell_m.append(stell_m[m])
        #print("low met")
    elif two_third_met <= met[m]:
        high_met_x.append(x_kpc[m])
        high_met_y.append(y_kpc[m])
        high_met_z.append(z_kpc[m])
        high_stell_m.append(stell_m[m])
        #print("high met")
    else:
        mid_met_x.append(x_kpc[m])
        mid_met_y.append(y_kpc[m])
        mid_met_z.append(z_kpc[m])
        mid_stell_m.append(stell_m[m])
        #print("mid met")
    

heatmap, xedges, yedges = np.histogram2d(low_met_x, low_met_y, weights=low_stell_m ,bins=150)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
f = plt.figure()
f.set_figwidth(9)
f.set_figheight(9)
plt.clf()
with np.errstate(divide='ignore'):plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('viridis')
plt.colorbar(label='log10(Stellar Mass)')
plt.title('Low met, x vs. y, Stellar Mass weighted')
plt.axis('auto')
plt.xlabel('x [kpc]')
plt.ylabel('y [kpc]')
plt.savefig('x_y_low.png')
plt.close()


heatmap, xedges, yedges = np.histogram2d(low_met_x, low_met_z, weights=low_stell_m ,bins=150)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
f = plt.figure()
f.set_figwidth(9)
f.set_figheight(9)
plt.clf()
with np.errstate(divide='ignore'):plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('viridis')
plt.colorbar(label='log10(Stellar Mass)')
plt.title('Low met, x vs. z, Stellar Mass weighted')
plt.axis('auto')
plt.xlabel('x [kpc]')
plt.ylabel('z [kpc]')
plt.savefig('x_z_low.png')
plt.close()

heatmap, xedges, yedges = np.histogram2d(mid_met_x, mid_met_y, weights=mid_stell_m ,bins=150)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
f = plt.figure()
f.set_figwidth(9)
f.set_figheight(9)
plt.clf()
with np.errstate(divide='ignore'):plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('viridis')
plt.colorbar(label='log10(Stellar Mass)')
plt.title('Mid met, x vs. y, Stellar Mass weighted')
plt.axis('auto')
plt.xlabel('x [kpc]')
plt.ylabel('y [kpc]')
plt.savefig('x_y_mid.png')
plt.close()

heatmap, xedges, yedges = np.histogram2d(mid_met_x, mid_met_z, weights=mid_stell_m ,bins=150)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
f = plt.figure()
f.set_figwidth(9)
f.set_figheight(9)
plt.clf()
with np.errstate(divide='ignore'):plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('viridis')
plt.colorbar(label='log10(Stellar Mass)')
plt.title('Mid met, x vs. z, Stellar Mass weighted')
plt.axis('auto')
plt.xlabel('x [kpc]')
plt.ylabel('z [kpc]')
plt.savefig('x_z_mid.png')
plt.close()

heatmap, xedges, yedges = np.histogram2d(high_met_x, high_met_y, weights=high_stell_m ,bins=150)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
f = plt.figure()
f.set_figwidth(9)
f.set_figheight(9)
plt.clf()
with np.errstate(divide='ignore'):plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('viridis')
plt.colorbar(label='log10(Stellar Mass)')
plt.title('High met, x vs. y, Stellar Mass weighted')
plt.axis('auto')
plt.xlabel('x [kpc]')
plt.ylabel('y [kpc]')
plt.savefig('x_y_high.png')
plt.close()

heatmap, xedges, yedges = np.histogram2d(high_met_x, high_met_z, weights=high_stell_m ,bins=150)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
f = plt.figure()
f.set_figwidth(9)
f.set_figheight(9)
plt.clf()
with np.errstate(divide='ignore'):plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('viridis')
plt.colorbar(label='log10(Stellar Mass)')
plt.title('High met, x vs. z, Stellar Mass weighted')
plt.axis('auto')
plt.xlabel('x [kpc]')
plt.ylabel('z [kpc]')
plt.savefig('x_z_high.png')
plt.close()

print('low met:', min(met), 'one_third:', one_third_met, 'two_third:', two_third_met, 'high_met:', max(met))

