# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 16:18:28 2023

@author: alexi
"""

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
    #print(f_key) ##shows what keys we have to work with
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
solar_frac_met = []

for j in range(len(met)):
    solar_frac_met.append(np.log10( met[j] / 0.0196))##makes metallicity in solar fraction 
    
rad = np.sqrt(((x) ** 2) + ((y) ** 2) + ((z) ** 2))
##print('x:', x)
##print('y:', y)
##print('z:', z)
rad_kpc = rad * 1000
x_kpc = x * 1000
y_kpc = y * 1000
z_kpc = z * 1000


num = 20 #number of graphs you want
k = 1
#graph_range = [(((k-1) / num) * max(solar_frac_met)), ((k/num) * max(solar_frac_met))]##FIXME something went wrong 
for i in range(num):
    plot_x = []
    plot_y = []
    plot_z = []
    plot_stell_m = []
    graph_range = [((k / num) * max(solar_frac_met)), (((k+1)/num) * max(solar_frac_met))]
    for m in range(len(solar_frac_met)):
        #print("high:", graph_range[1], "low:", graph_range[0], "value:", solar_frac_met[m])
        if graph_range[1] < solar_frac_met[m] < graph_range[0]:
            plot_x.append(x_kpc[m])
            print('x:', x_kpc[m])
            plot_y.append(y_kpc[m])
            plot_z.append(z_kpc[m])
            plot_stell_m.append(stell_m[m])
            
    # if graph_range[0] > 0.0001:
    #     title_lower = np.log10(graph_range[0]/0.0196)
    # else:
    #     title_lower = -2.36967467411
    # title_upper = np.log10(graph_range[1]/0.0196)
    title_lower = graph_range[1]
    title_upper = graph_range[0]
    
    heatmap, xedges, yedges = np.histogram2d(plot_x, plot_y, weights=plot_stell_m ,bins=150)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    f = plt.figure()
    f.set_figwidth(9)
    f.set_figheight(9)
    plt.clf()
    with np.errstate(divide='ignore'):plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
    plt.set_cmap('viridis')
    plt.colorbar(label='log10(Stellar Mass)')
    plt.title( 'x vs. y, Stellar Mass weighted, met range:' + str(title_lower)+ ', ' + str(title_upper))
    plt.axis('auto')
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.savefig('x_y_' + str(k) + '.png')
    plt.close()
    
    heatmap, xedges, yedges = np.histogram2d(plot_x, plot_z, weights=plot_stell_m ,bins=150)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    f = plt.figure()
    f.set_figwidth(9)
    f.set_figheight(9)
    plt.clf()
    with np.errstate(divide='ignore'):plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
    plt.set_cmap('viridis')
    plt.colorbar(label='log10(Stellar Mass)')
    plt.title( 'x vs. z, Stellar Mass weighted, met range:' + str(title_lower)+ ', ' + str(title_upper))
    plt.axis('auto')
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    plt.savefig('x_z_' + str(k) + '.png')
    plt.close()
    
    k += 1


