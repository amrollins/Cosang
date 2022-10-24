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
import os
import math
from numpy import genfromtxt

###
#IMPORT/CLEAN TAGS
###
#each tag is a csv without the main halo
#write this section to read in hdf5 tags

rawID = []
rawMet = []
rawStelm = []

for i in range(153,264):
    if i==231 or i==237 or i==238 or i==239:            #rewrite this to deal directly with the hdf5 files, maybe use the output of TagMix from Sharaham?
        continue
    file = 'tags'+str(i)+'.csv'
    tag = genfromtxt(file, delimiter=',')
    if len(tag) != 0:
        rawID = np.concatenate((rawID,tag[:,0]))
        rawMet = np.concatenate((rawMet,tag[:,1]))
        rawStelm = np.concatenate((rawStelm,tag[:,2]))
        
   
#remove nan & zero

ID = []
Met = []
Stelm = []

for i in range(len(rawID)):
    if math.isnan(rawStelm[i])==False and math.isnan(rawMet[i])==False and rawStelm[i] != math.inf:
        ID.append(rawID[i])
        Met.append(rawMet[i])
        Stelm.append(rawStelm[i])


###
#IMPORT POS DATA FROM LAST SNAPSHOT
###
#index of positions is the same as their ID, so in this format IDlast isn't needed

pos264 = genfromtxt('264pos.csv',delimiter=',') ##write this to get the data needed
#IDlast = pos264[1:,1]
Xlast = pos264[1:,2]
Ylast = pos264[1:,3]
Zlast = pos264[1:,4]

###
#GET DATA WITHIN VIRIAL RADIUS
###

maingal = [51.0098,53.4828,47.2304]                     #Pos of main galaxy, came from 0.0.h5 of 264
Mv = 8.49862*10**11                                     #virial mass
Rv = 0.154039                                         #virial radius
#radius=0.2                                              #Radius around center, 200 kpc
                                                        #actually creates box sides=0.4 around center

rMet = []
rStelm = []
rX = []
rY = []
rZ = []
rR = []

for i in range(len(ID)):
    r = math.sqrt((Xlast[int(ID[i])]-maingal[0])**2 +(Ylast[int(ID[i])]-maingal[1])**2 +(Zlast[int(ID[i])]-maingal[2])**2)
    if r <= Rv:
        rMet.append(Met[i])
        rStelm.append(Stelm[i])
        rX.append(Xlast[int(ID[i])])
        rY.append(Ylast[int(ID[i])])
        rZ.append(Zlast[int(ID[i])])
        rR.append(r)


# =============================================================================
# for i in range(len(ID)):
#     if Xlast[int(ID[i])]>= (maingal[0]-radius) and Xlast[int(ID[i])]<=(maingal[0]+radius):
#         if Ylast[int(ID[i])]>= (maingal[1]-radius) and Ylast[int(ID[i])]<=(maingal[1]+radius):
#             if Zlast[int(ID[i])]>= (maingal[2]-radius) and Zlast[int(ID[i])]<=(maingal[2]+radius):
#                 rMet.append(Met[i])
#                 rStelm.append(Stelm[i])
#                 rX.append(Xlast[int(ID[i])])
#                 rY.append(Ylast[int(ID[i])])
#                 rZ.append(Zlast[int(ID[i])])
# =============================================================================

#%%  
###
#ANALYSIS VALUES
###

#average metallicity, stellar mass weighted
avg_met = np.sum(np.array(rMet)*np.array(rStelm))/np.sum(np.array(rStelm))
fixed_avg_met = np.log10(avg_met*(10**10)/(0.015))

print('Average metallicity (simulation value): '+str(avg_met))
print('Average metallicity (physical value): '+str(fixed_avg_met)+'\n')

#correct metallicity
fixedMet = np.log10(np.array(rMet)*(10**10)/(0.015))

print('Maximum metallicity (simulation value): '+str(max(rMet)))
print('Maximum metallicity (physical value): '+str(max(fixedMet))+'\n')

print('Minimum metallicity (simulation value): '+str(min(rMet)))
print('Minimum metallicity (physical value): '+str(min(fixedMet))+'\n')

#avgerage metallicity at 30kpc
rR = np.array(rR)
rMet = np.array(rMet)
rStelm = np.array(rStelm)
Met30 = rMet[np.where(np.logical_and(rR>=0.03, rR<=0.04))]      #set correct range
Stelm30 = rStelm[np.where(np.logical_and(rR>=0.03, rR<=0.04))]
avg_met30 = np.sum(np.array(Met30)*np.array(Stelm30))/np.sum(np.array(Stelm30))
fixed_avg_met30 = np.log10(avg_met30*(10**10)/(0.015))
fixedMet30 = np.log10(np.array(Met30)*(10**10)/(0.015))

print('Average metallicity at 30 kpc (simulation value): '+str(avg_met30))
print('Average metallicity at 30 kpc (physical value): '+str(fixed_avg_met30)+'\n')

print('Maximum metallicity at 30 kpc (simulation value): '+str(max(Met30)))
print('Maximum metallicity at 30 kpc (physical value): '+str(max(fixedMet30))+'\n')

print('Minimum metallicity at 30 kpc (simulation value): '+str(min(Met30)))
print('Minimum metallicity at 30 kpc (physical value): '+str(min(fixedMet30))+'\n')



#%%
###
#PLOT R v. MET
###
#how to weight with stellar mass?

Rkpc = rR*1000                              #axes into kpc

#points
plt.figure(0)
plt.plot(Rkpc,fixedMet,'.',markersize=1)
m, b = np.polyfit(Rkpc, fixedMet, 1)
plt.plot(Rkpc, m*Rkpc + b)
plt.title('Metallicity v. Radius')
plt.xlabel('Radius [kpc]')
plt.ylabel('Metallicity [Fe/H]')
plt.axis('auto')

print('Gradient/Slope NOT weighted: '+str(m)+' [Fe/H]/kpc')


#stellar mass weighted histogram
plt.figure(1)
heatmap, xedges, yedges = np.histogram2d(Rkpc,fixedMet,weights=rStelm,bins=50)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.clf()
plt.imshow(np.log10(heatmap.T), extent=extent, origin='lower')
plt.set_cmap('gist_heat_r')
plt.colorbar(label='log10(Stellar Mass)')
plt.axis('auto')
plt.title('Metallicity v. Radius, Stellar Mass weighted')
plt.xlabel('Radius [kpc]')
plt.ylabel('Metallicity [Fe/H]')
plt.show()

#%%
#average bins into scatter plot

from scipy import stats
#average metallicity in bins
heatmap, xedges, yedges = np.histogram2d(Rkpc,rMet,weights=rMet,bins=50)

ret = stats.binned_statistic_2d(Rkpc, rMet,rMet, statistic='mean', bins=50)

x_med = []
for i in range(len(ret.x_edge)-1):
    max = ret.x_edge[i+1]
    min = ret.x_edge[i]
    med = (min+max)/2
    x_med.append(med)
#%%    
print(ret.statistic)
ret_stat = np.nan_to_num(ret.statistic)
print(ret_stat)
#%%

#make in single array
ret_stat_single = []
for i in ret_stat:
    for x in i:
        ret_stat_single.append(x)
ret_stat_single = np.array(ret_stat_single)
    
fixed_medMet = np.log10(ret_stat_single*(10**10)/(0.015))
#print(fixed_medMet)

print('.......')
print(fixed_medMet)
print(x_med*50)
#x_med_repeat = np.repeat(x_med,50)
#print(x_med_repeat)
plt.plot(x_med*50,fixed_medMet,'.')

#%%
m, b = np.polyfit(x_med*50, fixed_medMet, 1)
m = np.nan_to_num(m,neginf=0.0)
b = np.nan_to_num(b,neginf=0.0)
print(m)
plt.plot(x_med*50, m*fixed_medMet + b)
