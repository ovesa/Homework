#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 23 13:56:58 2018

ASTR 555: Homework #3

study the nearby galaxy population and construct a rough galaxy luminosity 
function for nearby galaxies using SDSS.

@author: oanavesa

"""
#%%

# importing the necessary imports
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import scipy as sp
from scipy import stats
import pandas as pd
#%%

# defining constants
c = 300000 #(km/s) # speed of light
H0 = 70 #(km s^-1 Mpc^-1) =  Hubble Constant

''' Problem 1 part b) part i '''

def distance(c,z,H0):
    '''Calculating the distance to the galaxies.
    Input: c = speed of light (km/s); z = redshift; H0 = Hubble Constant  (km s^-1 Mpc^-1)
    Output: d = distance (pc)'''
    d = (c*z)/H0
    d=d*1000000
    return d

''' Problem 1 part b) part ii '''
def absolutemagnitude(m,d):
    '''Calculating the absolute magnitudes of the galaxies.
    Input: m = apparent magnitude; d = distance (pc)
    Output: M = absolute magnitude'''
    M = m - 5*np.log10(d/10)
    return M

# defining the path
path = '/Users/oanavesa/Desktop/GradSchool/FirstYear/Galaxies1/Homework/HW3/ew2.xlsx'


#%%

# importing the columns
df = pd.read_excel(path)
ra  = df['ra']
dec  = df['dec']
petroMag_g = df['petroMag_g']
petroMag_u = df['petroMag_u']
petroMag_i = df['petroMag_i']
petroMag_r = df['petroMag_r']
petroMag_z = df['petroMag_z']
corrected_g = df['corrected_g']
corrected_u = df['corrected_u']
corrected_z = df['corrected_z']
corrected_i = df['corrected_i']
corrected_r = df['corrected_r']
redshift = df['redshift']

# empty lists for the absolute magnitudes and distance
absolute_g = []
absolute_u = []
distances = []

# calculating the distances and appending them
for i in df.index:
    dis = distance(c,redshift[i],H0)
    distances.append(dis)
    
# calculating the absolute magnitudes in the G-band and appending them  
for j in df.index:
    mag  = absolutemagnitude(corrected_g[j],distances[j])
    absolute_g.append(mag)

# calculating the absolute magnitudes in the U-band and appending them  
for k in df.index:
    mag  = absolutemagnitude(corrected_u[k],distances[k])
    absolute_u.append(mag)
    
# Calculating the difference between the U-G bands and appending them
diff = []
for t in range(0,len(absolute_u)):
    diff.append(absolute_u[t] - absolute_g[t])

#%%

''' Problem 1 part b) part iii '''

# plotting the galaxy color vs absolute magnitude
plt.scatter(diff,absolute_g)
plt.ylim([-25,-1])
plt.title('Galaxy Color vs. Absolute Magnitude', weight='bold',fontsize=14)
plt.xlabel('U-G', weight='bold',fontsize=12)
plt.ylabel('Absolute Magnitude of G-Filter', weight='bold',fontsize=12)
plt.savefig('galaxycolour.png')
plt.show()

# plotting a histogram of the raw luminosity function
plt.hist(absolute_u, range=(-22,-10))
plt.title('Raw Luminosity Function', weight='bold',fontsize=14)
plt.ylabel('Number in Each Bin', weight='bold',fontsize=12)
plt.xlabel('Absolute Magnitude ', weight='bold',fontsize=12)
plt.savefig('rawlumfunction.png')
plt.show()


    
#%%
''' Problem 1 c) part ii: calculatimg the effective volume of the sample as a 
    function of absolute magnitude'''


# apparent magnitude cut-off
sdss_m = 17.77

# calculating the maximum distance for the maximum volume calculation
distancemax =[]
for i in range(0,len(absolute_g)):
   d = (10**((sdss_m - absolute_g[i])/5))*10
   d = d/1000000
   distancemax.append(d)
   
# calculating the minimum distance 
dmin = distance(c,0.015,70)
dmin_mpc = dmin/1000000

# calculating the volume to the most distant object
v_max=[]
for t in range(0,len(distancemax)):
    if distancemax[t] < dmin_mpc:
        pass
    else:
        v = ((4*np.pi)/3)*((distancemax[t])**3-(dmin_mpc)**3)
    v_max.append(v)


# doing 1 / the volume
one_divide_vmax = []
for r in range(0,len(v_max)):
    divide = 1/v_max[r]
    one_divide_vmax.append(divide)
    

#%%

''' Problem 1 c) part iii: create a luminosity function'''

# creates the luminosity function
lum, mag, v3 = sp.stats.binned_statistic(absolute_g, one_divide_vmax, statistic='sum', bins=100, range=(-22,-10) )


# plotting the luminosity function
plt.semilogy(mag[:100],lum)
plt.xlim([-16.6,-23])
plt.title('Luminosity Selection Function', fontsize=14, weight='bold')
plt.ylabel('Luminosity Function', weight='bold',fontsize=12)
plt.xlabel('Absolute Magnitude', weight='bold',fontsize=12)
#plt.ylim(top=10**-3)
plt.savefig('luminosityfunction.png')
plt.show()



#%%

''' Problem 1 c) part iv: calculate the ratio of the volume within the distance
    to the object to the volume within the maximum distance '''

# calculating the volume within the distance to the object
actual_volume = []
for q in range(0,len(distances)):
    di = distances[q]/1000000
    actual = (9200/41250)*((4*np.pi)/3)*di**3
    actual_volume.append(actual)
    
# calculating the ratio
volume_ratio = []
for t in range(0,len(actual_volume)):
    ratio = actual_volume[t]/v_max[t]
    volume_ratio.append(ratio)
    

plt.hist(volume_ratio,bins=100, range=(0,1))
plt.title('V/V_max', fontsize=14, weight='bold')
plt.xlabel('V/V_max', weight='bold',fontsize=12)
plt.ylabel('Number in Each Bin', weight='bold',fontsize=12)
#plt.savefig('volumeratio.png')
plt.show()
