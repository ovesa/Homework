#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 11:13:45 2018

Assignment: Write a computer function that returns a Sersic surface 
brightness profile given input parameters (effective radius, 
surface brightness at Re, and Sersic index) for the model and an array 
of input radii. The function should have an option to return the profile 
in either surface brightness or magnitude units. Use your function to make a 
plot of a set of various Sersic profiles that have the same effective radius 
and effective surface brightness, but different Sersic indices (0.5 < n < 6).

@author: oanavesa
"""
# import necessary imports
import numpy as np
import matplotlib.pyplot as plt
import math
from cycler import cycler


# defining my array of radii
input_radii = np.arange(0,25,.01)
# sersix index array (0.5 < n < 6) in increments of 0.5
sersic_index = np.arange(0.5,6.5,0.5)

''' Returns a Sersic surface brightness profile in units of flux with the 
       given paraments of surface brightness at the effect radius (brightre), 
       the effective radius (re), input radii (rad), and Sersic indices (n) '''
def surface_brightness(brightre,rad,re,n):
    bn = 2*n- 0.324
    exponent = np.exp(-bn*((rad/re)**(1/n)-1))
    Ir = brightre*exponent
    return Ir


''' Returns a Sersic surface brightness profile in units of magnitudes with 
      the given paraments of surface brightness at the effect radius (brightre), 
      the effective radius (re), input radii (rad), and Sersic indices (n) '''
def magnitude_units(mu_e,rad,re,n):
    bn = 2*n- 0.324
    mu = mu_e + 2.5*bn*math.log(np.exp(1),10)*((rad/re)**(1/n)-1)
    return mu


# takes the surface_brightess function and creates an array for each surface 
# brightness profile for each Sersic index
surface = []
for i in sersic_index:
    t = surface_brightness(25,input_radii,8,i)
    surface.append(t)

# using this for plotting the sersic indices in the legend of the plot
lmk =np.where(sersic_index)[0]
p = sersic_index[lmk]

# plots the surface brightness profile in units of flux
plt.figure()
plt.subplots(figsize=(5,5))
for j in surface:
    plt.rc('axes', prop_cycle=(cycler('color', ['blue', 'green', 'red', 'darkred', 'magenta', 'brown', 'black', 'purple', 'pink', 'orange', 'teal', 'coral', 'lightblue', 'lime', 'lavender', 'turquoise', 'darkgreen', 'tan', 'salmon', 'gold', 'darkblue'])))  
    line = plt.plot(input_radii,j)
    plt.yscale('log')
    #plt.xscale('log')
    plt.title('Radial Surface Brightness Profile',fontweight='bold')
    plt.xlabel('Radius (kpc)',fontweight='bold')
    plt.ylabel('Surface Brightness',fontweight='bold')
    for k in range(len(p)):
        plt.legend((p),title="Sersic Index (n)",fancybox = True,loc='upper right',
          ncol=2)
plt.savefig('Vesa_Sersic_Brightness.png')
plt.show()

# takes the magnitude_units function and creates an array for each surface 
# brightness profile for each Sersic index
magnitudes = []
for d in sersic_index:
    w = magnitude_units(25,input_radii,8,d)
    magnitudes.append(w)


# plots the surface brightness profile in units of magnitude
plt.figure()
plt.subplots(figsize=(5,5))
for q in magnitudes:
    plt.rc('axes', prop_cycle=(cycler('color', ['blue', 'green', 'red', 'darkred', 'magenta', 'brown', 'black', 'purple', 'pink', 'orange', 'teal', 'coral', 'lightblue', 'lime', 'lavender', 'turquoise', 'darkgreen', 'tan', 'salmon', 'gold', 'darkblue'])))  
    plt.plot(input_radii,q)
    #plt.gca().invert_yaxis()
    #plt.xscale('log')
    #plt.yscale('log')
    plt.title('Sersic Brightness Profile in Magnitude Units', fontweight='bold')
    plt.xlabel('Radius (kpc)',fontweight='bold')
    plt.ylabel('Magnitude Brightness (mag*kpc^-2)',fontweight='bold')
    for k in range(len(p)):
        plt.legend((p),title="Sersic Index (n)",fancybox = True,loc='upper right',
          ncol=2)
ax = plt.gca()
ax.invert_yaxis()
plt.savefig('Vesa_Sersic_Brightness_Mag.png')
plt.show()




