#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 22:04:59 2018

ASTR 555 Homework 4.

Problem 1: Normalize the Salpeter and Kroupa et al. (1993) IMFs above 1 solar
mass. This is done in the latex document.

Problem 2:  Calculate the total integrated luminosity and to the total 
integrated mass as a function of stellar mass for each age using a chosen IMF.

Problem 3: Create a predicted spectrum for these stellar populations using
the Planck Function and a chosen IMF.

I elected to use the Salpeter IMF. The isochrones are from the 
YaPSI (YALE-POTSDAM) group for a solar metallicity population at ages of 
roughly 100 Myr, 1 Gyr, and 5 Gyr.

@author: oanavesa
"""

# import necessary modules
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

# path files
# used isochrones from YAPSI
path = '/Users/oanavesa/Desktop/GradSchool/FirstYear/Galaxies1/Homework/'
gyr1 = np.loadtxt(path+'gyr1.txt', unpack=True)
gyr5 = np.loadtxt(path+'gyr5.txt', unpack=True)
myr100 = np.loadtxt(path+'myr1.txt', unpack=True)

# masses
Mass_gyr1 = gyr1[1] 
Mass_gyr5 = gyr5[1] 
Mass_myr100 = myr100[1]

# temperatures
Temperature_gyr1 = gyr1[2] 
Temperature_gyr5 = gyr5[2] 
Temperature_myr100 = myr100[2]

#constants
N0 = 0.03726060064088233
c = 3.00*10**(8) # m/s
h = 6.63*10**(-34) # J*s
k = 1.38*10**(-23)  # J/K


'''Problem 2: Calculate the total integrated luminosity and to the total 
integrated mass as a function of stellar mass for each age. '''


# IMF = Salpeter IMF
# get normalized IMF = N0 from Problem 1 (k in Problem 1)
def IMF(M):
    '''Salpeter IMF (initial mass function);
    Input: N0 = constant that was normalized in Problem 1,
    M = mass in solar masses;
    Output: return the IMF '''
    imf = N0*M**(-2.35)
    return imf




# total mass integral: IMF*M
def mass(M):
    '''Total mass integral;
    Input: N0 = constant that was normalized in Problem 1,
    M = mass;
    Out: returns the total mass '''
    integral =  N0*(M**(-2.35))*M
    return integral

# total luminosity integral: IMF *M^{3.5}
# assuming a main-sequence stellar population
def luminosity(M):
    ''' Total luminosity integral;l
    Input:  N0 = constant that was normalized in Problem 1,
    M = mass;
    Output: returns the total luminosity
    I assumed that the stars in these stellar populations are main-sequence stars'''
    integral = N0*M**(-2.35)*M**(3.5)
    return integral
    
# finds the total integrated mass and luminosity

# Gyr =1 
mass_list_gyr1 = []
luminosity_list_gyr1 = []
for i in range(0, len(Mass_gyr1)-1):
    # upper mass limit
    mass_bin_1 = Mass_gyr1[i+1] 
    # lower mass limit
    mass_bin_2 = Mass_gyr1[i]
    # integral for mass
    integrate = quad(mass,mass_bin_2,mass_bin_1)[0]
    # integral for luminosity
    integrate2 = quad(luminosity,mass_bin_2,mass_bin_1)[0]
    mass_list_gyr1.append(integrate)
    luminosity_list_gyr1.append(integrate2)
    
    
#Gyr = 5
mass_list_gyr5 = []
luminosity_list_gyr5 = []
for i in range(0, len(Mass_gyr5)-1):
    # upper mass limit
    mass_bin_1 = Mass_gyr5[i+1] 
    # lower mass limit
    mass_bin_2 = Mass_gyr5[i]
    # integral for mass
    integrate = quad(mass,mass_bin_2,mass_bin_1)[0]
    # integral for luminosity
    integrate2 = quad(luminosity,mass_bin_2,mass_bin_1)[0]
    mass_list_gyr5.append(integrate)
    luminosity_list_gyr5.append(integrate2)
    
    
#Myr = 100
mass_list_myr100 = []
luminosity_list_myr100 = []
for i in range(0, len(Mass_myr100)-1):
    # upper mass limit
    mass_bin_1 = Mass_myr100[i+1] 
    # lower mass limit
    mass_bin_2 = Mass_myr100[i]
    # integral for mass
    integrate = quad(mass,mass_bin_2,mass_bin_1)[0]
    # integral for luminosity
    integrate2 = quad(luminosity,mass_bin_2,mass_bin_1)[0]
    mass_list_myr100.append(integrate)
    luminosity_list_myr100.append(integrate2)


# to find cumulative sums of the masses for each stellar population age
# https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.cumsum.html
    
# Gyr = 1
cumulative_distribution_mass_gyr1 = np.cumsum(mass_list_gyr1)
cumulative_distribution_luminosity_gyr1 = np.cumsum(luminosity_list_gyr1)
# Gyr = 5
cumulative_distribution_mass_gyr5 = np.cumsum(mass_list_gyr5)
cumulative_distribution_luminosity_gyr5 = np.cumsum(luminosity_list_gyr5)
# Myr = 100
cumulative_distribution_mass_myr100 = np.cumsum(mass_list_myr100)
cumulative_distribution_luminosity_myr100 = np.cumsum(luminosity_list_myr100)


# need to normalize the cumulative sums to 1
# used https://stackoverflow.com/questions/29661574/normalize-numpy-array-columns-in-python

# Gyr = 1
mass_gyr1_normalized = (cumulative_distribution_mass_gyr1 - cumulative_distribution_mass_gyr1.min()) / cumulative_distribution_mass_gyr1.ptp()
luminosity_gyr1_normalized = (cumulative_distribution_luminosity_gyr1 - cumulative_distribution_luminosity_gyr1.min()) / cumulative_distribution_luminosity_gyr1.ptp()
# Gyr = 5
mass_gyr5_normalized = (cumulative_distribution_mass_gyr5 - cumulative_distribution_mass_gyr5.min()) / cumulative_distribution_mass_gyr5.ptp()
luminosity_gyr5_normalized = (cumulative_distribution_luminosity_gyr5 - cumulative_distribution_luminosity_gyr5.min()) / cumulative_distribution_luminosity_gyr5.ptp()
# Myr = 1
mass_myr100_normalized = (cumulative_distribution_mass_myr100 - cumulative_distribution_mass_myr100.min()) / cumulative_distribution_mass_myr100.ptp()
luminosity_myr100_normalized =  (cumulative_distribution_luminosity_myr100 - cumulative_distribution_luminosity_myr100.min()) / cumulative_distribution_luminosity_myr100.ptp()

# Plots the Total Mass and Luminosity Contributions

plt.figure()
f, (ax1, ax2,ax3) = plt.subplots(1, 3, sharey=True)
f.suptitle('Total Mass and Luminosity Contributions',fontsize=14,fontweight='bold',y=.99)

ax1.plot(Mass_gyr1[1:], mass_gyr1_normalized, color= 'k', label = 'Mass Contribution',linestyle='solid')
ax1.plot(Mass_gyr1[1:], luminosity_gyr1_normalized, color = 'r', label =' Luminosity Contribution',linestyle='solid')
ax1.set_title('Gyr = 1')
ax1.set_xlabel('Mass (Msun)',fontweight='bold')
ax1.set_ylabel('Cumulative Contribution',fontweight='bold')

ax2.plot(Mass_gyr5[1:], mass_gyr5_normalized, color = 'k', label = 'Mass Contribution',linestyle='solid')
ax2.plot(Mass_gyr5[1:], luminosity_gyr5_normalized, color = 'r', label =' Luminosity Contribution',linestyle='solid')
ax2.set_title('Gyr = 5')
ax2.set_xlabel('Mass (Msun)',fontweight='bold')

ax3.set_title('Myr = 100')
ax3.plot(Mass_myr100[1:],mass_myr100_normalized, color = 'k', label = 'Mass Contribution',linestyle='solid')
ax3.set_xlabel('Mass (Msun)',fontweight='bold')
ax3.plot(Mass_myr100[1:],luminosity_myr100_normalized, color = 'r', label =' Luminosity Contribution',linestyle='solid')

plt.legend(loc='best', bbox_to_anchor=(.5, -.15),
          fancybox=True, shadow=True, ncol=2,fontsize = 'small')
plt.show()




''' Problem 3: Predicted Spectrum'''

# wavelength range
wavelenth = np.arange(1,19500, 100) # angstroms
wavelenth= wavelenth*1e-10 # m
# Planck Function
def Planck_function(lam,T):
    '''Planck Function.; 
    Input: h = Planck's constant, c = speed of light, T = temperature, 
    lam = wavelength, k = Boltzmann's constant.
    Output: the flux at a given wavelength and temperature. '''
    B_numerator = (2*h*c*c)/(lam**5)
    B_denominator = np.exp((h*c)/(lam*k*T))-1
    planck = B_numerator/B_denominator
    return planck
#    
    

# Flux = integral f_Star * IMF f_star = black body, IMF = original imf (from definition)
# obtain ttmperure range range
 # two for loops, one trhorugh the wavelgntha dn one through the tmeprature, inlude mass rnage
 # do the quad stuff  limits of integrations is mass range as in first for loop

# Flux for Gyr 1
flux_gyr1 = []
sum_gyr1 = []
for i in range(0,len(wavelenth)):
    expon =  10**(Temperature_gyr1[i])
    planck_values = Planck_function(wavelenth[i],expon)
    flux_gyr1.append(planck_values)
    f = planck_values*IMF(Mass_gyr1[i])
    sum_gyr1.append(np.sum(f))
    
# Flux for Gyr 5  
flux_gyr5 = []
sum_gyr5 = []
for i in range(0,len(wavelenth)):
    expon =  10**(Temperature_gyr5[i])
    planck_values = Planck_function(wavelenth[i],expon)
    flux_gyr5.append(planck_values)
    f = planck_values*IMF(Mass_gyr5[i])
    sum_gyr5.append(np.sum(f))
    

# Flux for Myr 100   
flux_myr500 = []
sum_myr100 = []
for i in range(0,len(wavelenth)):
    expon =  10**(Temperature_myr100[i])
    planck_values = Planck_function(wavelenth[i],expon)
    flux_myr500.append(planck_values)
    f = planck_values*IMF(Mass_myr100[i])
    sum_myr100.append(np.sum(f))
#    sum_gry1.append(flux_gyr1*IMF(Mass_gyr1[i]))
    

# Plots the predicted spectrums
plt.figure()
plt.plot(wavelenth/1e-10, sum_gyr1, label='Gyr 1')
plt.plot(wavelenth/1e-10, sum_gyr5, label='Gyr 5')
plt.plot(wavelenth/1e-10,sum_myr100, label='Myr 100')
plt.title('Predicted Spectra',fontsize=14,fontweight='bold')
plt.xlabel('Wavelength (Å)',fontweight='bold')
plt.ylabel('Relative Flux',fontweight='bold')
plt.legend(loc='best',
          fancybox=True, shadow=True)

plt.tight_layout()
plt.show()


    
    
    
    
    
    
    
    
