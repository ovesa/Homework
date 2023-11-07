t#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 14:07:21 2018

These python routines are for ASTR 545 Homework 1. It is due on Monday 28, 2018.

@author: oanavesa
"""
#%%
# import necessary modules
import numpy as np
import matplotlib.pyplot as plt
import glob

#%% Problem 2

''' Problem 2 '''

# defining path and text files 
path = '/Users/oanavesa/Desktop/GradSchool/FirstYear/StellarSpectroscopy/HW/HW1/'
star1a = 'star1-Ha.txt'
star1b= 'star1-Hb.txt'
star2a = 'star2-Ha.txt'
star2b= 'star2-Hb.txt'
star3a = 'star3-Ha.txt'
star3b= 'star3-Hb.txt'
star4a = 'star4-Ha.txt'
star4b= 'star4-Hb.txt'

# upacking all of the columns
wavelengthvaluestar1a, fluxvaluestar1a, uncertaintyvaluestar1a = np.loadtxt(path+star1a, unpack=True,usecols=[0,1,2])
wavelengthvaluestar1b, fluxvaluestar1b, uncertaintyvaluestar1b = np.loadtxt(path+star1b, unpack=True,usecols=[0,1,2])
wavelengthvaluestar2a, fluxvaluestar2a, uncertaintyvaluestar2a = np.loadtxt(path+star2a, unpack=True,usecols=[0,1,2])
wavelengthvaluestar2b, fluxvaluestar2b, uncertaintyvaluestar2b = np.loadtxt(path+star2b, unpack=True,usecols=[0,1,2])
wavelengthvaluestar3a, fluxvaluestar3a, uncertaintyvaluestar3a = np.loadtxt(path+star3a, unpack=True,usecols=[0,1,2])
wavelengthvaluestar3b, fluxvaluestar3b, uncertaintyvaluestar3b = np.loadtxt(path+star3b, unpack=True,usecols=[0,1,2])
wavelengthvaluestar4a, fluxvaluestar4a, uncertaintyvaluestar4a = np.loadtxt(path+star4a, unpack=True,usecols=[0,1,2])
wavelengthvaluestar4b, fluxvaluestar4b, uncertaintyvaluestar4b = np.loadtxt(path+star4b, unpack=True,usecols=[0,1,2])

# cycles through all of text files in the directory
fname = glob.glob(path+'*.txt')

# plotting the absorption line of all the stars
for i in range(len(fname)):
    lam,flux,unc = np.genfromtxt(fname[i], unpack=True,usecols=[0,1,2])
    plt.title(fname[i][72:80])
    plt.plot(lam,flux)
    plt.ylabel(fname[i][78:80])
    plt.xlabel('$\AA$')
    plt.show()
    
    
    
''' 
    Choosing the absorption profile and obtaining the equivalent widths and 
    uncertainties.
    
    Looping through all of the uncertainty values per star that have an uncertainty
    smaller than 0.0196 Å--where it seems that the uncertainty starts dropping by a lot
    in the text files--indicating where the absorption was occurring. I did not just choose 
    0.0199 Å just because I thought that that was too close to the wings. I also looked at 
    the approximate wavelength on the plots where the dip started happening and chose a point
    around there.
    
'''    


# define the equivalent width equation
# col1 = the wavelength in pixel j
# col2 = the relative flux in pixel j
# col3 = the uncertainty at pixel j
# Notes: just want a single value for the equivalent width and uncertainty. Don't use array/list



''' Star 1-Hα '''

# obtaining equivalent width
width = 0
for i in np.where(uncertaintyvaluestar1a<0.0196):
    W = (1-fluxvaluestar1a[i])*(wavelengthvaluestar1a[2]-wavelengthvaluestar1a[1])
    width = np.sum(W)
#print(width)
     
# obtaining uncertainty   
uncertainty = 0
for j in range(len(uncertaintyvaluestar1a[np.where(uncertaintyvaluestar1a<0.0196)])):
    uc = (uncertaintyvaluestar1a[i]*(wavelengthvaluestar1a[2]-wavelengthvaluestar1a[1]))**2
    uncertainty  = np.sum(uc)
    uncertsq = np.sqrt(uncertainty)
#print(uncertsq)



''' Srar 1-Hβ '''

# obtaining equivalent width
width = 0
for i in np.where(uncertaintyvaluestar1b<.0196):
    W = (1-fluxvaluestar1b[i])*(wavelengthvaluestar1b[2]-wavelengthvaluestar1b[1])
    width = np.sum(W)
#print(width)
     
# obtaining uncertainty    
uncertainty = 0
for j in range(len(uncertaintyvaluestar1b[np.where(uncertaintyvaluestar1b<.0196)])):
    uc = (uncertaintyvaluestar1b[i]*(wavelengthvaluestar1b[2]-wavelengthvaluestar1b[1]))**2
    uncertainty  = np.sum(uc)
    uncertsq = np.sqrt(uncertainty)
#print(uncertsq)



''' Star 2-Hα '''

# obtaining equivalent width
width = 0
for i in np.where(uncertaintyvaluestar2a<.01979):
    W = (1-fluxvaluestar2a[i])*(wavelengthvaluestar2a[2]-wavelengthvaluestar2a[1])
    width = np.sum(W)
#print(width)
     
# obtaining uncertainty       
uncertainty = 0
for j in range(len(uncertaintyvaluestar2a[np.where(uncertaintyvaluestar2a<.01979)])):
    uc = (uncertaintyvaluestar2a[i]*(wavelengthvaluestar2a[2]-wavelengthvaluestar2a[1]))**2
    uncertainty  = np.sum(uc)
    uncertsq = np.sqrt(uncertainty)
print(uncertsq)



''' Srar 2-Hβ '''

# obtaining equivalent width
width = 0
for i in np.where(uncertaintyvaluestar2b<.0196):
    W = (1-fluxvaluestar2b[i])*(wavelengthvaluestar2b[2]-wavelengthvaluestar2b[1])
    width = np.sum(W)
#print(width)
     
# obtaining uncertainty       
uncertainty = 0
for j in range(len(uncertaintyvaluestar2b[np.where(uncertaintyvaluestar2b<.0196)])):
    uc = (uncertaintyvaluestar2b[i]*(wavelengthvaluestar2b[2]-wavelengthvaluestar2b[1]))**2
    uncertainty  = np.sum(uc)
    uncertsq = np.sqrt(uncertainty)
#print(uncertsq)



''' Star 3-Hα '''

# obtaining equivalent width
width = 0
for i in np.where(uncertaintyvaluestar3a<.0196):
    W = (1-fluxvaluestar3a[i])*(wavelengthvaluestar3a[2]-wavelengthvaluestar3a[1])
    width = np.sum(W)
#print(width)
     
# obtaining uncertainty       
uncertainty = 0
for j in range(len(uncertaintyvaluestar3a[np.where(uncertaintyvaluestar3a<.0196)])):
    uc = (uncertaintyvaluestar3a[i]*(wavelengthvaluestar3a[2]-wavelengthvaluestar3a[1]))**2
    uncertainty  = np.sum(uc)
    uncertsq = np.sqrt(uncertainty)
#print(uncertsq)



''' Srar 3-Hβ '''

# obtaining equivalent width
width = 0
for i in np.where(uncertaintyvaluestar3b<.0196):
    W = (1-fluxvaluestar3b[i])*(wavelengthvaluestar3b[2]-wavelengthvaluestar3b[1])
    width = np.sum(W)
#print(width)
     
# obtaining uncertainty        
uncertainty = 0
for j in range(len(uncertaintyvaluestar3b[np.where(uncertaintyvaluestar3b<.0196)])):
    uc = (uncertaintyvaluestar3b[i]*(wavelengthvaluestar3b[2]-wavelengthvaluestar3b[1]))**2
    uncertainty  = np.sum(uc)
    uncertsq = np.sqrt(uncertainty)
#print(uncertsq)



''' Star 4-Hα '''

# obtaining equivalent width
width = 0
for i in np.where(uncertaintyvaluestar4a<.0196):
    W = (1-fluxvaluestar4a[i])*(wavelengthvaluestar4a[2]-wavelengthvaluestar4a[1])
    width = np.sum(W)
#print(width)
     
# obtaining uncertainty        
uncertainty = 0
for j in range(len(uncertaintyvaluestar4a[np.where(uncertaintyvaluestar4a<.0196)])):
    uc = (uncertaintyvaluestar4a[i]*(wavelengthvaluestar4a[2]-wavelengthvaluestar4a[1]))**2
    uncertainty  = np.sum(uc)
    uncertsq = np.sqrt(uncertainty)
#print(uncertsq)



''' Srar 4-Hβ '''

# obtaining equivalent width
width = 0
for i in np.where(uncertaintyvaluestar4b<.0196):
    W = (1-fluxvaluestar4b[i])*(wavelengthvaluestar4b[2]-wavelengthvaluestar4b[1])
    width = np.sum(W)
#print(width)
     
# obtaining uncertainty        
uncertainty = 0
for j in range(len(uncertaintyvaluestar4b[np.where(uncertaintyvaluestar4b<.0196)])):
    uc = (uncertaintyvaluestar4b[i]*(wavelengthvaluestar4b[2]-wavelengthvaluestar4b[1]))**2
    uncertainty  = np.sum(uc)
    uncertsq = np.sqrt(uncertainty)
#print(uncertsq)

#%% Problem 3

''' Problem 3 '''

# defining constants
T = 10000 # (K)
k = 1.38e-16 # (erg/K)
h = 6.626e-27 # (erg*s)
c = 3e10 # (cm/s)
lam = 160 # (cm)
lamang = lam/1e-8 # lam in (Å)
dA = 0.1 # (Å/pixel) = spacing
q = 8*np.pi*h*c*c # (erg*cm^2/s)
hcvalue = h*c # (erg*cm)
lam1 = 2600 # (Å)
lam2 = 7100 # (Å)


# wavelength range: 2600–7100 Å
wavang= np.arange(lam1,lam2,dA) # (Å)
wavcm = np.arange(lam1*1e-8,lam2*1e-8,dA*1e-8) # (cm)
# creates an empty list for the incident spectrum
incidentblackbody = []
# creates a second empty list the size of the wavelength array in angstroms in order to get the shapes to match later on
incidentblackbody_sameshape = np.zeros(wavang.shape) 


''' Part 2a: Computing the incident spectrum over the wavelength range 2600–7100 Å '''


# calculate incident spectrum using Blackbody flux given
# split up units for it to work: lam^4 in cm and lam in Å
for i in range(len(wavang)): 
    blackbodyspecfront = (q/(wavcm[i]**4))*(1/(wavang[i]))
    blackbodyspecsecond = np.exp(hcvalue/(wavcm[i]*k*T))-1
    blackbody = blackbodyspecfront*(1/blackbodyspecsecond)
   #print(blackbody)
    incidentblackbody = np.append(incidentblackbody,blackbody)
# copies the values calculated above from incidentblackbody to a different array whose size matches the size of the
# wavelength array in angstroms, wavang
incidentblackbody_sameshape[:] = incidentblackbody[:]
  

''' Part 2b: Computing the transmitted spectrum for a specific cross section and for wavelength values less than 3646 Å '''


# defining more constants
numdensity =  3e12 #(cm^(-3))
lamnot = 3646 # (Å)


alpha = []
# finds the indices in the wavelength array in angstroms where the value is less than 3636 Å
specificvalues = np.where(wavang <= lamnot)[0]
# uses the above code to find the exact values that match those indices
wav = wavang[specificvalues]
# computes the cross section given in the homework assignment
alphas = (2.1e-15)*((wavang[specificvalues])/lamnot)**3


''' Need to use the equation of radiative transfer: I = I0*exp(-tau), where tau = n*L*α '''


# computing arrays to be the same size overall so that I can plot them


# creates an array of 0's the same size as the difference between the wavelength array in angstroms and the
# cross section array for later use
differencearray = np.zeros((wavang.shape[0]-alphas.shape[0]))
# concatenates together the cross section array and the above to get the cross section array to match the wavelength array in angstrojs
# so that it can be plotted
newalphas = np.concatenate((alphas, differencearray))
# creates a 0's array the same size as the cross section array
tau_real = np.zeros(newalphas.shape)
# computes an array of 0's the same shape as the cross section array
exptausamesize = np.zeros((newalphas.shape))


# creates an empty list to store that computed tau values in
tau2 = []
inc = []
# for loop computes the tau values over the new cross section for the radiation transfer equation
for j in range(len(newalphas)):
    tau = numdensity*newalphas[j]*((lam)) #alpha*n*l
    tau2 = np.append(tau2,tau)
# stores the values calculated above into an array of 0's to get the shapes to match up in the end
tau_real[:] = tau2[:]
# gives negative values of the tau values to store into the radiative transfer equation later
negtau = -1*tau_real


# an empty array to store the values obtained by taking the exponent on the tau values
exptaurealvalues = []
# for loop takes the exponent of each tau values
for a in range(len(negtau)):
    exptau = np.exp(negtau[a])
    exptaurealvalues = np.append(exptaurealvalues,exptau)
# stores the values obtained above into an array that matches the size of the cross section array
exptausamesize[:] =  exptaurealvalues[:]


# calculates the transmitted spectrum (the radiative transfer equation)
transmittedspectrum = incidentblackbody_sameshape*exptausamesize


''' Plotting the incident spectrum (solid curve) and the transmitted spectrum (short-dash curve) '''
   
 
plt.figure()
plt.plot(wavang,incidentblackbody_sameshape,label='Incident Spectrum')
plt.plot(wavang,transmittedspectrum[:],'--', label='Transmitted Spectrum')
plt.xlabel('Wavelength (Å)', fontsize=12, fontweight='bold')
plt.ylabel(r'Blackbody Flux $(erg \cdot s^{-1} cm^{-1} Å^{-1})$', fontsize=12, fontweight='bold')
plt.title('Ionization Break', fontweight='bold',fontsize=13)
plt.legend(loc='upper right')
plt.plot()
plt.show()
