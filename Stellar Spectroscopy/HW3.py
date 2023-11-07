r#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ASTR 545: Homework 3

Using the four Hα and Hβ spetra from HW 1, find the best (b,N) pair that minimizes Chi^2. Plot
b vs Chi^2 as a function of the best log(N) for each star. Creating an absorption profile for
Na I D_2 without pressure broadening and with pressure broadening effects.

Created on Mon Sep 10 12:45:50 2018


@author: oanavesa
"""
#%%
# necessary imports
import numpy as np
import scipy as sp
from scipy import special
import matplotlib.pyplot as plt
#%%


# defining path and text files  from HW 1
path = '/Users/oanavesa/Desktop/GradSchool/FirstYear/StellarSpectroscopy/HW/HW1/'
star1a = 'star1-Ha.txt'
star1b= 'star1-Hb.txt'
star2a = 'star2-Ha.txt'
star2b= 'star2-Hb.txt'
star3a = 'star3-Ha.txt'
star3b= 'star3-Hb.txt'
star4a = 'star4-Ha.txt'
star4b= 'star4-Hb.txt'

# upacking all of the columns: wavelength, flux, and uncertainty values
wavelengthvaluestar1a, fluxvaluestar1a, uncertaintyvaluestar1a = np.loadtxt(path+star1a, unpack=True,usecols=[0,1,2])
wavelengthvaluestar1b, fluxvaluestar1b, uncertaintyvaluestar1b = np.loadtxt(path+star1b, unpack=True,usecols=[0,1,2])
wavelengthvaluestar2a, fluxvaluestar2a, uncertaintyvaluestar2a = np.loadtxt(path+star2a, unpack=True,usecols=[0,1,2])
wavelengthvaluestar2b, fluxvaluestar2b, uncertaintyvaluestar2b = np.loadtxt(path+star2b, unpack=True,usecols=[0,1,2])
wavelengthvaluestar3a, fluxvaluestar3a, uncertaintyvaluestar3a = np.loadtxt(path+star3a, unpack=True,usecols=[0,1,2])
wavelengthvaluestar3b, fluxvaluestar3b, uncertaintyvaluestar3b = np.loadtxt(path+star3b, unpack=True,usecols=[0,1,2])
wavelengthvaluestar4a, fluxvaluestar4a, uncertaintyvaluestar4a = np.loadtxt(path+star4a, unpack=True,usecols=[0,1,2])
wavelengthvaluestar4b, fluxvaluestar4b, uncertaintyvaluestar4b = np.loadtxt(path+star4b, unpack=True,usecols=[0,1,2])

# converting all the wavelengths from Å to cm
wavelengthvaluestar1a=wavelengthvaluestar1a/1e+8
wavelengthvaluestar1b=wavelengthvaluestar1b/1e+8
wavelengthvaluestar2a=wavelengthvaluestar2a/1e+8
wavelengthvaluestar2b=wavelengthvaluestar2b/1e+8
wavelengthvaluestar3a=wavelengthvaluestar3a/1e+8
wavelengthvaluestar3b=wavelengthvaluestar3b/1e+8
wavelengthvaluestar4a=wavelengthvaluestar4a/1e+8
wavelengthvaluestar4b=wavelengthvaluestar4b/1e+8
#%%

# Defining constants
me = 9.1e-28  # g = mass of electron
c = (3e8)*100 # cm
el = (4.803e-10) # ESU = g1/2 cm3/2 s-1

#atomic constants
#Hα
lam_alpha = 6564.623/1e+8 # cm = transition wavelength
f_alpha = 0.6958 #  oscillator strength
gamma_alpha = 6.465e7 # 1/s = damping constant
#Hβ
lam_beta = 4862.688/1e+8 # cm = transition wavelength
f_beta = 0.1218 # oscillator strength
gamma_beta = 2.062e7 # 1/s = damping constant

# arrays for the b and N search grid
b = np.arange(1e+6,7.5e+6,500000) # cm/s # Doppler parameter (13 values)
n = np.arange(12,22.0,.1) # log N = column density (100 values)
 
# total number of pixels applied in the sums 
M = wavelengthvaluestar1a.shape[0] + wavelengthvaluestar1a.shape[0]
# the number of free parameters to be estimated
m = 2
v = M-m
#%% Star 1

''' Problem 1: Finding the (b,N) pair that minimizes Chi^2. Part (a)
    Star 1'''

# empty lists to store b, N, and chi^2 values
bval = []
nval=[]
cb = []

# search grid for b values
for i in range(0,len(b)):
    # Doppler width
    cda = (lam_alpha/c)*b[i]
    cdb = (lam_beta/c)*b[i]
    
    # y values for α and β, respectively
    yvala = (1/cda)*((gamma_alpha*(lam_alpha)**2)/(4*np.pi*c))
    yvalb = (1/cdb)*((gamma_beta*(lam_beta)**2)/(4*np.pi*c))
    
    # search grid for N values
    for p in range(0,len(n)):
        # empty lists for the fluxes and sums
        sumalphabs = []
        sumbetabs = []
        fluxalpha1 = []   
        fluxbeta1 = []  
        
        # looping through the indices of the wavelength of Star 1
        for g in range(0,len(wavelengthvaluestar1a)):
            # x-values for Hα
            xvala = ((wavelengthvaluestar1a[g]-lam_alpha)/(cda))
            # voigt profile
            voigta = np.real(sp.special.wofz(xvala+1j*yvala))
            # Alpha calculation
            alpha = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_alpha*lam_alpha*lam_alpha)/(cda))*voigta)
            # Relative flux: I/I0 = exp(-tau), where tau = n*alpha*L
            relativefluxa = np.exp(-1*alpha*(10**n[p]))
            fluxalpha1.append(relativefluxa)
            
            # x-values for Hβ
            xvalb = ((wavelengthvaluestar1b[g] - lam_beta)/(cdb))
            # voigt profile
            voigtb = np.real(sp.special.wofz(xvalb+1j*yvalb))
            # Alpha calculation
            alphab = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_beta*lam_beta*lam_beta)/(cdb))*voigtb)
            # Relative flux: I/I0 = exp(-tau), where tau = n*alpha*L
            relativefluxb = np.exp(-1*alphab*(10**n[p]))
            fluxbeta1.append(relativefluxb)
            
            # Chi^2 calculations
            # Hα chis
            v2 = 1.0/v
            alphabs = ((fluxvaluestar1a[g]-relativefluxa)/uncertaintyvaluestar1a[g])
            sqralpha = alphabs**2
            sumalphabs.append(sqralpha)
            
            # Hβ chis
            betabs = ((fluxvaluestar1b[g]-relativefluxb)/uncertaintyvaluestar1b[g])
            sqrbeta = betabs**2
            sumbetabs.append(sqrbeta)

        # computing Chi^2
        chisq = v2*(np.sum(sumalphabs) +np.sum(sumbetabs))
        cb.append(chisq)
        bval.append(b[i])
        nval.append(n[p])

# obtaining the minimum value for chi^2
minchistar1 = min(cb)
# find that corresponding index
indstar1 = np.where(cb==minchistar1)[0][0]
# obtaining the b value for the minimum chi^2 value
bvalstar1 = bval[indstar1]
# obtaining the N value for the minimum chi^2 value
nvalstar1 = nval[indstar1]
print(bval[indstar1]/100000)
print(nval[indstar1])
# The (b,N) pair is (30,12.5)
#%%

''' Star 2'''

# empty lists to store b, N, and chi^2 values
bval2 = []
nval2=[]
cb2 = []

# search grid for b values
for i in range(0,len(b)):
    # Doppler width
    cda = (lam_alpha/c)*b[i]
    cdb = (lam_beta/c)*b[i]
    
    # y values for α and β, respectively
    yvala = (1/cda)*((gamma_alpha*(lam_alpha)**2)/(4*np.pi*c))
    yvalb = (1/cdb)*((gamma_beta*(lam_beta)**2)/(4*np.pi*c))
    
    # search grid for N values
    for p in range(0,len(n)):
        # empty lists for the fluxes and sums
        sumalphabs2 = []
        sumbetabs2 = []
        fluxalpha2 = []   
        fluxbeta2 = []  
    
        # looping through the indices of the wavelength of Star 1
        for g in range(0,len(wavelengthvaluestar2a)):
            # x-values for Hα
            xvala = ((wavelengthvaluestar2a[g]-lam_alpha)/(cda))
            # voigt profile
            voigta = np.real(sp.special.wofz(xvala+1j*yvala))
            # Alpha calculation
            alpha = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_alpha*lam_alpha*lam_alpha)/(cda))*voigta)
            # Relative flux: I/I0 = exp(-tau), where tau = n*alpha*L
            relativefluxa = np.exp(-1*alpha*(10**n[p]))
            fluxalpha2.append(relativefluxa)
            
            # x-values for Hβ
            xvalb = ((wavelengthvaluestar2b[g] - lam_beta)/(cdb))
            # voigt profile
            voigtb = np.real(sp.special.wofz(xvalb+1j*yvalb))
            # Alpha calculation
            alphab = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_beta*lam_beta*lam_beta)/(cdb))*voigtb)
            # Relative flux: I/I0 = exp(-tau), where tau = n*alpha*L
            relativefluxb = np.exp(-1*alphab*(10**n[p]))
            fluxbeta2.append(relativefluxb)
            
            # Chi^2 calculations
            # Hα chis
            v2 = 1.0/v
            alphabs = ((fluxvaluestar2a[g]-relativefluxa)/uncertaintyvaluestar2a[g])
            sqralpha = alphabs**2
            sumalphabs2.append(sqralpha)
            
            # Hβ chis
            betabs = ((fluxvaluestar2b[g]-relativefluxb)/uncertaintyvaluestar2b[g])
            sqrbeta = betabs**2
            sumbetabs2.append(sqrbeta)

        # Computing chi^2
        chisq = v2*(np.sum(sumalphabs2) +np.sum(sumbetabs2))
        cb2.append(chisq)
        bval2.append(b[i])
        nval2.append(n[p])

# obtaining the minimum value for chi^2
minchistar2 = min(cb2)
# find that corresponding index
indstar2 = np.where(cb2==minchistar2)[0][0]
# obtaining the b value for the minimum chi^2 value
bvalstar2 = bval2[indstar2]
# obtaining the N value for the minimum chi^2 value
nvalstar2 = nval2[indstar2]
print(bval2[indstar2]/100000)
print(nval2[indstar2])
# The (b,N) pair is (55,13.7)
#%%


''' Star 3 '''

# empty lists to store b, N, and chi^2 values
bval3 = []
nval3=[]
cb3 = []

# search grid for b values
for i in range(0,len(b)):
    # Doppler width
    cda = (lam_alpha/c)*b[i]
    cdb = (lam_beta/c)*b[i]
    
    # y values for α and β, respectively
    yvala = (1/cda)*((gamma_alpha*(lam_alpha)**2)/(4*np.pi*c))
    yvalb = (1/cdb)*((gamma_beta*(lam_beta)**2)/(4*np.pi*c))
    
    # search grid for N values
    for p in range(0,len(n)):
        # empty lists for the fluxes and sums
        sumalphabs3 = []
        sumbetabs3 = []
        fluxalpha3 = []   
        fluxbeta3 = []  
    
        for g in range(0,len(wavelengthvaluestar3a)):
            # x-values for Hα
            xvala = ((wavelengthvaluestar3a[g]-lam_alpha)/(cda))
            # voigt profile
            voigta = np.real(sp.special.wofz(xvala+1j*yvala))
            # Alpha calculation
            alpha = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_alpha*lam_alpha*lam_alpha)/(cda))*voigta)
            # Relative flux: I/I0 = exp(-tau), where tau = n*alpha*L 
            relativefluxa = np.exp(-1*alpha*(10**n[p]))
            fluxalpha3.append(relativefluxa)
            
            # x-values for Hβ
            xvalb = ((wavelengthvaluestar3b[g] - lam_beta)/(cdb))
            # voigt profile
            voigtb = np.real(sp.special.wofz(xvalb+1j*yvalb))
            # Alpha calculation
            alphab = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_beta*lam_beta*lam_beta)/(cdb))*voigtb)
            # Relative flux: I/I0 = exp(-tau), where tau = n*alpha*L 
            relativefluxb = np.exp(-1*alphab*(10**n[p]))
            fluxbeta3.append(relativefluxb)
            
            # Chi^2 calculations
            # Hα chis
            v2 = 1.0/v
            alphabs = ((fluxvaluestar3a[g]-relativefluxa)/uncertaintyvaluestar3a[g])
            sqralpha = alphabs**2
            sumalphabs3.append(sqralpha)
            
            # Hβ chis
            betabs = ((fluxvaluestar3b[g]-relativefluxb)/uncertaintyvaluestar3b[g])
            sqrbeta = betabs**2
            sumbetabs3.append(sqrbeta)

        # Computing chi^2
        chisq = v2*(np.sum(sumalphabs3) +np.sum(sumbetabs3))
        cb3.append(chisq)
        bval3.append(b[i])
        nval3.append(n[p])

# obtaining the minimum value for chi^2
minchistar3 = min(cb3)
# find that corresponding index
indstar3 = np.where(cb3==minchistar3)[0][0]
# obtaining the b value for the minimum chi^2 value
bvalstar3 = bval3[indstar3]
# obtaining the N value for the minimum chi^2 value
nvalstar3 = nval3[indstar3]
print(bval3[indstar3]/100000)
print(nval3[indstar3])
# The (b,N) pair is (25,14.2)
#%%

''' Star 4 '''

# empty lists to store b, N, and chi^2 values
bval4 = []
nval4=[]
cb4 = []
emptyl = []
# search grid for b values
for i in range(0,len(b)):
    # Doppler width
    cda = (lam_alpha/c)*b[i]
    cdb = (lam_beta/c)*b[i]
    
    # y values for α and β, respectively
    yvala = (1/cda)*((gamma_alpha*(lam_alpha)**2)/(4*np.pi*c))
    yvalb = (1/cdb)*((gamma_beta*(lam_beta)**2)/(4*np.pi*c))
    
    # search grid for N values
    for p in range(0,len(n)):
        # empty lists for the fluxes and sums
        sumalphabs4 = []
        sumbetabs4 = []
        fluxalpha4 = []   
        fluxbeta4 = []  
    
        for g in range(0,len(wavelengthvaluestar4a)):
            # x-values for Hα
            xvala = ((wavelengthvaluestar4a[g]-lam_alpha)/(cda))
#            print(xvala)
            # voigt profile
            voigta = np.real(sp.special.wofz(xvala+1j*yvala))
            # Alpha calculation
            alpha = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_alpha*lam_alpha*lam_alpha)/(cda))*voigta)
            taua = -1*alpha*(10**n[p])
            #print(taua)
            # Relative flux: I/I0 = exp(-tau), where tau = n*alpha*L 
            relativefluxa = np.exp(-1*alpha*(10**n[p]))
            
#            print(relativefluxa)
            fluxalpha4.append(relativefluxa)
            
            # x-values for Hβ
            xvalb = ((wavelengthvaluestar4b[g] - lam_beta)/(cdb))
            # voigt profile
            voigtb = np.real(sp.special.wofz(xvalb+1j*yvalb))
            # Alpha calculation
            alphab = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_beta*lam_beta*lam_beta)/(cdb))*voigtb)
            # Relative flux: I/I0 = exp(-tau), where tau = n*alpha*L 
#            print(taua)
            relativefluxb = np.exp(-1*alphab*(10**n[p]))
            fluxbeta4.append(relativefluxb)
            
            # Chi^2 calculations
            # Hα chis
            v2 = 1.0/v
            alphabs = ((fluxvaluestar4a[g]-relativefluxa)/uncertaintyvaluestar4a[g])
            sqralpha = alphabs**2
            sumalphabs4.append(sqralpha)
            
            # Hβ chis
            betabs = ((fluxvaluestar4b[g]-relativefluxb)/uncertaintyvaluestar4b[g])
            sqrbeta = betabs**2
            sumbetabs4.append(sqrbeta)

        # Computing chi^2
        chisq = v2*(np.sum(sumalphabs4) +np.sum(sumbetabs4))
        cb4.append(chisq)
        bval4.append(b[i])
        nval4.append(n[p])

# obtaining the minimum value for chi^2
minchistar4 = min(cb4)
# find that corresponding index
indstar4 = np.where(cb4==minchistar4)[0][0]
# obtaining the b value for the minimum chi^2 value
bvalstar4 = bval4[indstar4]
# obtaining the N value for the minimum chi^2 value
nvalstar4 = nval4[indstar4]
print(bval4[indstar4]/100000)
print(nval4[indstar4])
# The (b,N) pair is (60,17.3)
#%%
''' Problem 1 part b: computing the gas temperature of each star '''

mion = 1.67e-24 # g # hydrogen mass
k = 1.38e-16 # cm2 g s-2 K-1 # Boltzmann constant
c = (3e8)*100 # cm/s

def temp(b,mion,k):
    T = ((b**2)*mion)/(2*k)
    return T

tempstar1 = temp(bvalstar1,mion,k)
tempstar2 = temp(bvalstar2,mion,k)
tempstar3 = temp(bvalstar3,mion,k)
tempstar4 = temp(bvalstar4,mion,k)
print(tempstar1,tempstar2,tempstar3,tempstar4)
#%%

''' Problem 1 part c: For each star, over plot your χ2ν curves for each b as 
a function of your best log(N) '''

''' Star 4 - recopied code from part b. I just got rid of the n loop (column density
because I am just looking at the best log(N) for each star. Therefore, will only get
13 log(N) values'''

# empty lists
cb4 = []

# going through the b array (13 values)
for i in range(0,len(b)):
    cda = (lam_alpha/c)*b[i]
    cdb = (lam_beta/c)*b[i]
    
    yvala = (1/cda)*((gamma_alpha*(lam_alpha)**2)/(4*np.pi*c))
    yvalb = (1/cdb)*((gamma_beta*(lam_beta)**2)/(4*np.pi*c))
    
    # creating empty lists for the chi^2 calculations
    sumalphabs4 = []
    sumbetabs4 = []
    fluxalpha4 = []   
    fluxbeta4 = []  

    for g in range(0,len(wavelengthvaluestar4a)):
        xvala = ((wavelengthvaluestar4a[g]-lam_alpha)/(cda))
        voigta = np.real(sp.special.wofz(xvala+1j*yvala))
        alpha = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_alpha*lam_alpha*lam_alpha)/(cda))*voigta)
        relativefluxa = np.exp(-1*alpha*(10**17.3))
        fluxalpha4.append(relativefluxa)
                
        xvalb = ((wavelengthvaluestar4b[g] - lam_beta)/(cdb))
        voigtb = np.real(sp.special.wofz(xvalb+1j*yvalb))
        alphab = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_beta*lam_beta*lam_beta)/(cdb))*voigtb)
        relativefluxb = np.exp(-1*alphab*(10**17.3))
        fluxbeta4.append(relativefluxb)

        v2 = 1.0/v
        alphabs = ((fluxvaluestar4a[g]-relativefluxa)/uncertaintyvaluestar4a[g])
        sqralpha = alphabs**2
        sumalphabs4.append(sqralpha)
        
        betabs = ((fluxvaluestar4b[g]-relativefluxb)/uncertaintyvaluestar4b[g])
        sqrbeta = betabs**2
        sumbetabs4.append(sqrbeta)

    chisq = v2*(np.sum(sumalphabs4) +np.sum(sumbetabs4))
    cb4.append(chisq)

plt.figure()
plt.title('Star 4: b vs. $\chi^2$ for $\log(N) = 17.3$', fontweight='bold' )
plt.plot(b/100000,cb4)
plt.xlabel('Doppler Parameter (km/s)', fontweight='bold')
plt.ylabel(r'$\chi^2$ Values', fontweight='bold' )
plt.savefig('Star4_b_vs._best_log_N=17.3.png')
plt.show()


#%%

''' Star 3 - recopied code from part b. I just got rid of the n loop (column density
because I am just looking at the best log(N) for each star. Therefore, will only get
13 log(N) values'''


cb3 = []
for i in range(0,len(b)):
    cda = (lam_alpha/c)*b[i]
    cdb = (lam_beta/c)*b[i]
    
    yvala = (1/cda)*((gamma_alpha*(lam_alpha)**2)/(4*np.pi*c))
    yvalb = (1/cdb)*((gamma_beta*(lam_beta)**2)/(4*np.pi*c))
    
    sumalphabs3 = []
    sumbetabs3 = []
    fluxalpha3 = []   
    fluxbeta3 = []  

    for g in range(0,len(wavelengthvaluestar3a)):
        xvala = ((wavelengthvaluestar3a[g]-lam_alpha)/(cda))
        voigta = np.real(sp.special.wofz(xvala+1j*yvala))
        alpha = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_alpha*lam_alpha*lam_alpha)/(cda))*voigta)
        relativefluxa = np.exp(-1*alpha*(10**14.2))
        fluxalpha3.append(relativefluxa)  
        
        xvalb = ((wavelengthvaluestar3b[g] - lam_beta)/(cdb))
        voigtb = np.real(sp.special.wofz(xvalb+1j*yvalb))
        alphab = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_beta*lam_beta*lam_beta)/(cdb))*voigtb)
        relativefluxb = np.exp(-1*alphab*(10**14.2))
        fluxbeta3.append(relativefluxb)
        
        v2 = 1.0/v
        alphabs = ((fluxvaluestar3a[g]-relativefluxa)/uncertaintyvaluestar3a[g])
        sqralpha = alphabs**2
        sumalphabs3.append(sqralpha)
        
        betabs = ((fluxvaluestar3b[g]-relativefluxb)/uncertaintyvaluestar3b[g])
        sqrbeta = betabs**2
        sumbetabs3.append(sqrbeta)

    chisq = v2*(np.sum(sumalphabs3) +np.sum(sumbetabs3))
    cb3.append(chisq)


plt.figure()
plt.title('Star 3: b vs. $\chi^2$ for $\log(N) = 14.2$', fontweight='bold' )
plt.plot(b/100000,cb3)
plt.xlabel('Doppler Parameter (km/s)', fontweight='bold')
plt.ylabel(r'$\chi^2$ Values', fontweight='bold' )
plt.savefig('Star3_b_vs._best_log_N=14.2.png')
plt.show()

#%%


''' Star 2 - recopied code from part b. I just got rid of the n loop (column density
because I am just looking at the best log(N) for each star. Therefore, will only get
13 log(N) values'''

cb2 = []
for i in range(0,len(b)):  
    cda = (lam_alpha/c)*b[i]
    cdb = (lam_beta/c)*b[i]
    
    yvala = (1/cda)*((gamma_alpha*(lam_alpha)**2)/(4*np.pi*c))
    yvalb = (1/cdb)*((gamma_beta*(lam_beta)**2)/(4*np.pi*c))
    
    sumalphabs2 = []
    sumbetabs2 = []
    fluxalpha2 = []   
    fluxbeta2 = []  
    
    for g in range(0,len(wavelengthvaluestar2a)):
        xvala = ((wavelengthvaluestar2a[g]-lam_alpha)/(cda))
        voigta = np.real(sp.special.wofz(xvala+1j*yvala))
        alpha = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_alpha*lam_alpha*lam_alpha)/(cda))*voigta)
        relativefluxa = np.exp(-1*alpha*(10**13.7))
        fluxalpha2.append(relativefluxa)
        
        xvalb = ((wavelengthvaluestar2b[g] - lam_beta)/(cdb))
        voigtb = np.real(sp.special.wofz(xvalb+1j*yvalb))
        alphab = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_beta*lam_beta*lam_beta)/(cdb))*voigtb)
        relativefluxb = np.exp(-1*alphab*(10**13.7))
        fluxbeta2.append(relativefluxb)
        
        v2 = 1.0/v
        alphabs = ((fluxvaluestar2a[g]-relativefluxa)/uncertaintyvaluestar2a[g])
        sqralpha = alphabs**2
        sumalphabs2.append(sqralpha)
        
        betabs = ((fluxvaluestar2b[g]-relativefluxb)/uncertaintyvaluestar2b[g])
        sqrbeta = betabs**2
        sumbetabs2.append(sqrbeta)

    chisq = v2*(np.sum(sumalphabs2) +np.sum(sumbetabs2))
    cb2.append(chisq)


plt.figure()
plt.title('Star 2: b vs. $\chi^2$ for $\log(N) = 13.7$', fontweight='bold' )
plt.plot(b/100000,cb2)
plt.xlabel('Doppler Parameter (km/s)', fontweight='bold')
plt.ylabel(r'$\chi^2$ Values', fontweight='bold' )
plt.savefig('Star2_b_vs._best_log_N=13.7.png')
plt.show()
#%%

''' Star 1 - recopied code from part b. I just got rid of the n loop (column density
because I am just looking at the best log(N) for each star. Therefore, will only get
13 log(N) values'''


cb=[]
for i in range(0,len(b)):
    cda = (lam_alpha/c)*b[i]
    cdb = (lam_beta/c)*b[i]
    
    yvala = (1/cda)*((gamma_alpha*(lam_alpha)**2)/(4*np.pi*c))
    yvalb = (1/cdb)*((gamma_beta*(lam_beta)**2)/(4*np.pi*c))
    
    sumalphabs = []
    sumbetabs = []
    fluxalpha1 = []   
    fluxbeta1 = []  

    for g in range(0,len(wavelengthvaluestar1a)):
        xvala = ((wavelengthvaluestar1a[g]-lam_alpha)/(cda))
        voigta = np.real(sp.special.wofz(xvala+1j*yvala))
        alpha = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_alpha*lam_alpha*lam_alpha)/(cda))*voigta)
        relativefluxa = np.exp(-1*alpha*(10**12.5))
        fluxalpha1.append(relativefluxa)
         
        xvalb = ((wavelengthvaluestar1b[g] - lam_beta)/(cdb))
        voigtb = np.real(sp.special.wofz(xvalb+1j*yvalb))
        alphab = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f_beta*lam_beta*lam_beta)/(cdb))*voigtb)
        relativefluxb = np.exp(-1*alphab*(10**12.5))
        fluxbeta1.append(relativefluxb) 
        
        v2 = 1.0/v
        alphabs = ((fluxvaluestar1a[g]-relativefluxa)/uncertaintyvaluestar1a[g])
        sqralpha = alphabs**2
        sumalphabs.append(sqralpha)
        
        betabs = ((fluxvaluestar1b[g]-relativefluxb)/uncertaintyvaluestar1b[g])
        sqrbeta = betabs**2
        sumbetabs.append(sqrbeta)

    chisq = v2*(np.sum(sumalphabs) +np.sum(sumbetabs))
    cb.append(chisq)

plt.figure()
plt.title('Star 1: b vs. $\chi^2$ for $\log(N) = 12.5$', fontweight='bold' )
plt.plot(b/100000,cb)
plt.xlabel('Doppler Parameter (km/s)', fontweight='bold')
plt.ylabel(r'$\chi^2$ Values', fontweight='bold' )
plt.savefig('Star1_b_vs._best_log_N=12.5.png')
plt.show()

#%%

''' Problem 2 '''

# constants for 
lam0 = 5891.583/1e+8 # cm = transition wavelength
nij = 0.1591*lam0 # phase shift (cm)
KIJ = 6.76e-16
ma = 22.98977 #(amu) mass of sodium atom 
mag = 3.817540787e-23 # g = mass of sodium atom
gammaa = 6.30e7 #s^-1 damping constant
In = np.pi/2
Tsol = 6000 # K = solar atmosphere temperature
Pe = 50 # dynes cm^-2 = pressure
bturb = 200000 #cm/s # b value
k = 1.38e-16 # cm2 g s-2 K-1
f3 = 0.6550 # oscillator strength
NaN = 13.0 # best column density value


''' Problem 2 part a '''
def findb(k,Tsol,mag,bturb):
    ''' Input: k = Boltzmann constant (cm^2 g s^-2 K^-1), Tsol = solar atmosphere temperature (K),
    mag = mass of a sodium atom (g), bturb = turbulent atmospheric velocity (cm/s)
    
        Output: bline = Doppler parameter of the Na I D_2 line (cm/s) '''
    bthermal = np.sqrt((2*k*Tsol)/mag)
    bline = np.sqrt(bturb**2+bthermal**2)
    bline = bline*1e-5
    return bline


''' Problem 2 part b '''

def gammapressurebroadening(Pe,Tsol,k,KIJ,nij,mag,me):
    ''' Input: Pe = pressure (dynes cm^−2),k = Boltzmann constant (cm^2 g s^-2 K^-1),
        Tsol = solar atmosphere temperature (K),mag = mass of a sodium atom (g),me = mass
        of an electron (g), KIJ = some constant in the pressure broadening equation, 
        nij = phase shift (cm)
    
        Output: total = pressure broadening gamma (1/s)'''
    term1 = (Pe*Tsol**(-(5/6))/k)
    term1 = 2*np.pi*term1
    term2 = (KIJ*In/nij)
    term2 = term2**(2/3)
    term3 = term3 = (8*k*((1/me)+(1/mag))/np.pi)
    term3 = term3**(1/6)
    total = term1*term2*term3
    return total


def totalgamma(gammapressurebroadening,gammaa):
    ''' Input: gammapressurebroadening = definition from above, gammaa = gamma for sodium atom 
        Output: tgam = total damping constant (1/s)'''
    tgam = gammapressurebroadening(Pe,Tsol,k,KIJ,nij,mag,me)**2+gammaa**2
    return np.sqrt(tgam)

#%%
''' Problem 2 part c: Na1  D_2 Line: No Pressure Broadening and Pressured Broadened '''

# wavelength of NaI
wavena = np.arange(5887,5897,.1) # 	Å
wavena = wavena/1e+8 # cm

# empty lists
fluxalphana = []
fluxbetana = []
    
for g in range(0,len(wavena)):
    # Doppler width parameter
    cda = (lam0/c)*288000
    cdb = (lam0/c)*288000
    
    # NaI without pressure broadening
    yvala = (1/cda)*((gammaa*(lam0)**2)/(4*np.pi*c))
    xvala = ((wavena[g]-lam0)/(cda))
    # Voigt profile
    voigta = np.real(sp.special.wofz(xvala+1j*yvala))
    # Alpha calculation
    alpha = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f3*lam0*lam0)/(cda))*voigta)
    # Relative flux: I/I0 = exp(-tau), where tau = n*alpha*L
    relativefluxa = np.exp(-1*alpha*(10**13))
    fluxalphana.append(relativefluxa)
    
    # NaI with pressure broadening
    yvalb = (1/cdb)*((totalgamma(gammapressurebroadening,gammaa)*(lam0)**2)/(4*np.pi*c))
    xvalb = ((wavena[g] - lam0)/(cdb))
    voigtb = np.real(sp.special.wofz(xvalb+1j*yvalb))
    alphab = (((np.sqrt(np.pi)*el*el)/(me*c*c))*((f3*lam0*lam0)/(cdb))*voigtb)
    relativefluxb = np.exp(-1*alphab*(10**13))
    fluxbetana.append(relativefluxb)
    


plt.figure()
plt.title(r"Na I $D_{2}$ Line", fontweight='bold')
plt.plot(wavena*1e+8,fluxalphana, color='red', label='No Pressure Broadening')
plt.plot(wavena*1e+8,fluxbetana, label='Pressure Broadening')
plt.xlabel('Wavelength ($\AA$)', fontweight='bold')
plt.ylabel('Relative Flux', fontweight='bold')
plt.legend(fancybox = True,loc='lower right', shadow='True')
plt.savefig('Na1.png')
plt.show()


