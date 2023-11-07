#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 12:57:40 2018

ASTR 545: Homework 4

1) Neutral Hydrogen Bound-Free Absorption
2) H− Ion Bound-Free Absorption
3) Plotting the “Total” Cross Section and Discussing the Balmer Decrement
    

@author: oanavesa
"""
#%%
# importing the necessary imports
import numpy as np
import matplotlib.pyplot as plt

#%%


# Defining the temperatures and electron densities for each star
tau = 2/3
TK5V = 4470 # (K)
neK5V = 2.6*10**(14.) # (cm^-3)
TG5V = 5930 # (K)
neG5V = 1.2*10**(14.) # (cm^-3)
TA5V = 8430 # (K)
neA5V = 3.9*10**(14.) # (cm^-3)
TB5V = 13980 # (K)
neB5V = 4.4*10**(13.) # (cm^-3)


# defining the wavelength
wavelength = np.arange(1000,20001,1) # (A)
wavelength_cm=wavelength*1e-8

# defining constants
h = 6.62606885*10**(-27) # ergs*s
c = 30000000000 #cm/s
k = 1.380658e-16 # cm2 g s-2 K-1
Aoprime = 1.0449*10**(-26) # (cm^2 A^-3)
RH = 13.598 # (eV)
RH_cgs = 2.1786397e-11 #(ergs)
lam_n = (h*c)/RH_cgs

#%%

def gaunt_factor(n,lamrange):
    '''Calculates the bound-free Gaunt factor used in the alphaBFHI(wavelength, T) definition.
    Input: n = energy level; lamrange = wavelength range (Å).
    Output: returns g = Gaunt factor '''
    en = (lam_n*n**2/lamrange) - 1
    if 1<en<0:
        top = (121/700)*(1-en**2)
        bottom  = n**(2/3)*(1+en**2)**(2/3)
        g = 1 - (top/bottom)
    elif en>1:
        kn = n*en**(-1/2)
        top2 = (121/700)*(1-(1/en**2))
        bottom2   = kn**(2/3)*(1+(1/en**2))**(2/3)
        g = 1+ (top2/bottom2)
    else:
        g = 0
    return g


#%%

'''Problem 1 part (b) '''

# list of n energy levels 1-6
nlist = np.arange(1,7)
# 7th excitation stage
x7 = RH_cgs*(1-(1/7**2))
# 1st excitation stage
x1 = RH_cgs

# Neutral Hydrogen Bouns-Free Absorption
def alphaBFHI(wavelength, T):
    ''' Calculates the α^bf_HI(λ) values to plot it.
    Inputs: wavelength = wavelength range (Å); T = temperature (K) '''
    flu = []
    for j in range(0,len(wavelength)):
        front = Aoprime*wavelength[j]**3
        #Unsold Contribution
        ktterm = ((k*T)/(2*x1))*(np.exp(-x7/(k*T))-np.exp(-x1/(k*T)))
        problem=0
        # loops through energy levels 1-6
        for i in range(0,len(nlist)): 
            xn = RH_cgs*(1-(1/nlist[i]**2))
            frontterm = (gaunt_factor(nlist[i],wavelength_cm[j])/nlist[i]**3)*np.exp(-xn/(k*T))
            problem += frontterm
        probssss = front*(ktterm + problem)
        flu.append(probssss)
    return flu

alphaK = alphaBFHI(wavelength, TK5V)
alphaG = alphaBFHI(wavelength, TG5V)
alphaA = alphaBFHI(wavelength, TA5V)
alphaB = alphaBFHI(wavelength, TB5V)


# plotting the Neutral Hydrogen Bouns-Free Absorption for Sanity Checks
fig= plt.figure()
plt.subplot(2, 2, 1)
plt.plot(wavelength,alphaK)
plt.xlim([1000,16000])

plt.subplot(2, 2, 2)
plt.plot(wavelength,alphaG)
plt.xlim([1000,16000])

plt.subplot(2, 2, 3)
plt.plot(wavelength,alphaA)
plt.xlim([1000,16000])

plt.subplot(2, 2, 4)
plt.plot(wavelength,alphaB)
plt.xlim([1000,16000])

plt.show()

#%%
def gaunt_factor2(n,lamrange):
    ''' Calculates the bound-free Gaunt factor for each level n.
    Input: n = energy level; lamrange: wavelength range (cm)
    Output: returns g = the bound-free Gaunt factor.
    Used to calculate the individual Gaunt factors for n=1-6 below'''
    for j in range(0,len(lamrange)):
        en = (lam_n*n**2/lamrange[j]) - 1
        if np.where(en <= 1 and en >= 0):
            top = (121/700)*(1-en**2)
            bottom  = n**(2/3)*(1+en**2)**(2/3)
            g = 1 - (top/bottom)
        elif np.where(en>1):
            kn = n*en**(-1/2)
            top2 = (121/700)*(1-(1/en**2))
            bottom2   = kn**(2/3)*(1+(1/en**2))**(2/3)
            g = 1+ (top2/bottom2)
        elif np.where(en <0):
            g = 0
        else:
            g = 0
    return g

# empty list to append the bound-free Gaunt factors
gaunt_values = []
for t in range(0,len(nlist)):
    value = gaunt_factor2(nlist[t],wavelength_cm)
    num = nlist[t]
    gaunt_values.append(value)
    print("n value= " + str(num), "Gaunt value= " + str(value))
#%%
# calculates the values of α^bf_HI(λ) at λ_n for n=2,3,4

def return_alphabfHI_values_for_lam_n(num,numb,T,wavelength2):
    ''' Calculates the α^bf_HI(λ) at λ_n for n=1,2,3,4.
    Inputs: num = index from list gaunt_values where I stored the bound-free Gaunt 
    values calculated in the above for loop using the definition gaunt_factor2;
    numb = 1,2,3,4 = energy levels n, T = temperature (K), and wavelength (Å).
    Output: returns the α^bf_HI(λ) value '''
    front = Aoprime*wavelength2**3
    ktterm = ((k*T)/(2*x1))*(np.exp(-x7/(k*T))-np.exp(-x1/(k*T)))
    xn = RH_cgs*(1-(1/numb**2))
    frontterm = (gaunt_values[num]/numb**3)*np.exp(-xn/(k*T))
    probssss = front*(frontterm)
    return probssss

def Unsold_contribution(T):
    '''Calculates the Unsold Contribution from the bound-free Gaunt factor.
    Input: T = temperature (K). Based on the Boltzmann's cosntant = k,
    the first excitation stage = Ryderg's constant, and the 7th excitation stage 
    calculated at the beginning of the code'''
    ktterm = ((k*T)/(2*x1))*(np.exp(-x7/(k*T))-np.exp(-x1/(k*T)))
    return ktterm


print("Star: " + str(TK5V) + " K ")
TK5V_n2 = return_alphabfHI_values_for_lam_n(1,2,TK5V,1e8*lam_n*2**2)
TK5V_n3 = return_alphabfHI_values_for_lam_n(2,3,TK5V,1e8*lam_n*3**2)
TK5V_n4 = return_alphabfHI_values_for_lam_n(3,4,TK5V,1e8*lam_n*4**2)
TK5V_unsold =  Unsold_contribution(TK5V)
print("The alpha values are: \n" + "n=2 :" + str(TK5V_n2), "\n n=3 :" + str(TK5V_n3), "\n n=4 :" + str(TK5V_n4))
print("The Unsold Contribution is " + str(TK5V_unsold))

print(" \n Star: " + str(TG5V) + " K ")
TG5V_n2 = return_alphabfHI_values_for_lam_n(1,2,TG5V,1e8*lam_n*2**2)
TG5V_n3 = return_alphabfHI_values_for_lam_n(2,3,TG5V,1e8*lam_n*3**2)
TG5V_n4 = return_alphabfHI_values_for_lam_n(3,4,TG5V,1e8*lam_n*4**2)
TG5V_unsold =  Unsold_contribution(TG5V)
print("The alpha values are: \n" + "n=2 :" + str(TG5V_n2), "\n n=3 :" + str(TG5V_n3), "\n n=4 :" + str(TG5V_n4))
print("The Unsold Contribution is " + str(TG5V_unsold))


print(" \n Star: " + str(TA5V) + " K ")
TA5V_n2 = return_alphabfHI_values_for_lam_n(1,2,TA5V,1e8*lam_n*2**2)
TA5V_n3 = return_alphabfHI_values_for_lam_n(2,3,TA5V,1e8*lam_n*3**2)
TA5V_n4 = return_alphabfHI_values_for_lam_n(3,4,TA5V,1e8*lam_n*4**2)
TA5V_unsold =  Unsold_contribution(TA5V)
print("The alpha values are: \n" + "n=2 :" + str(TA5V_n2), "\n n=3 :" + str(TA5V_n3), "\n n=4 :" + str(TA5V_n4))
print("The Unsold Contribution is " + str(TA5V_unsold))

print(" \n Star: " + str(TB5V) + " K ")
TB5V_n2 = return_alphabfHI_values_for_lam_n(1,2,TB5V,1e8*lam_n*2**2)
TB5V_n3 = return_alphabfHI_values_for_lam_n(2,3,TB5V,1e8*lam_n*3**2)
TB5V_n4 = return_alphabfHI_values_for_lam_n(3,4,TB5V,1e8*lam_n*4**2)
TB5V_unsold =  Unsold_contribution(TB5V)
print("The alpha values are: \n" + "n=2 :" + str(TB5V_n2), "\n n=3 :" + str(TB5V_n3), "\n n=4 :" + str(TB5V_n4))
print("The Unsold Contribution is " + str(TB5V_unsold))


print("\n Lambda_n for n=2 is " + str(1e8*lam_n*2**2) + " Å" )
print("Lambda_n for n=3 is " + str(1e8*lam_n*3**2) + " Å" )
print("Lambda_n for n=4 is " + str(1e8*lam_n*4**2) + " Å" )



#%%

''' Problem 2: H- Ion Bound-Free Absorption '''

# new wavelength range defined
wavelength2 = np.arange(1000.,16001.,1.) #(Å)


# a coefficients from Polynomial fits of Wishort 
ai_values = [1.99654,-1.18267*(10**(-5.)),2.64243*(10**(-6.)),-4.40524*(10**(-10.)),3.23992*(10**(-14.)), -1.39568*(10**(-18.)),2.78701*(10**(-23.))]



def HIonAbs(T,ne):
    ''' Calculates the H- Ion Bound-Free Absorption.
    Input: T = temperature (K), ne = electron densities (cm^-3)
    Output: alpp = H- Ion Bound-Free Absorption '''
    alpp = []
    for j in range(0,len(wavelength2)):
        most = 0
        ph=(9.66*(10**(15.)))*((T)**(3.0/2.0))*(np.exp((-1.20964*(10**(-12.)))/(k*T)))
        wave = wavelength2[j]
        for i in range(0,7):
            aiterm = ai_values[i]
            st = aiterm*(wave**(1.*i))
            most += st
        alp = (1*10**(-18.))*most*(ne/ph)h
        alpp.append(alp)
    return alpp

#%%

# calculates the transition wavelength for n=2 - Balmer decrement
m=2
lam_n_full_value = 1e8*((h*c)/RH_cgs)*m**2
print(lam_n_full_value)

# calculates the value of alpha at n=2
def alphaH_HI_value_at_lam_2(T,ne):
    ''' Calculates the alpha value for the H- ion bound-free absorption at n=2
    which is the Balmer decrement.
    Input: T = temperature (K); ne = electron density (cm^-3)
    Output: alpp2 = H- ion bound-free absorption at n=2 '''
    alpp2 = []
    most = 0
    ph=(9.66*10**(15))*((T)**(3.0/2.0))*(np.exp((-1.20964*10**(-12))/(k*T)))
    for i in range(0,7):
        aiterm = ai_values[i]
        st = aiterm*(lam_n_full_value**(1.*i))
        most += st
    alp = 1e-18*most*(ne/ph)
    alpp2.append(alp)
    return alpp2


HIonboundK = HIonAbs(TK5V,neK5V)
HIonboundG = HIonAbs(TG5V,neG5V)
HIonboundA = HIonAbs(TA5V,neA5V)
HIonboundB = HIonAbs(TB5V,neB5V)


print("Values for α^bf_H-/HI(λ) at Balmer decrement, λ_2 = " + str(lam_n_full_value) + " Å: \n" )
lam2_K5V = alphaH_HI_value_at_lam_2(TK5V,neK5V)
lam2_G5V = alphaH_HI_value_at_lam_2(TG5V,neG5V)
lam2_A5V = alphaH_HI_value_at_lam_2(TA5V,neA5V)
lam2_B5V = alphaH_HI_value_at_lam_2(TB5V,neB5V)
print("Star K5V:" + str(lam2_K5V))
print("Star G5V:" + str(lam2_G5V))
print("Star A5V:" + str(lam2_A5V))
print("Star B5V:" + str(lam2_B5V))


#%%

# Plotting the H- Ion Bound-Free Absorption for sanity checks
fig= plt.figure()
plt.subplot(2, 2, 1)
plt.plot(wavelength2,HIonAbs(TK5V,neK5V))

plt.subplot(2, 2, 2)
plt.plot(wavelength2,HIonAbs(TG5V,neG5V))

plt.subplot(2, 2, 3)
plt.plot(wavelength2,HIonAbs(TA5V,neA5V))

plt.subplot(2, 2, 4)
plt.plot(wavelength2,HIonAbs(TB5V,neB5V))

plt.show()

#%%
'''Problem 3 part a: n_H/n_HI ratio'''

# phi function
def phi_T(T):
    ''' Calculates the Phi function
    Input: T = temperature (K)
    Output: phi '''
    phi = (9.66*10**(15))*T**(3/2)*np.exp(-(1.20964*10**(-12))/(k*T))
    return phi

# calculates ratio
def alp_ratio(ne,T):
    ''' Calculates the n_H/n_HI ratio = ne/phi
    Input: ne = electron density (cm^-3); T = temperature (K)
    Output: ratio = ne/phi'''
    ratio = ne/phi_T(T)
    return ratio

# appends the ratios to a list
ro = [alp_ratio(neK5V,TK5V),alp_ratio(neG5V,TG5V),alp_ratio(neA5V,TA5V),alp_ratio(neB5V,TB5V)]
label = ['K5V','G5V','A5V','B5V']


#plots the ratio
plt.figure()
plt.plot(ro)
plt.title('$n_{H}/n_{HI}$', fontsize=14, weight='bold')
#plt.savefig('ratio2,png')
plt.show()


#%%
'''Problem 3 part b: total absorption'''

# adds up both alphas to get the total absorption
def tota(alphaneutral,alphaion):
    total2 = []
    for y in range(0,len(HIonboundK)):
        total = alphaneutral[y] + alphaion[y]
        total2.append(total)
    return total2    

# plots the total alpha for sanity checks
fig= plt.figure()
plt.subplot(2, 2, 1)
plt.plot(wavelength2,tota(alphaK,HIonboundK))
plt.xlim([1000,16000])
 
plt.subplot(2, 2, 2)
plt.plot(wavelength2,tota(alphaG,HIonboundG))
plt.xlim([1000,16000])

plt.subplot(2, 2, 3)
plt.plot(wavelength2,tota(alphaA,HIonboundA))
plt.xlim([1000,16000])

plt.subplot(2, 2, 4)
plt.plot(wavelength2,tota(alphaB,HIonboundB))
plt.xlim([1000,16000])

plt.show()

#%%


fig= plt.figure()
fig.subplots_adjust(hspace=0.4)
plt.suptitle('Cross Sections', weight='bold', fontsize='14')
ax=plt.subplot(2, 2, 1)
plt.plot(wavelength,alphaK, '--')
plt.plot(wavelength2,HIonAbs(TK5V,neK5V),':',color='red')
plt.plot(wavelength2,tota(alphaK,HIonboundK), 'k')
plt.xlim([1000,16000])
plt.title('Star K5V', weight='bold')
plt.ylim(bottom=0)


plt.subplot(2, 2, 2)
plt.plot(wavelength,alphaG, '--')
plt.title('Star G5V', weight='bold')
plt.plot(wavelength2,HIonAbs(TG5V,neG5V),':',color='red')
plt.plot(wavelength2,tota(alphaG,HIonboundG),'k')
plt.xlim([1000,16000])
plt.ylim(bottom=0)


plt.subplot(2, 2, 3)
plt.title('Star A5V', weight='bold')
plt.plot(wavelength,alphaA, '--')
plt.plot(wavelength2,HIonAbs(TA5V,neA5V),':',color='red')
plt.plot(wavelength2,tota(alphaA,HIonboundA), 'k')
plt.xlim([1000,16000])
plt.ylim(bottom=0)


plt.subplot(2, 2, 4)
plt.title('Star B5V', weight='bold')
plt.plot(wavelength,alphaB, '--', label = 'Neutral Hydrogren')
plt.plot(wavelength2,HIonAbs(TB5V,neB5V),':',color='red',label = '$H^-$ Bound-Free Hydrogren')
plt.plot(wavelength2,tota(alphaB,HIonboundB),'k',label = 'Total Bound-Free')
plt.xlim([1000,16000])
plt.legend(loc='upper center', bbox_to_anchor=(-.2, -0.44),
          fancybox=True, shadow=True, ncol=3)
plt.ylim(bottom=0)

fig.text(0.06, 0.5, 'Absorption $(cm^{2}/HI)$ ', ha='center', va='center', rotation='vertical', weight='bold', fontsize=12)
fig.text(0.5, 0.04, 'Wavelength $(\AA)$', ha='center', va='center', weight='bold', fontsize='12')
plt.savefig('totl.png')


plt.show()
