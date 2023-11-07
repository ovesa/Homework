#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 7 14:26:39 2018

ASTR 545: Homework 5B. 
1) Ionization Fractions of Calcium
2) Excitation and Ionization– the “Balmer Hydrogen”.
3) Equilibrium Gas Physics.
4) Equilibrium Gas Physics as a Function of Temperature

@author: oanavesa
"""

# import necessary packages
import numpy as np
from scipy import optimize
from astropy.io import ascii
import matplotlib.pyplot as plt 


# ionization stages for Calcium:
# n1,20 for neutral, n2,20 for singly ionized, and n3,20 for doubly ionized
# total calcium number density is n20 = n1,20 + n2,20 + n3,20
# text files
path = '/Users/oanavesa/Desktop/GradSchool/FirstYear/StellarSpectroscopy/HW/HW5/PartB/'
U_pot = 'PartFuncs.txt'
data = ascii.read(path+U_pot,guess=False, format='basic')

# Chi values
chi_values = [13.5984340*1.6022e-12,24.5873876*1.6022e-12,54.4177630*1.6022e-12,6.113158*1.6022e-12,11.87172*1.6022e-12] 
# hard coded partition functions
u13=10**0.181
u23=10**0.391
u33 =10**0.

# defining constants
T_G3V = 6000. # (K)
ne_G3V = 2*10**(13.) # (cm^-3)
T_G3Ia = 6000. # (K)
ne_G3Ia = 8*10**(11.) #(cm^-3)
C_phi = 4.83*10**(15.)
 R = 13.598*1.60218e-12 #ergs
f_11 = 0.15
T_problem2 = 10000 #K
U11_prob2 =  0.302
P_prob2 = 100e0 #dynes cm^-3
A_H = 1.00794e0
A_He = 4.002602e0
x_H = 0.7
x_He = 0.3


T_prob3 = 6000 #K
P_g = 100e0 #dynes cm−2 
X = 0.71e0 # mass fraction of H
Y = 0.27e0 # mass fraction of He
Z = 0.02e0 # mass fraction of Ca

A_Ca = 40.078e0 # atomic weights
m_h =1.6737236*10**(-24) # mass of hydrogen
m_he = 6.6464764*10**(-24) # mass of helium
me = 9.10938*10**(-28) # mass of electron
m_ca = 6.6551079*10**(-23) # mass of calcium
ma = 1.66054e-24

# definitions start here and are interspered throughout
def phi_jk(C_phi,T,k,U_jk1,U_jk,XI):
   first = C_phi*T**(3/2)
   pot = ((U_jk1)/U_jk)
   expot = np.exp(-XI/(k*T))
   total = first*pot*expot
   return total

# ionization stages

def f13(Y13,Y23):
    f13 = 1 + Y13 +Y13*Y23
    return 1/f13

def f23(f13,Y13):
    f23 = f13*Y13
    return f23

def f33(f23,Y23):
    f33 = f23*Y23
    return f33


def f11(Y11):
    f11 = 1+Y11
    return 1/f11


def f21(f11,Y11):
    f21 = f11*Y11
    return f21

def f12(Y12,Y22):
    f12 = 1+Y12+Y22*Y12
    return 1/f12

def f22(f12,Y12):
    f22 = f12*Y12
    return f22

def f32(f22,Y22):
    f32 = f22*Y22
    return f32


def phi_jk2(C_phi,T,k,U_jk1,U_jk,XI):
   first = C_phi*T**(3/2)
   pot = ((U_jk1)/U_jk)
   expot = np.exp(-XI/(k*T))
   total = first*pot*expot
   return total
             


def n211_n1(n,U_11,R,k,T,f11):
    g = 2*n**2
    U = U_11
    expot = np.exp(-(R*(1-(1/n**2)))/(k*T))
    total = (g/U)*expot*f11
    return total
    

def solve_for_nuclear_particle_density(P,k,T):
    n = (P)/(k*T)
    return n

def abundance_fraction(x_1,A_1,x_2,A_2):
    alpha_1 = x_1/A_1
    sum_alpha = alpha_1 + (x_2/A_2) 
    abundance_1 = alpha_1/sum_alpha
    return abundance_1

def solve_for_n_k(alpha_k,n_N):
    n_k = alpha_k*n_N
    return n_k

def transcential_eq(ne):
    '''Computes the transcential equation to solve for ne.
    Input: ne = electron number density and T = temperature
    Output: Used in conjunction with the brentq method to sovle for ne'''
    Y11 = (ne**-1)*(phi_jk2(C_phi,T_prob3,k,10**0,10**.301,chi_values[0]))
    f11val = f11(Y11)
    Y12 = (ne**-1)*(phi_jk2(C_phi,T_prob3,k,10**.301,10**0,chi_values[1]))
    Y22 = (ne**-1)*(phi_jk2(C_phi,T_prob3,k,10**0,10**.301,chi_values[2]))
    f12_val = f12(Y12,Y22)
    f22_val = f22(f12_val,Y12)
    Y13 = (ne**-1)*(phi_jk2(C_phi,T_prob3,k,10**.391,10**.181,chi_values[3]))
    Y23 = (ne**-1)*(phi_jk2(C_phi,T_prob3,k,10**0,10**.391,chi_values[4]))
    f13_val = f13(Y13,Y23)
    f23_val = f23(f13_val,Y13)


    first_term = solve_for_nuclear_particle_density(P_g,k,T_prob3)-ne
    abundances_H = abundance_fraction2(X,A_H,Y,A_He,Z,A_Ca)*( (0)*(f11(Y11))+(1)*(f21(f11val,Y11)) )
    abundances_He = abundance_fraction2(Y,A_He,X,A_H,Z,A_Ca)*( (0)*(f12(Y12,Y22))+(1)*(f22(f12_val,Y12))+(2)*(f32(f22_val,Y22)) )
    abundances_Ca =  abundance_fraction2(Z,A_Ca,X,A_H,Y,A_He)*( (0)*( f13(Y13,Y23))+(1)*(f23(f13_val,Y13))+(2)*(f33(f23_val,Y23)))
    second_term = abundances_H+abundances_He+abundances_Ca
    third_term = first_term*second_term
    fourth_term = ne
    final = third_term-fourth_term
    return final
    
def abundance_fraction2(x_1,A_1,x_2,A_2,x_3,A_3):
    alpha_1 = x_1/A_1
    sum_alpha = alpha_1 + (x_2/A_2) +(x_3/A_3)
    abundance_1 = alpha_1/sum_alpha
    return abundance_1

'''Problem 1': Computing the Calcium ionization fractions'''
G3Ia_phi_13 = phi_jk(C_phi,T_G3Ia,k,u23,u13,chi_values[3])
G3Ia_Saha_13 = G3Ia_phi_13*(1/ne_G3Ia)

G3Ia_phi_23 = phi_jk(C_phi,T_G3Ia,k,u33,u23,chi_values[4])
G3Ia_Saha_23 = G3Ia_phi_23*(1/ne_G3Ia)

G3Ia_f13 = f13(G3Ia_Saha_13,G3Ia_Saha_23)
G3Ia_f23 = f23(G3Ia_f13,G3Ia_Saha_13)
G3Ia_f33 = f33(G3Ia_f23,G3Ia_Saha_23)

G3V_phi_13 = phi_jk(C_phi,T_G3V,k,u23,u13,chi_values[3])
G3V_Saha_13 = G3V_phi_13*(1/ne_G3V)

G3V_phi_23 = phi_jk(C_phi,T_G3V,k,u33,u23,chi_values[4])
G3V_Saha_23 = G3V_phi_23*(1/ne_G3V)

G3V_f13 = f13(G3V_Saha_13,G3V_Saha_23)
G3V_f23 = f23(G3V_f13,G3V_Saha_13)
G3V_f33  = f33(G3V_f23,G3V_Saha_23)

n211_n1_value = n211_n1(2,U11_prob2,R,k,T_problem2,f_11)
#part a) nuclear particle number density
nuclear_particle_density = solve_for_nuclear_particle_density(P_prob2,k,T_problem2)
# part b) abundance fractions
abundance_h = abundance_fraction(x_H,A_H,x_He,A_He)
# part c) number densitiesn_va
n_1 = solve_for_n_k(abundance_h,nuclear_particle_density)
solve_n_211 = n211_n1_value*n_1

#%%
'''Problem 3'''
# solving for ne using the transcential_eq above
root = optimize.brentq(transcential_eq,6e11,7e15,xtol= 10**-5)


# For table 1
abundance_H = abundance_fraction2(X,A_H,Y,A_He,Z,A_Ca)
abundance_He = abundance_fraction2(Y,A_He,X,A_H,Z,A_Ca)
abundance_Ca = abundance_fraction2(Z,A_Ca,X,A_H,Y,A_He)
n_value = solve_for_nuclear_particle_density(P_g,k,T_prob3)
electron_density_val = root
n_total = n_value-electron_density_val


# For Table 2
def electron_pressure(k,T,ne):n
    e_pressure = (ne)*(k)*(T)
    return e_pressure

def mu_n(x_1,A_1,x_2,A_2,x_3,A_3):
    alpha_1 = x_1/A_1
    alpha_2 = x_2/A_2
    alpha_3 = x_3/A_3
    mu = alpha_1+alpha_2+alpha_3
    mu_n = (mu)**(-1)
    return mu_n

# Saha Equation values and ionization fractions
Y11 = (electron_density_val**-1)*(phi_jk2(C_phi,T_prob3,k,10**0,10**.301,chi_values[0]))
f11val = f11(Y11)
f21_val= f21(f11val,Y11)
Y12 = (electron_density_val**-1)*(phi_jk2(C_phi,T_prob3,k,10**.301,10**0,chi_values[1]))
Y22 = (electron_density_val**-1)*(phi_jk2(C_phi,T_prob3,k,10**0,10**.301,chi_values[2]))
f12_val = f12(Y12,Y22)
f22_val = f22(f12_val,Y12)
Y13 = (electron_density_val**-1)*(phi_jk2(C_phi,T_prob3,k,10**.391,10**.181,chi_values[3]))
Y23 = (electron_density_val**-1)*(phi_jk2(C_phi,T_prob3,k,10**0,10**.391,chi_values[4]))
f13_val = f13(Y13,Y23)
f23_val = f23(f13_val,Y13)
f32_value = f32(f22_val,Y22)
f33_val = f33(f23_val,Y23)
f22_val = f22(f12_val,Y12)

# definitions
def mu_e(x_1,A_1,x_2,A_2,x_3,A_3,f11,f21,f12,f22,f32,f13,f23,f33):
    alpha_1 = x_1/A_1
    alpha_2 = x_2/A_2
    alpha_3 = x_3/A_3
    mu_H = alpha_1*(0)*f11 + alpha_1*1*f21
    mu_He = alpha_2*(0)*f12 + alpha_2*1*f22+alpha_2*2*f32
    mu_Ca = alpha_3*(0)*f13 + alpha_3*(1)*f23+alpha_3*(2)*f33
    total = mu_H +mu_He +mu_Ca
    total2 = total**(-1)
    return total2

def total_mu(mu_n_val,mu_e_val):
    inverse_mun = 1/mu_n_val
    inverse_mue = 1/mu_e_val
    inverse_mu_tot = inverse_mun + inverse_mue
    mu_tot = 1/inverse_mu_tot
    return mu_tot

def electron_number_density(n_1,n_2,n_3,f11,f21,f12,f22,f32,f13,f23,f33):
    electron_num_dens_H = n_1*(0)*f11 + n_1*1*f21
    electron_num_dens_He = n_2*(0)*f12 + n_2*1*f22+n_2*2*f32
    electron_num_dens_Ca= n_3*(0)*f13 + n_3*1*f23+n_3*2*f33
    total = electron_num_dens_H+electron_num_dens_He +electron_num_dens_Ca
    return total

def electron_density(me,n_h,n_he,n_Ca):
    e_dens = me*(1*n_h  + 2*n_he + 3*n_Ca)
    return e_dens

def solve_for_n_k(alpha_k,n_N):
    n_k = alpha_k*n_N
    return n_k

def nuclear_mass_density(m_H,m_He,m_Ca,n_H,n_He,n_Ca):
    p_n = m_H*n_H + m_He*n_He + m_Ca*n_Ca
    return p_n

def fraction_of_electrons(j,n_jk,ne):
    f_ne_j_k = (j-1)*(n_jk/ne)
    return f_ne_j_k

def solve_for_number_densities_of_ionization_stages(f_jk,n_k):
    n_jk = f_jk*n_k
    return n_jk



# calculating pressure
electron_pressure_value = electron_pressure(k,T_prob3,electron_density_val)
ratio_pe_pg = electron_pressure_value/P_g

# calculating molecular weight
molecular_n_value = mu_n(X,A_H,Y,A_He,Z,A_Ca)
mu_e_value = mu_e(X,A_H,Y,A_He,Z,A_Ca,f11val,f21_val,f12_val,f22_val,f32_value,f13_val,f23_val,f33_val)
total_mu_value = total_mu(molecular_n_value,mu_e_value)
# solving for n_k
n_H = solve_for_n_k(abundance_H,n_total)
n_He = solve_for_n_k(abundance_He,n_total)
n_Ca = solve_for_n_k(abundance_Ca,n_total)
# solving for densities
electron_density_value = me*electron_density_val
nuclear_particle_mass_density =  nuclear_mass_density(m_h,m_he,m_ca,n_H,n_He,n_Ca)
total_density = (1 + (electron_density_value/nuclear_particle_mass_density))*nuclear_particle_mass_density



#table 3
n_H = solve_for_n_k(abundance_H,n_total) # H density
n_He = solve_for_n_k(abundance_He,n_total) #He density
n_Ca = solve_for_n_k(abundance_Ca,n_total)#Ca density
n11 = solve_for_number_densities_of_ionization_stages(f11val,n_H) #HI density
n21 = solve_for_number_densities_of_ionization_stages(f21_val,n_H) #HII density
n12 = solve_for_number_densities_of_ionization_stages(f12_val,n_He) #HeI density
n22 = solve_for_number_densities_of_ionization_stages(f22_val,n_He) #HeII density
n32  = solve_for_number_densities_of_ionization_stages(f32_value,n_He) #HeIII density
n13 = solve_for_number_densities_of_ionization_stages(f13_val,n_Ca) #CaI density
n23 = solve_for_number_densities_of_ionization_stages(f23_val,n_Ca) #CaII density
n33 = solve_for_number_densities_of_ionization_stages(f33_val,n_Ca) #CaIII density

# phi(T)
phi11 = (phi_jk2(C_phi,T_prob3,k,10**0,10**.301,chi_values[0])) #HI
phi12 = (phi_jk2(C_phi,T_prob3,k,10**.301,10**0,chi_values[1])) # HeI
phi22 = (phi_jk2(C_phi,T_prob3,k,10**0,10**.301,chi_values[2])) #HeII
phi13 = (phi_jk2(C_phi,T_prob3,k,10**.391,10**.181,chi_values[3])) #CaI
phi23 = (phi_jk2(C_phi,T_prob3,k,10**0,10**.391,chi_values[4])) #CaII

# partition function
logu11 = 10**0/10**.301
logu12 = 10**.301/10**0
logu22 = 10**0/10**.301
logu13 =10**.391/10**.181
logu23 = 10**0/10**.391



#fraction Contributions to Free Electrons
f_ne_11 = fraction_of_electrons(1,n11,electron_density_val) #HI
f_ne_21 = fraction_of_electrons(2,n21,electron_density_val) #HII
f_ne_12 = fraction_of_electrons(1,n12,electron_density_val) #HeI
f_ne_22 = fraction_of_electrons(2,n22,electron_density_val) #HeII
f_ne_32 = fraction_of_electrons(3,n32,electron_density_val) #HeIII
f_ne_13 = fraction_of_electrons(1,n13,electron_density_val) #CaI
f_ne_23 = fraction_of_electrons(2,n23,electron_density_val) #CaII
f_ne_33 = fraction_of_electrons(3,n33,electron_density_val) #CaIII




#%%
''' Problem 4'''

path = '/Users/oanavesa/Desktop/GradSchool/FirstYear/StellarSpectroscopy/HW/HW5/PartB/'
U_pot = 'PartFuncs.txt'
same = 'same.txt'
data = ascii.read(path+U_pot,guess=False, format='basic')
upot=np.loadtxt(path+U_pot,unpack=True) #partition function
upotvalues = np.loadtxt(path+same,unpack=True) #partition function
theta = [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,9999]
# temperature list obtained from theta values from text file above
Templist = [5040/0.2,5040/0.4,5040/0.6,5040/0.8,5040/1.0,5040/1.2,5040/1.4,5040/1.6,5040/1.8,5040/20,5040/9999] 

#Chi values
chi_values = [13.5984340*1.6022e-12,24.5873876*1.6022e-12,54.4177630*1.6022e-12,6.113158*1.6022e-12,11.87172*1.6022e-12] 

# temperature range of 1500 - 25000 K.
list_of_temp = np.arange(1500,25000,250)


def general_phi(j,kg,T):
    ''' The generalized Phi(T) value for a given ionization stage (j), species
    (kg), and temperature (T).
    Output: the general Phi(T) value'''
    # gives you index for column for which theta matches
    minlist = []
    for i in range(0,len(Templist)):
        pick = Templist[i]
        difference = np.abs(pick-T)
        minlist.append(difference)
        minvalue = minlist.index(min(minlist))
        
    
    first = C_phi*T**(3/2)
    
    
    if kg==1:
        upot_val1 = upotvalues[minvalue][j]
        upot_val2 = upotvalues[minvalue][j+1]
        potential_value = (10**upot_val2)/(10**upot_val1)
        chival = chi_values[j-1]
        finalval = first*potential_value*np.exp(-chival/(k*T))
        
    elif kg==2:
        upot_val1 = upotvalues[minvalue][j+2]
        upot_val2 = upotvalues[minvalue][j+3]
        potential_value = (10**upot_val2)/(10**upot_val1)
        chival = chi_values[j]
        finalval = first*potential_value*np.exp(-chival/(k*T))


    elif kg==3:
        upot_val1 = upotvalues[minvalue][j+5]
        upot_val2 = upotvalues[minvalue][j+6]
        potential_value = (10**upot_val2)/(10**upot_val1)
        chival = chi_values[j+2]
        finalval = first*potential_value*np.exp(-chival/(k*T))
    
    return finalval




def transcential_eq_prob4(ne,T):
    '''Computes the transcential equation to solve for ne.
    Input: ne = electron number density and T = temperature
    Output: Used in conjunction with the brentq method to sovle for ne'''
    Y11 = (ne**-1)*general_phi(1,1,T)
    f11val = f11(Y11)
    Y12 = (ne**-1)*general_phi(1,2,T)
    Y22 = (ne**-1)*general_phi(2,2,T)
    f12_val = f12(Y12,Y22)
    f22_val = f22(f12_val,Y12)
    Y13 = (ne**-1)*general_phi(1,3,T)
    Y23 = (ne**-1)*general_phi(2,3,T)
    f13_val = f13(Y13,Y23)
    f23_val = f23(f13_val,Y13)
    

    nt = solve_for_nuclear_particle_density(P_g,k,T)
    first_term = (nt-ne)
    abundances_H = abundance_fraction2(X,A_H,Y,A_He,Z,A_Ca)*( (0)*(f11(Y11))+(1)*(f21(f11val,Y11)) )
    abundances_He = abundance_fraction2(Y,A_He,X,A_H,Z,A_Ca)*( (0)*(f12(Y12,Y22))+(1)*(f22(f12_val,Y12))+(2)*(f32(f22_val,Y22)) )
    abundances_Ca =  abundance_fraction2(Z,A_Ca,X,A_H,Y,A_He)*( (0)*( f13(Y13,Y23))+(1)*(f23(f13_val,Y13))+(2)*(f33(f23_val,Y23)))
    second_term = abundances_H+abundances_He+abundances_Ca
    third_term = first_term*second_term
    fourth_term = ne
    final = third_term-fourth_term
    return final

# empty lists for all of the variables computed for Problem 4
ne_list = []
n_list=[]
n_total_list = []
rho_e_list = []
rho_list = []
rho_N_list=[]
mu_e_list = []
mu_list = []
mu_N_list = []
pressure_ratio_list = []
f11_hydrogen_list =[]
f21_hydrogen_list =[]
n_H_list = []
n11_hydrogen_list =[]
n21_hydrogen_list =[]

f12_helium_list =[]
f22_helium_list =[]
f32_helium_list =[]
n_He_list = []
n12_helium_list =[]
n22_helium_list =[]
n32_helium_list = []

f13_calcium_list =[]
f23_calcium_list =[]
f33_calcium_list =[]
n_Ca_list = []
n13_calcium_list =[]
n23_calcium_list =[]
n33_calcium_list =[]

f_ne_n11_list = []
f_ne_n21_list = []
f_ne_n12_list = []
f_ne_n22_list = []
f_ne_n32_list = []
f_ne_n13_list = []
f_ne_n23_list = []
f_ne_n33_list = []

# for loop that loops through temperature range to compute
# the equalibrium conditions of the gas structure
for i in range(0,len(list_of_temp)):
    call = list_of_temp[i]
    # root = ne
    root = optimize.brentq(transcential_eq_prob4, 1e5,4e22, xtol = 10**(-5), args=(call))
    ne_list.append(root)
    #solving for nuclear particle density
    n_values =  solve_for_nuclear_particle_density(P_g,k,call)
    n_list.append(n_values)
    # total particle density
    ntotal = n_values-root
    n_total_list.append(ntotal)
    
    # solving for rho
    # density equiliibrium values
    electron_density_value2 = me*root
    rho_e_list.append(electron_density_value2)
    
    # solving for abundances
    abundance_H = abundance_fraction2(X,A_H,Y,A_He,Z,A_Ca)
    abundance_He = abundance_fraction2(Y,A_He,X,A_H,Z,A_Ca)
    abundance_Ca = abundance_fraction2(Z,A_Ca,X,A_H,Y,A_He)
    
    # solving for n_k
    n_H = solve_for_n_k(abundance_H,n_values)
    n_He = solve_for_n_k(abundance_He,n_values)
    n_Ca = solve_for_n_k(abundance_Ca,n_values)
    nuclear_particle_mass_density2 =  nuclear_mass_density(m_h,m_he,m_ca,n_H,n_He,n_Ca)
    rho_N_list.append(nuclear_particle_mass_density2)
    total_density2 = (1 + (electron_density_value2/nuclear_particle_mass_density2))*nuclear_particle_mass_density2
    rho_list.append(total_density2)
    
    # Saha equations and ionization fractions
    Y11 = (root**-1)*general_phi(1,1,call)
    f11val = f11(Y11)
    Y12 = (root**-1)*general_phi(1,2,call)
    Y22 = (root**-1)*general_phi(2,2,call)
    f12_val = f12(Y12,Y22)
    f22_val = f22(f12_val,Y12)
    Y13 = (root**-1)*general_phi(1,3,call)
    Y23 = (root**-1)*general_phi(2,3,call)
    f13_val = f13(Y13,Y23)
    f23_val = f23(f13_val,Y13)
    f33_val = f33(f23_val,Y23)
    f32_value = f32(f22_val,Y22)
    f21_val = f21(f11val,Y11)
    
    # solving for molecular weights
    mu_e_value2 = mu_e(X,A_H,Y,A_He,Z,A_Ca,f11val,f21_val,f12_val,f22_val,f32_value,f13_val,f23_val,f33_val)
    mu_e_list.append(mu_e_value2)
    molecular_n_value2 = mu_n(X,A_H,Y,A_He,Z,A_Ca)
    mu_N_list.append(molecular_n_value2)
    total_mu_value2 = total_mu(molecular_n_value2,mu_e_value2)
    mu_list.append(total_mu_value2)
    
    # solving for pressure
    electron_pressure_value = electron_pressure(k,call,root)
    pressure_ratio_list.append(electron_pressure_value/P_g)

    # Ionization fractions for hydrogen
    f11_hydrogen_list.append(f11val)
    f21_hydrogen_list.append(f21_val)
    n_H_list.append(n_H)
    
    # atom/ion density for hydrogen
    n11_v = solve_for_number_densities_of_ionization_stages(f11val,n_H) #HI density
    n21_v = solve_for_number_densities_of_ionization_stages(f21_val,n_H) #HII density

    n11_hydrogen_list.append(n11_v)
    n21_hydrogen_list.append(n21_v)
    
    # Ionization fractions for helium
    f12_helium_list.append(f12_val)
    f22_helium_list.append(f22_val)
    f32_helium_list.append(f32_value)
    
    # atom/ion density for helium
    n12_v = solve_for_number_densities_of_ionization_stages(f12_val,n_He) #HeI density
    n22_v = solve_for_number_densities_of_ionization_stages(f22_val,n_He) #HeII density
    n32_v = solve_for_number_densities_of_ionization_stages(f32_value,n_He) #HeIII density

    n12_helium_list.append(n12_v)
    n22_helium_list.append(n22_v)
    n32_helium_list.append(n32_v)
    
    n_He_list.append(n_He)
    
    # Ionization fractions for calcium
    f13_calcium_list.append(f13_val)
    f23_calcium_list.append(f23_val)
    f33_calcium_list.append(f33_val)
    
    n_Ca_list.append(n_Ca)
    
    # atom/ion density for calcium
    n13_v = solve_for_number_densities_of_ionization_stages(f13_val,n_Ca) #HeI density
    n23_v = solve_for_number_densities_of_ionization_stages(f23_val,n_Ca) #HeII density
    n33_v = solve_for_number_densities_of_ionization_stages(f33_val,n_Ca) #HeIII density

    n13_calcium_list.append(n13_v)
    n23_calcium_list.append(n23_v)
    n33_calcium_list.append(n33_v)
    
    # fraction contributions to free electrons
    f_ne_11v = fraction_of_electrons(1,n11_v,root) #HI
    f_ne_21v = fraction_of_electrons(2,n21_v,root) #HII
    f_ne_12v = fraction_of_electrons(1,n12_v,root) #HeI
    f_ne_22v = fraction_of_electrons(2,n22_v,root) #HeII
    f_ne_32v = fraction_of_electrons(3,n32_v,root) #HeIII
    f_ne_13v = fraction_of_electrons(1,n13_v,root) #CaI
    f_ne_23v = fraction_of_electrons(2,n23_v,root) #CaII
    f_ne_33v = fraction_of_electrons(3,n33_v,root) #CaIII
    f_ne_n11_list.append(f_ne_11v)
    f_ne_n21_list.append(f_ne_21v)
    f_ne_n12_list.append(f_ne_12v)
    f_ne_n22_list.append(f_ne_22v)
    f_ne_n32_list.append(f_ne_32v)
    f_ne_n13_list.append(f_ne_13v)
    f_ne_n23_list.append(f_ne_23v)
    f_ne_n33_list.append(f_ne_33v)

    
    
'''Plotting: Plot 1 '''
plt.figure()
plt.title('Number Densities',size=14,fontweight='bold')
plt.plot(list_of_temp,np.log10(ne_list), label='Electron', color='g')
plt.plot(list_of_temp,np.log10(n_list),label='Total Particle', color='k')
plt.plot(list_of_temp,np.log10(n_total_list), label='Nuclear', color='r')
plt.ylabel('log($n_{tot}$,$n_N$,$n_e$) [$cm^{-3}]$',size=12)
plt.xlabel('T (K)',size=12)
plt.ylim([10,15])
plt.legend(loc='lower right', fancybox=True, shadow=True, ncol=1)
plt.savefig('Number Densities.png')
plt.show()
#
plt.figure()
plt.title('Mass Densities',size=14,fontweight='bold')
plt.plot(list_of_temp,np.log10(rho_e_list),label='Electron', color='g')
plt.plot(list_of_temp,np.log10(rho_N_list),label='Total Particle', color='k')
plt.plot(list_of_temp,np.log10(rho_list),label='Nuclear', color='r')
plt.ylabel('log($ρ_{tot}$,$ρ_N$,$ρ_e$) [$cm^{-3}]$',size=12)
plt.xlabel('T (K)',size=12)
plt.legend(loc='lower right', fancybox=True, shadow=True, ncol=1)
plt.savefig('Mass Densities.png')
plt.show()
#
plt.figure()
plt.title('Mean Molecular Weights',size=14,fontweight='bold')
plt.plot(list_of_temp,np.log10(mu_e_list),label='Electron', color='g')
plt.plot(list_of_temp,np.log10(mu_N_list),label='Nuclear', color='r')
plt.plot(list_of_temp,np.log10(mu_list),label='Total Particle', color='k')
plt.ylim([-.5,4])
plt.ylabel('log($µ_{tot}$,$µ_N$,$µ_e$) [$cm^{-3}]$',size=12)
plt.xlabel('T (K)',size=12)
plt.legend(loc='upper right', fancybox=True, shadow=True, ncol=1)
plt.savefig('Mean Molecular Weights.png')
plt.show()

plt.figure()
plt.title('Electron/Gas Pressure Ratio',size=14,fontweight='bold')
plt.plot(list_of_temp,pressure_ratio_list, color='g')
plt.axhline(y=0.50,linestyle=':',color='k')
plt.xlabel('T (K)',size=12)
plt.ylabel('$P_e/P_g$',size=12)
plt.ylim(top=1)
plt.savefig('ElectronGas Pressure Ratio.png')
plt.show()


'''Plotting: Plot 2 '''

plt.figure()
plt.title('Hydrogen Ionization Fractions',size=14,fontweight='bold')
plt.xlabel('T (K)',size=12)
plt.ylabel('f($H^0$,f($H^+$))',size=12)
plt.plot(list_of_temp,f11_hydrogen_list,label='$H^0$', color='r')
plt.plot(list_of_temp,f21_hydrogen_list,label='$H^+$', color='g')
plt.legend(loc='middle right', fancybox=True, shadow=True, ncol=1)
plt.savefig('Hydrogen Ionization Fractions.png')
plt.show()


plt.figure()
plt.title('Hydrogen Number Densities',size=14,fontweight='bold')
plt.xlabel('T (K)',size=12)
plt.ylabel('log(n) [$cm^{-3}$]',size=12)
plt.plot(list_of_temp,np.log10(n_H_list), label='Total Number Density', color='b')
plt.plot(list_of_temp,np.log10(n11_hydrogen_list), label='$H^0$ Number Density', color='r')
plt.plot(list_of_temp,np.log10(n21_hydrogen_list),label='$H^+$ Number Density', color='g')
plt.legend(loc='middle right', fancybox=True, shadow=True, ncol=1)
plt.savefig('Hydrogen Number Densities.png')
plt.ylim([9,15])
plt.show()

plt.figure()
plt.title('Helium Ionization Fractions',size=14,fontweight='bold')
plt.xlabel('T (K)',size=12)
plt.ylabel('f($He^0$),f($He^+$,f($He^{++}$))',size=12)
plt.plot(list_of_temp,f12_helium_list,label='$H^0$', color='k')
plt.plot(list_of_temp,f22_helium_list,label='$H^+$', color='r')
plt.plot(list_of_temp,f32_helium_list,label='$H^{++}$', color='g')
plt.legend(loc='middle right', fancybox=True, shadow=True, ncol=1)
plt.savefig('Helium Ionization Fractions.png')
plt.show()

plt.figure()
plt.title('Helium Number Densities',size=14,fontweight='bold')
plt.xlabel('T (K)',size=12)
plt.ylabel('log(n) [$cm^{-3}$]',size=12)
plt.plot(list_of_temp,np.log10(n_He_list), label='Total Number Density', color='b')
plt.plot(list_of_temp,np.log10(n12_helium_list), label='$H^0$ Number Density', color='k')
plt.plot(list_of_temp,np.log10(n22_helium_list), label='$H^+$ Number Density', color='r')
plt.plot(list_of_temp,np.log10(n32_helium_list), label='$H^{++}$ Number Density', color='g')
plt.ylim([9,14])
plt.legend(loc='middle left', fancybox=True, shadow=True, ncol=1)
plt.savefig('Helium Number Densities.png')
plt.show()

plt.figure()
plt.title('Calcium Ionization Fractions',size=14,fontweight='bold')
plt.xlabel('T (K)',size=12)
plt.ylabel('f($Ca^0$),f($Ca^+$,f($Ca^{++}$))',size=12)
plt.plot(list_of_temp,f13_calcium_list,label='$Ca^0$', color='k')
plt.plot(list_of_temp,f23_calcium_list,label='$Ca^+$', color='r')
plt.plot(list_of_temp,f33_calcium_list,label='$Ca^{++}$', color='g')
plt.legend(loc='middle left', fancybox=True, shadow=True, ncol=1)
plt.savefig('Calcium Ionization Fractions.png')
plt.show()

plt.figure()
plt.title('Calcium Number Densities',size=14,fontweight='bold')
plt.xlabel('T (K)',size=12)
plt.ylabel('log(n) [$cm^{-3}$]',size=12)
plt.plot(list_of_temp,np.log10(n_Ca_list), label='Total Number Density', color='b')
plt.plot(list_of_temp,np.log10(n13_calcium_list), label='$Ca^0$ Number Density', color='k')
plt.plot(list_of_temp,np.log10(n23_calcium_list), label='$Ca^+$ Number Density', color='r')
plt.plot(list_of_temp,np.log10(n33_calcium_list), label='$Ca^{++}$ Number Density', color='g')
plt.ylim([7,12])
plt.legend(loc='middle left', fancybox=True, shadow=True, ncol=1)
plt.savefig('Calcium Number Densities.png')
plt.show()


'''Plotting: Plot 3 '''

plt.figure()
plt.title('Fractional $e^-$ Donation',size=14,fontweight='bold')
plt.xlabel('T (K)',size=12)
plt.ylabel('log(f($n_e$))',size=12)
plt.plot(list_of_temp,np.log10(f_ne_n11_list), label='$H^0$', color='orange')
plt.plot(list_of_temp,np.log10(f_ne_n21_list), label='$H^+$ ', color='r')
plt.plot(list_of_temp,np.log10(f_ne_n12_list), label='$He^0$ ', color='brown')
plt.plot(list_of_temp,np.log10(f_ne_n22_list), label='$He^+$ ', color='k')
plt.plot(list_of_temp,np.log10(f_ne_n32_list), label='$He^{++}$ ', color='k', linestyle=':')
plt.plot(list_of_temp,np.log10(f_ne_n13_list), label='$Ca^0$ ', color='b')
plt.plot(list_of_temp,np.log10(f_ne_n23_list), label='$Ca^+$ ', color='g')
plt.plot(list_of_temp,np.log10(f_ne_n33_list), label='$Ca^{++}$ ', color='g', linestyle=':')
plt.legend(loc='upper center', bbox_to_anchor=(.5, -.2),
          fancybox=True, shadow=True, ncol=3)
plt.ylim([-5,1])
plt.savefig('Fractional $e^-$ Donation.png')
plt.show()
