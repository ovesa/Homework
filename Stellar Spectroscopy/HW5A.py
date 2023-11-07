#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 30 20:26:25 2018

ASTR 545: Homework 5a. The calculations for the problems.

@author: oanavesa
"""
# imports
import numpy as np
from scipy import optimize

# constants
T = 10000 #K
P = 100 #dynes cm^-3
k = 1.380658e-16 # cm^2 g s^-2 K^-1
A_H = 1.00794
A_He = 4.002602
x_H = 0.7
x_He = 0.3
f_11 = 0.15
f_21 = 0.85
f_12 = 0.01
f_22 = 0.89
f_32 = 0.10
a = 7.56*10**(-15) # erg cm^-3 K^-4
m_h =1.6737236*10**(-24)
m_he = 6.6464764*10**(-24)
me = 9.10938*10**(-28)

# equations


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

def solve_for_number_densities_of_ionization_stages(f_jk,n_k):
    n_jk = f_jk*n_k
    return n_jk

def electron_number_density(n_1,n_2,f11,f21,f12,f22,f32):
    electron_num_dens_H = n_1*(0)*f11 + n_1*1*f21
    electron_num_dens_He = n_2*(0)*f12 + n_2*1*f22+n_2*2*f32
    total = electron_num_dens_H+electron_num_dens_He
    return total

def electron_pressure(k,T,ne):
    e_pressure = ne*k*T
    return e_pressure

def gas_pressure(k,T,n):
    g_pressure = k*T*n
    return g_pressure

def radiation_pressure(a,T):
    rad_pressure = (a/3)*T**4
    return rad_pressure

def mu_n(x_1,A_1,x_2,A_2):
    alpha_1 = x_1/A_1
    alpha_2 = x_2/A_2
    mu = alpha_1+alpha_2
    mu_n = (mu)**(-1)
    return mu_n

def mu_e(x_1,A_1,x_2,A_2,f11,f21,f12,f22,f32):
    alpha_1 = x_1/A_1
    alpha_2 = x_2/A_2
    mu_H = alpha_1*(0)*f11 + alpha_1*1*f21
    mu_He = alpha_2*(0)*f12 + alpha_2*1*f22+alpha_2*2*f32
    total = mu_H +mu_He
    total2 = total**(-1)
    return total2

def total_mu(mu_n_val,mu_e_val):
    inverse_mun = 1/mu_n_val
    inverse_mue = 1/mu_e_val
    inverse_mu_tot = inverse_mun + inverse_mue
    mu_tot = 1/inverse_mu_tot
    return mu_tot
    

def nuclear_mass_density(m_H,m_He,n_H,n_He):
    p_n = m_H*n_H + m_He*n_He
    return p_n

def   _density(me,n_h,n_he):
    e_dens = me*(1*n_h  + 2*n_he)
    return e_dens
    

#Problem 2
#part a) nuclear particle number density
nuclear_particle_density = solve_for_nuclear_particle_density(P,k,T)
# part b) abundance fractions
abundance_h = abundance_fraction(x_H,A_H,x_He,A_He)
abundance_he = abundance_fraction(x_He,A_He,x_H,A_H)
# part c) number densities
n_1 = solve_for_n_k(abundance_h,nuclear_particle_density)
n_2 = solve_for_n_k(abundance_he,nuclear_particle_density)
# part d) number densities of all ionization stages
n_11 = solve_for_number_densities_of_ionization_stages(f_11,n_1)
n_12 = solve_for_number_densities_of_ionization_stages(f_12,n_2)
n_22 = solve_for_number_densities_of_ionization_stages(f_22,n_2)
n_21 = solve_for_number_densities_of_ionization_stages(f_21,n_1)
n_32 = solve_for_number_densities_of_ionization_stages(f_32,n_2)
# part e) electron number density
total_electron_density = electron_number_density(n_1,n_2,f_11,f_21,f_12,f_22,f_32)
# part f) total particle number density
total_particle_density = nuclear_particle_density + total_electron_density
# part g) electron pressure
electron_pressure_value = electron_pressure(k,T,total_electron_density)
# part g) gas pressure
gas_pressure_value = gas_pressure(k,T,total_particle_density)
# part h) radiation pressure
radiation_pressure_value = radiation_pressure(a,T)
# part h) total pressure
total_pressure = radiation_pressure_value + gas_pressure_value
# part h) ratio of radiation and total pressure
Beta = gas_pressure_value/total_pressure
# part i) molecular weights of nuclear particles
mu_n_value =  mu_n(x_H,A_H,x_He,A_He)
# part i) molecular weights of electron
mu_e_value = mu_e(x_H,A_H,x_He,A_He,f_11,f_21,f_12,f_22,f_32)
# part i) total molecular weight
total_mu_value = total_mu(mu_n_value,mu_e_value)
# part j) nuclear particle mass density
nuclear_particle_mass_density = nuclear_mass_density(m_h,m_he,n_1,n_2)
# part j) electron mass density
electron_density_value= electron_density(me,n_1,n_2)
# part j) total mass density
total_density = (1 + (electron_density_value/nuclear_particle_mass_density))*nuclear_particle_mass_density
# part k)
percentage = (n_11/total_electron_density)*100

