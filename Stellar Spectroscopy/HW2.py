#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 14:43:32 2018

ASTR 545: Stellar Spectroscopy HW2


@author: oanavesa
"""
# import necessary imports
import matplotlib.pyplot as plt
import numpy as np
#%% Upload Txt Files

path = '/Users/oanavesa/Desktop/GradSchool/FirstYear/StellarSpectroscopy/HW/HW2/'
star3 = 'star3-spec.txt'

# wavelength (nm)
# flux density, f_lam, (erg*s^-1*cm^-2*Å^-1)
wavelength, flux = np.loadtxt(path+star3, unpack=True,usecols=[0,1])

#%% Problem 2 Part a

''' Problem 2; Part a: A plot of the six stellar spectra is also provided 
(“stars1-6.pdf”). On this plot 12 vertical dashed lines are drawn. 
Identify the twelve lines/features indicated by the vertical dashed lines. 
Plot “star3” and clearly label your line identifications on the plot.'''

# Creating the vertical lines for the features
plt.figure(figsize=(10,10))
plt.plot(wavelength,flux)
plt.axvline(x=430, color='red', linestyle='-.',label='G-Band') # G Band (CH)
plt.axvline(x=393, color='darkgreen', linestyle='-.',label='CaII K') # Ca II K
plt.axvline(x=396, color='black', linestyle='-.',label='CaII H') # Ca II H
plt.axvline(x=486, color='green', linestyle='-.', label='Hβ') # H beta
plt.axvline(x=656, color='darkred', linestyle='-.',label='Hα' ) # H alpha
plt.axvline(x=434, color='skyblue', linestyle='-.',label='Hγ') # H gamma
plt.axvline(x=410, color='brown', linestyle='-.',label='Hδ') # H delta
plt.axvline(x=388, color='purple', linestyle='-.',label='H8') # H 8
plt.axvline(x=383, color='orange', linestyle='-.',label='H9') # H 9
plt.axvline(x=364, color='gold', linestyle='-.',label='H∞') # H infinity
plt.axvline(x=589, color='pink', linestyle='-.',label='Na D') # Na D
plt.axvline(x=379, color='darkblue', linestyle='-.',label='H10') # H10

# plotting
plt.xlim((345,670))
plt.title('Spectral Classification of Star 3', fontweight='bold')
plt.ylabel('Flux Density (erg*s^-1*cm^-2*Å^-1)', fontweight='bold')
plt.xlabel('Wavelength (nm)', fontweight='bold')
plt.legend(title="Absorption Lines",fancybox = True,loc='upper right',
          ncol=2)
plt.savefig('Vesa_Star3_Lines.png')
plt.show()
#%% Problem 2; part c


''' Problem 2; Part c '''

''' Normalizing Star 3's Spectrum so that fν = fλ at λ = 550 nm  '''




# defining the constants and initial parameters
c = 3e8 # (m/s)
lam0 = 550 #(nm)
lam_m = lam0*1e-9 #(m)



# fν = fλ at λ = 550 nm 
# multiplying the flux values by c/5.5e-7 = constant
# fλ = fν *(c/λ^2)
fn = []
for i in range(len(wavelength)):
    fn2 = (flux[i]*((c)/(lam_m**2))) 
    fn.append(fn2)


# calculating a list of fν values index by index usng fn above and wavelength
# fν = fλ *(λ^2/c)
fnu = []
for j in range(len(wavelength)):
    fnucalc = fn[j] * (((wavelength[j]*(1e-9))**2)/c)
    fnu.append(fnucalc)



# plotting the original spectrum (fλ) and then the converted spectrum (fν) 
plt.figure()
plt.plot(wavelength,flux,label='fλ')
plt.plot(wavelength,fnu,label='fν')
plt.legend(labels=['fλ','fν'],fancybox = True,loc='upper right',
          ncol=2)
plt.xlim(350,670)
plt.title('Convert Spectrum Normalization',fontweight='bold')
plt.ylabel('Flux (erg s^−1 cm^−2  ̊A−1)',fontweight='bold')
plt.xlabel('Wavelength (nm)',fontweight='bold')
plt.savefig('Vesa_converted_spectrum.png')
plt.show()



#%% Upload files for Problem 3

# define path variables 
path = '/Users/oanavesa/Desktop/GradSchool/FirstYear/StellarSpectroscopy/HW/HW2/'
Vega = 'alpha_lyr-spec.txt'
cd = 'cd_34d241.txt'
Bband = 'filter-Bband.txt'
Vband = 'filter-Vband.txt'
hr15 = 'hr1544.txt'
hr34 = 'hr3454.txt'
Itt10 = 'ltt1020.txt'

# wavelength (nm)
# flux density, f_lam, (erg*s^-1*cm^-2*Å^-1)
# Johnson-Cousins B filter
waveA_Bband, R_lamB = np.loadtxt(path+Bband, unpack=True,usecols=[0,1])
# Johnson-Cousins V filter
waveA_Vband, R_lamV = np.loadtxt(path+Vband, unpack=True,usecols=[0,1])
# CD 34d241
wavecd, fluxcd = np.loadtxt(path+cd, unpack=True,usecols=[0,1])
# HR 1544
wavehr15, fluxhr15 = np.loadtxt(path+hr15, unpack=True,usecols=[0,1])
# HR 3454
wavehr34, fluxhr34 = np.loadtxt(path+hr34, unpack=True,usecols=[0,1])
# LTT 1020
waveItt10, fluxItt10 = np.loadtxt(path+Itt10, unpack=True,usecols=[0,1])
# Vega
wavevega, fluxvega = np.loadtxt(path+Vega, unpack=True,usecols=[0,1])

#%%
''' Problem 3: Compute the Johnson-Cousins B and V magnitudes on the Vega 
        system and then the B − V colors of of the fours star, CD 34d241, 
        HR 1544, HR 3454, and LTT 1020. '''



# creating the magnitude equation
def magnitudes(Fband,Fband_standard):
    '''Magnitude equation: inputs are the Fband calculated for each target
    star using numerical integration and Fband_standard calculated using numerical
    integration for Vega in the B and V band, separately '''
    M = -2.5*np.log10(Fband/Fband_standard)
    return M



# normalizing the Johnson-Cousins B and V filters
bband = np.trapz(R_lamB,waveA_Bband)
normb = R_lamB/bband
vband = np.trapz(R_lamV,waveA_Vband)
normv = R_lamV/vband



''' Computing flux for Vega = Fband_standard in the magnitudes equation '''



# Vega B filter
vegab = []
for k in range(len(wavevega)):
    #interpolating the values for Vega in B filter
    interpvega = np.interp(wavevega[k],waveA_Bband,normb)
    # R(λ)*f(λ)
    multt = interpvega*fluxvega[k]
    vegab.append(multt)
# Fband value Vega in B filer
vegaBFband = np.trapz(vegab,wavevega)



# Vega V filter
vegav = []
for k in range(len(wavevega)):
    #interpolating the values for Vega in V filter
    interpvega = np.interp(wavevega[k],waveA_Vband,normv)
# R(λ)*f(λ)
    multt = interpvega*fluxvega[k]
    vegav.append(multt)
# Fband value Vega in V filer
vegaVFband = np.trapz(vegav,wavevega)


''' Computing Fband for target stars '''



# LTT 1020 star
Itt10 = []
fluxItt10=fluxItt10/10**16
for i in range(len(waveItt10)):
    #interpolating the values for LTT 1020 in B filter
    interpltt = np.interp(waveItt10[i],waveA_Bband,normb)
    multt = interpltt*fluxItt10[i]
    Itt10.append(multt)
# Fband value LTT 1020 in B filer
Ltt1020BFband = np.trapz(Itt10,waveItt10)
  
Itt10v = []
for i in range(len(waveItt10)):
    #interpolating the values for LTT 1020 in V filter
    interpltt = np.interp(waveItt10[i],waveA_Vband,normv)
    multt = interpltt*fluxItt10[i]
    Itt10v.append(multt)
# Fband value LTT 1020 in V filer  
Ltt1020VFband = np.trapz(Itt10v,waveItt10)
  
    
    
# HR 1544
fluxhr15=fluxhr15/10**16
hr15 = []
for j in range(len(wavehr15)):
    #interpolating the values for HR1544 in B filter
    interphr15 = np.interp(wavehr15[j],waveA_Bband,normb)
    multt = interphr15*fluxhr15[j]
    hr15.append(multt)
# Fband value HR 1455 in B filer  
hr1544BFband = np.trapz(hr15,wavehr15)

hr15v = []
for j in range(len(wavehr15)):
    #interpolating the values for HR 1544 in V filter
    interphr15 = np.interp(wavehr15[j],waveA_Vband,normv)
    multt = interphr15*fluxhr15[j]
    hr15v.append(multt)
# Fband value HR 1455 in V filer  
hr1544VFband = np.trapz(hr15v,wavehr15)



# HR 3454
fluxhr34=fluxhr34/10**16
hr34 = []
for a in range(len(wavehr34)):
    #interpolating the values for HR 3454 in B filter
    interphr34 = np.interp(wavehr34[a],waveA_Bband,normb)
    multt = interphr34*fluxhr34[a]
    hr34.append(multt)
# Fband value HR 3454 in B filer  
hr3454BFband = np.trapz(hr34,wavehr34)

hr34v = []
for k in range(len(wavehr34)):
    #interpolating the values for HR 3454 in V filter
    interphr34 = np.interp(wavehr34[k],waveA_Vband,normv)
    multt = interphr34*fluxhr34[k]
    hr34v.append(multt)
# Fband value HR 3454 in V filer  
hr3454VFband = np.trapz(hr34v,wavehr34)



# CD 34d241
fluxcd=fluxcd/10**16
cd = []
for k in range(len(wavecd)):
    #interpolating the values for CD 34d241 in B filter
    interpcd = np.interp(wavecd[k],waveA_Bband,normb)
    multt = interpcd*fluxcd[k]
    cd.append(multt)
# Fband value CD 34d241 in B filer  
cd34d241BFband = np.trapz(cd,wavecd)

cdv = []
for k in range(len(wavecd)):
    #interpolating the values for CD 34d241 in V filter
    interpcd = np.interp(wavecd[k],waveA_Vband,normv)
    multt = interpcd*fluxcd[k]
    cdv.append(multt)
# Fband value CD 34d241 in V filer  
cd34d241VFband = np.trapz(cdv,wavecd)


''' Printing B and V magnitude values for each star '''


# Print magnitudes values for each star for B filter
print('\n')
print('B values are:')
print('LTT1020:', magnitudes(Ltt1020BFband,vegaBFband))
print('HR 1544:', magnitudes(hr1544BFband,vegaBFband))
print('HR 3454:',magnitudes(hr3454BFband,vegaBFband))
print('CD 34d21:',magnitudes(cd34d241BFband,vegaBFband))

# Print magnitudes values for each star for V filter
print('\n')
print('V values are:')
print('LTT1020:',magnitudes(Ltt1020VFband,vegaVFband))
print('HR 1544:',magnitudes(hr1544VFband,vegaVFband))
print('HR 3454:',magnitudes(hr3454VFband,vegaVFband))
print('CD 34d21:',magnitudes(cd34d241VFband,vegaVFband))


# Print B-V values for above
print('\n')
print('B-V values are:')
print('LTT1020:',magnitudes(Ltt1020BFband,vegaBFband) - magnitudes(Ltt1020VFband,vegaVFband))
print('HR 1544:',magnitudes(hr1544BFband,vegaBFband) - magnitudes(hr1544VFband,vegaVFband))
print('HR 3454:',magnitudes(hr3454BFband,vegaBFband) - magnitudes(hr3454VFband,vegaVFband))
print('CD 34d21:',magnitudes(cd34d241BFband,vegaBFband) - magnitudes(cd34d241VFband,vegaVFband))

