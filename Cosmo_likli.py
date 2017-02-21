import matplotlib
import sys, platform, os
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import colors, ticker
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import camb
from camb import model, initialpower

# In this script we calculate the Likelihood for varying cosmology parameters w.r.t a fiducial cosmology given by Planck. 



# 1) Defining detector noise

#Using the defintiion of detector noise used in "Forecasting isocurvature models with CMB lensing information" by Santos et al (2012). See Eq A9 in this paper.
# Defining the noise paramaters - all quantities taken from the paper given above Table VII - we only take the 143 GHz channel

thetaarcmin = 7.1
thetaarcmin100 = 9.5
thetaarcmin217 = 5.0
thetaarcmin353 = 5.0
thetarad = thetaarcmin/3437.75
thetarad100 = thetaarcmin100/3437.75
thetarad217 = thetaarcmin217/3437.75
thetarad353 = thetaarcmin353/3437.75


sigmaT = 6.0016*(10**(-6))
sigmaT100 = 6.82*(10**(-6))
sigmaT217 = 13.0944*(10**(-6))
sigmaT353 = 40.1016*(10**(-6))

sigmaP = 11.4576*(10**(-6))
sigmaP100 = 10.9120*(10**(-6))
sigmaP217 = 26.7644*(10**(-6))
sigmaP353 = 81.2944*(10**(-6))

#Temperature noise
def dnoise(l):
    return ( ((thetaarcmin*sigmaT)**(-2))*np.exp(-l*(l+1)*(thetarad**2)/(8*np.log(2)))
    + ((thetaarcmin100*sigmaT100)**(-2))*np.exp(-l*(l+1)*(thetarad100**2)/(8*np.log(2)))
    + ((thetaarcmin217*sigmaT217)**(-2))*np.exp(-l*(l+1)*(thetarad217**2)/(8*np.log(2)))
    + ((thetaarcmin353*sigmaT353)**(-2))*np.exp(-l*(l+1)*(thetarad353**2)/(8*np.log(2))))**(-1)
    
#Polarization noise
def dnoiseP(l):
    return (  ((thetaarcmin*sigmaP)**(-2))*np.exp(-l*(l+1)*(thetarad**2)/(8*np.log(2)))
    + ((thetaarcmin100*sigmaP)**(-2))*np.exp(-l*(l+1)*(thetarad100**2)/(8*np.log(2)))
    + ((thetaarcmin217*sigmaP217)**(-2))*np.exp(-l*(l+1)*(thetarad217**2)/(8*np.log(2)))
    + ((thetaarcmin353*sigmaP353)**(-2))*np.exp(-l*(l+1)*(thetarad353**2)/(8*np.log(2))))**(-1)

    
# 2) Setting up CAMB to obtain the C_ls

#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()

# Here we set the initial condition for the Fluctuations.

# Setting this to 0 tells it to look at the vector in the following line.
pars.scalar_initial_condition = 1

#This vector follows the same notation as CAMB in terms of its components.
pars.InitialConditionVector = (0.,1.,1., 1., 1., 0., 0., 0., 0.)

#epsX here is the amplitude of the power we are adding to the power spectrum at a given k and kXpivot is the initial k 
pars.InitPower.set_params(ns=0.965, r=0, kXpivot= 0, epsX = 0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);

#calculate transfer functions for these parameters
results = camb.get_results(pars)

#get dictionary of CAMB power spectra
powers =results.get_cmb_power_spectra(pars)
for name in powers: print name
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print totCL.shape

ls = np.arange(totCL.shape[0])
lmax = 2000

#Setting up intial ks

ksmin = 10**(-2)
ksmax = 5*10**(-1)
numks = 10
xs = np.linspace(ksmin,ksmax,numks)
ks = xs
Allcls = np.zeros((lmax+1, numks))
AllclsEE = np.zeros((lmax+1, numks))
AllclsBB = np.zeros((lmax+1,numks))
AllclsTE = np.zeros((lmax+1, numks))
epspower = 10**(-1)

#The adiabtic, unmodified, C_ls 

Unmod_ad_totCl = np.genfromtxt('/Users/darshkodwani/Documents/Darsh/Toronto/Research/CAMB-May2016/Unmod_ad_totCl.dat')
Trans = np.genfromtxt('/Users/darshkodwani/Documents/Darsh/Toronto/Research/CAMB-May2016/test_transfer_out.dat')

Norm = 1/(7.42835025e12)

#Creating a set of Cls from normal/fiducial power spectrum 

add_initial_power1 = initialpower.InitialPowerParams()
add_initial_power1.set_params(kXpivot=0, epsX = 0)
Unmodcls = results.get_total_cls(lmax)
#Unmodcls = np.zeros((2001,4))

#Creating a set of Cls from modified power spectrum 

#Modifying H0

#H0likli = []
#H0s = np.linspace(60,70,numks)
#
#for i in H0s:
#    new_cosmo = pars.set_cosmology()
#    new_cosmo.set_cosmology(H0 = i)
#    new_results = camb.get_results(pars)
#    cl = new_results.get_total_cls(lmax)
#    tempH0 = []
#    for l in range(2, lmax):
#        ttempH0 = ((cl[l,0] - Norm*Unmod_ad_totCl[l,1])**(2))/((Norm*Unmod_ad_totCl[l,1] + dnoise(l))**(2))
#        tempH0.append(ttempH0)
#    H0likli.append(sum(tempH0))
#    #print pars
#    plt.loglog(np.arange(lmax+1),cl[:,0])
#    
#plt.xlim([2,lmax])
#plt.legend(H0s, loc='lower right')
#plt.show() 

#Modifying OmB

OmB = np.linspace(0.01,0.03,numks)
H0s = np.linspace(60,70,numks)
OmBlikli = []



for k in ks:
    add_initial_power = initialpower.InitialPowerParams()
    add_initial_power.set_params(kXpivot=k, epsX = epspower)
    results.power_spectra_from_transfer(add_initial_power)
    for j in H0s:    
        for i in OmB:
            new_cosmo = pars.set_cosmology()
            new_cosmo.set_cosmology(ombh2 = i, H0 = j)
            new_results = camb.get_results(pars)
            cl = new_results.get_total_cls(lmax)
            tempOmB = []
            for l in range(2, lmax):
                ttempOmB = ((cl[l,0] - Norm*Unmod_ad_totCl[l,1])**(2))/((Norm*Unmod_ad_totCl[l,1] + dnoise(l))**(2))
                tempOmB.append(ttempOmB)
            OmBlikli.append(sum(tempOmB))
    #plt.loglog(np.arange(lmax+1),cl[:,0])
    
plt.xlim([2,lmax])
plt.legend(OmB, loc='lower right')
plt.show() 

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    