import matplotlib
import sys, platform, os
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.cm as cm
#uncomment this if you are running remotely and want to keep in synch with repo changes
#if platform.system()!='Windows':
#    !cd $HOME/git/camb; git pull github master; git log -1
import camb
from camb import model, initialpower

print('Using CAMB installed at '+ os.path.realpath(os.path.join(os.getcwd(),'..')))
sys.path.insert(0,os.path.realpath(os.path.join(os.getcwd(),'..')))


# This script is an attempt to find the Fisher info kernel for adiabatic/non-adiabatic modes 

# 1) Importing Planck data

full_planck = fits.open('/Users/darshkodwani/Documents/Darsh/Toronto/Research/Sound_modes/COM_PowerSpect_CMB_R202.fits')

#Low l TT spectrum

TTLOLUNB = full_planck[1].data

xTTLOLUNB = TTLOLUNB.field(0)
yTTLOLUNB = TTLOLUNB.field(1)
ypeTTLOLUNB = TTLOLUNB.field(2)
ymeTTLOLUNB = TTLOLUNB.field(3)
yTTLOLUNBerr =  [ymeTTLOLUNB, ypeTTLOLUNB]
####################

#High l TT spectrum 

TTHILUNB = full_planck[8].data

xTTHILUNB = TTHILUNB.field(0)
yTTHILUNB = TTHILUNB.field(1)
ypeTTHILUNB = TTHILUNB.field(2)
ymeTTHILUNB = -TTHILUNB.field(2)
yTTHILUNBerr =  [ymeTTHILUNB, ypeTTHILUNB]

#Combining the full TT spectrum

xTTFULLUNB = np.append(xTTLOLUNB, xTTHILUNB)
yTTFULLUNB = (np.append(yTTLOLUNB, yTTHILUNB))*(10**(-12))
ypeTTFULLUNB = np.append(ypeTTLOLUNB, ypeTTHILUNB)
ymeTTFULLUNB = np.append(ymeTTLOLUNB, ymeTTHILUNB)
yTTFULLUNBerr =  [ymeTTFULLUNB, ypeTTFULLUNB]


# 2) Defining detector noise

#Using the defintiion of detector noise used in "Forecasting isocurvature models with CMB lensing information" by Santos et al (2012). See Eq A9 in this paper.
# Defining the noise paramaters - all quantities taken from the paper given above Table VII - we only take the 143 GHz channel

thetaarcmin = 7.1
thetarad = thetaarcmin/3437.75
sigmaT = 6.0016*(10**(-6))


def dnoise(l):
    return ((thetaarcmin*sigmaT)**2)*np.exp(l*(l+1)*thetarad**2/(8*np.log(2)))
    
# 3) Setting up CAMB to obtain the C_ls

#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
lSampleBoost = 1
#epsX here is the amplitude of the power we are adding to the power spectrum at a given k and kXpivot is the initial k 
pars.InitPower.set_params(ns=0.965, r=0, kXpivot= 0, epsX = 0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);

#calculate results for these parameters
results = camb.get_results(pars)

#get dictionary of CAMB power spectra
powers =results.get_cmb_power_spectra(pars)
for name in powers: print name

#plot the total lensed CMB power spectra versus unlensed, and fractional difference
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print totCL.shape

ls = np.arange(totCL.shape[0])
lmax = 2000

#Setting up intial ks

ksmin = 4
ksmax = 1
numks = 10
xs = np.linspace(ksmin,ksmax,numks)
ks = 10**(-xs)
Allcls = np.zeros((lmax+1, numks))
epspower = 1

#Creating a set of Cls from normal/fiducial power spectrum 

add_initial_power1 = initialpower.InitialPowerParams()
add_initial_power1.set_params(kXpivot=0, epsX = 0)
Unmodcls = results.get_total_cls(lmax)

#Creating a set of Cls from modified power spectrum 

count = 0 

for k in ks:
    add_initial_power = initialpower.InitialPowerParams()
    add_initial_power.set_params(kXpivot=k, epsX = epspower)
    results.power_spectra_from_transfer(add_initial_power)
    cl = results.get_total_cls(lmax)
    Allcls[:,count] = cl[:,0]
    count += 1
    plt.loglog(np.arange(lmax+1),cl[:,0])

plt.xlim([2,lmax])
plt.legend(ks, loc='lower right');'''''' 


# 4) Computing the fisher info kernel

Fishinfo = np.zeros((numks,numks))
Fishyinfo = []

count3 = 0 
    
for j in ks:
    countt = 0
    for k in ks:
        x = []
        for i in range(1, lmax):
            onedelta = ((2*i+1)/2)*((2*(Unmodcls[i,0])**2 - (Allcls[i,countt] - Allcls[i,count3])**2)*( -yTTFULLUNB[i]/((Unmodcls[i,0] + dnoise(i))**2) - 1/(Unmodcls[i,0] + dnoise(i))) 
            + (Unmodcls[i,0] - Allcls[i,countt])*(Unmodcls[i,0] - Allcls[i,count3])*(2*yTTFULLUNB[i]/((Unmodcls[i,0] + dnoise(i))**3) + 1/((Unmodcls[i,0] + dnoise(i))**2)))
            x.append(onedelta)
        Fishinfo[count3,countt] = sum(x) 
        countt += 1
    count3 += 1
    
plt.figure()
CS = plt.contourf(ks, ks, Fishinfo,cmap = plt.cm.bone)
cbar = plt.colorbar(CS)
plt.show()

