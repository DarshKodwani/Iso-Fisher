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



def dnoise(l):
    return ((thetaarcmin*sigmaT)**2)*np.exp(l*(l+1)*thetarad**2/(8*np.log(2)))
    
def dnoiseP(l):
    return ((thetaarcmin*sigmaP)**2)*np.exp(l*(l+1)*thetarad**2/(8*np.log(2)))



#def dnoise(l):
#    return ( ((thetaarcmin*sigmaT)**(-2))*np.exp(-l*(l+1)*(thetarad**2)/(8*np.log(2)))
#    + ((thetaarcmin100*sigmaT100)**(-2))*np.exp(-l*(l+1)*(thetarad100**2)/(8*np.log(2)))
#    + ((thetaarcmin217*sigmaT217)**(-2))*np.exp(-l*(l+1)*(thetarad217**2)/(8*np.log(2)))
#    + ((thetaarcmin353*sigmaT353)**(-2))*np.exp(-l*(l+1)*(thetarad353**2)/(8*np.log(2))))**(-1)
#    
#def dnoiseP(l):
#    return (  ((thetaarcmin*sigmaP)**(-2))*np.exp(-l*(l+1)*(thetarad**2)/(8*np.log(2)))
#    + ((thetaarcmin100*sigmaP)**(-2))*np.exp(-l*(l+1)*(thetarad100**2)/(8*np.log(2)))
#    + ((thetaarcmin217*sigmaP217)**(-2))*np.exp(-l*(l+1)*(thetarad217**2)/(8*np.log(2)))
#    + ((thetaarcmin353*sigmaP353)**(-2))*np.exp(-l*(l+1)*(thetarad353**2)/(8*np.log(2))))**(-1)
#    
    
# 2) Setting up CAMB to obtain the C_ls

#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()

# Here we set the initial condition for the Fluctuations. 0,1,2,3,4,5 follow same notation as CAMB the different types of fluctuations. 
 
pars.scalar_initial_condition = 0

pars.InitialConditionVector = (0.,0.,0., 0., 1., 0., 0., 0., 0.)

#epsX here is the amplitude of the power we are adding to the power spectrum at a given k and kXpivot is the initial k 
pars.InitPower.set_params(ns=0.965, r=0, kXpivot= 0, epsX = 0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);

#calculate results for these parameters
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

ksmin = 0.0001
ksmax = 0.1
numks = 100
xs = np.linspace(ksmin,ksmax,numks)
ks = xs
Allcls = np.zeros((lmax+1, numks))
AllclsEE = np.zeros((lmax+1, numks))
AllclsBB = np.zeros((lmax+1,numks))
AllclsTE = np.zeros((lmax+1, numks))
epspower = 10**(-1)

#The adiabtic, unmodified, C_ls 

Unmod_ad_totCl = np.genfromtxt('/Users/darshkodwani/Documents/Darsh/Toronto/Research/CAMB-May2016/Unmod_ad_totCl.dat')

Norm = 1/(7.42835025e12)

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


# 3) Computing the fisher info kernel

FFReal = np.zeros((numks, numks)) 

countreala = 0

#This is the full Fisher info with the full covariance

for k in ks: 
    countrealb = 0
    for j in ks:
        xtemp1 = []
        for i in range(2, lmax):
            tempd1 = ((2*i+1)/2)*((Unmodcls[i,0] - Allcls[i, countreala])*(Unmodcls[i,0] - Allcls[i,countrealb])*(((Norm*Unmod_ad_totCl[i,2] + dnoiseP(i))/((Norm*Unmod_ad_totCl[i,4])**2 - (Norm*Unmod_ad_totCl[i,2] + dnoiseP(i))*(Norm*Unmod_ad_totCl[i,1]+ dnoise(i))))**2)
            + ((Unmodcls[i,0] - Allcls[i, countreala])*(Unmodcls[i,1] - AllclsEE[i, countrealb])*(((Norm*Unmod_ad_totCl[i,4])**2 - (Norm*Unmod_ad_totCl[i,2] + dnoiseP(i))*(Norm*Unmod_ad_totCl[i,1]+ dnoise(i)))**2))
            + ((Unmodcls[i,0] - Allcls[i, countrealb])*(Unmodcls[i,1] - AllclsEE[i, countreala])*(((Norm*Unmod_ad_totCl[i,4])**2 - (Norm*Unmod_ad_totCl[i,2] + dnoiseP(i))*(Norm*Unmod_ad_totCl[i,1]+ dnoise(i)))**2))
            - (Unmodcls[i,0] - Allcls[i, countreala])*(Unmodcls[i,3] - AllclsTE[i, countrealb])*((2*(Norm*Unmod_ad_totCl[i,2]+ dnoiseP(i))*(Norm*Unmod_ad_totCl[i,4]))/(((Norm*Unmod_ad_totCl[i,4])**2 - (Norm*Unmod_ad_totCl[i,2]+ dnoiseP(i))*(Norm*Unmod_ad_totCl[i,1]+ dnoise(i)))**2))
            - (Unmodcls[i,0] - Allcls[i, countrealb])*(Unmodcls[i,3] - AllclsTE[i, countreala])*((2*(Norm*Unmod_ad_totCl[i,2]+ dnoiseP(i))*(Norm*Unmod_ad_totCl[i,4]))/(((Norm*Unmod_ad_totCl[i,4])**2 - (Norm*Unmod_ad_totCl[i,2]+ dnoiseP(i))*(Norm*Unmod_ad_totCl[i,1]+ dnoise(i)))**2))
            + (Unmodcls[i,1] - AllclsEE[i, countreala])*(Unmodcls[i,1] - AllclsEE[i, countrealb])*(((Norm*Unmod_ad_totCl[i,1] + + dnoise(i))/((Norm*Unmod_ad_totCl[i,4])**2 - (Norm*Unmod_ad_totCl[i,2] + + dnoiseP(i))*(Norm*Unmod_ad_totCl[i,1] + + dnoiseP(i))))**2)
            - (Unmodcls[i,1] - AllclsEE[i, countreala])*((Unmodcls[i,3] - AllclsTE[i, countrealb]))*((2*(Norm*Unmod_ad_totCl[i,4])*(Norm*Unmod_ad_totCl[i,1]+ + dnoise(i)))/(((Norm*Unmod_ad_totCl[i,4])**2 - (Norm*Unmod_ad_totCl[i,2] + + dnoiseP(i))*(Norm*Unmod_ad_totCl[i,1] + + dnoise(i)))**2))
            - (Unmodcls[i,1] - AllclsEE[i, countrealb])*((Unmodcls[i,3] - AllclsTE[i, countreala]))*((2*(Norm*Unmod_ad_totCl[i,4])*(Norm*Unmod_ad_totCl[i,1]+ + dnoise(i)))/(((Norm*Unmod_ad_totCl[i,4])**2 - (Norm*Unmod_ad_totCl[i,2] + + dnoiseP(i))*(Norm*Unmod_ad_totCl[i,1] + + dnoise(i)))**2))
            + (Unmodcls[i,3] - AllclsTE[i, countreala])*(Unmodcls[i,3] - AllclsTE[i, countrealb])*((2*((Norm*Unmod_ad_totCl[i,4])**2 - (Norm*Unmod_ad_totCl[i,2]+ + dnoiseP(i))*(Norm*Unmod_ad_totCl[i,1] + + dnoise(i))))/((Norm*Unmod_ad_totCl[i,4]**2 - (Norm*Unmod_ad_totCl[i,2]+ dnoiseP(i))*(Norm*Unmod_ad_totCl[i,1]+ dnoise(i))))))
            xtemp1.append(tempd1)
        FFReal[countreala, countrealb] = sum(xtemp1)/(4*epspower)
        countrealb += 1
    countreala += 1

    
plt.figure()
CS = plt.contourf(ks, ks, FFReal,cmap = plt.cm.bone)
cbar = plt.colorbar(CS)
plt.show()