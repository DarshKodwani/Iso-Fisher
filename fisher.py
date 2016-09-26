import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from math import *
from pylab import *

ax = subplot(1,3,1)
bx = subplot(1,3,2)
cx = subplot(1,3,3)

full_planck = fits.open('/Users/darshkodwani/Documents/Darsh/Toronto/Research/Sound_modes/COM_PowerSpect_CMB_R202.fits')

#Low l TT spectrum

TTLOLUNB = full_planck[1].data

xTTLOLUNB = TTLOLUNB.field(0)
yTTLOLUNB = TTLOLUNB.field(1)
ypeTTLOLUNB = TTLOLUNB.field(2)
ymeTTLOLUNB = TTLOLUNB.field(3)
yTTLOLUNBerr =  [ymeTTLOLUNB, ypeTTLOLUNB]

ax.plot(xTTLOLUNB, yTTLOLUNB, label = 'Temperature power spectrum for low l')

ax.legend(loc = 'upper right', prop={'size':10}, borderaxespad=0.)
ax.errorbar(xTTLOLUNB, yTTLOLUNB, yTTLOLUNBerr )
ax.set_xlabel("$l$")
ax.set_ylabel("$l(l+1)C^{TT}_l/2 \pi \ \ [\mu K^2]$")

#High l TT spectrum 

TTHILUNB = full_planck[8].data

xTTHILUNB = TTHILUNB.field(0)
yTTHILUNB = TTHILUNB.field(1)
ypeTTHILUNB = TTHILUNB.field(2)
ymeTTHILUNB = -TTHILUNB.field(2)
yTTHILUNBerr =  [ymeTTHILUNB, ypeTTHILUNB]

bx.plot(xTTHILUNB, yTTHILUNB, label = 'Temperature power spectrum for high l')

bx.legend(loc = 'upper right', prop={'size':10}, borderaxespad=0.)
bx.errorbar(xTTHILUNB, yTTHILUNB, yTTHILUNBerr )
bx.set_xlabel("$l$")
bx.set_ylabel("$l(l+1)C^{TT}_l/2 \pi \ \ [\mu K^2]$")

#Combining the full TT spectrum

xTTFULLUNB = np.append(xTTLOLUNB, xTTHILUNB)
yTTFULLUNB = np.append(yTTLOLUNB, yTTHILUNB)
ypeTTFULLUNB = np.append(ypeTTLOLUNB, ypeTTHILUNB)
ymeTTFULLUNB = np.append(ymeTTLOLUNB, ymeTTHILUNB)
yTTFULLUNBerr =  [ymeTTFULLUNB, ypeTTFULLUNB]

cx.plot(xTTFULLUNB, yTTFULLUNB, label = 'Temperature power spectrum for all l')

cx.legend(loc = 'upper right', prop={'size':10}, borderaxespad=0.)
cx.errorbar(xTTFULLUNB, yTTFULLUNB, yTTFULLUNBerr )
cx.set_xlabel("$l$")
cx.set_ylabel("$l(l+1)C^{TT}_l/2 \pi \ \ [\mu K^2]$")

show()

#Importing the theoreical Cls from Class (change to whatever we need for initial conditions)
ClassCls = np.fromfile('/Users/darshkodwani/Documents/Darsh/Toronto/Computation/class_public-2.5.0/output/Bi_cl.dat',dtype=float, sep=' ')

temp = np.reshape(ClassCls, (2999,8))

CTTth = temp[:,1] #TT Cls from Class for Baryon isocurvature (no baryon to isocruavture set- so its set to default for this one)

#Computing the naive noise

Noise = np.sqrt(np.square(ypeTTFULLUNB) + np.square(ymeTTFULLUNB))

#Using the defintiion of detector noise used in "Forecasting isocurvature models with CMB lensing information" by Santos et al (2012). See Eq A9 in this paper.


def dnoise(l):
    return ((np.exp(-l*(l+1)*16.2754)/(441.868)) + (np.exp(-l*(l+1)*9.09078)/(255.736)) + (np.exp(-l*(l+1)*4.50842)/(857.317) + (np.exp(-l*(l+1)*4.50842)/(8040.69)) ))**(-1)


#noise = subplot(1,1,1)    
#noise.plot(xTTFULLUNB, Noise)
#noise.plot(xTTFULLUNB, dnoise(xTTFULLUNB))


#Computing the likelihood function


#Test arrays created to test the function

#test1= np.array([1,2,3,4,5,6])
#test2= np.array([10,20,30,40,50,60])
#test3 = np.array([0.1,2,4,54,5,8])
#test4 = np.array([2,8,43,42,2,2,4,56,6, 3, 5, 3, 5, 6, 2])

def lnL(fs):
    x = []
    for i in range(1, len(xTTFULLUNB)):
        t = (2*xTTFULLUNB[i] + 1)*((yTTFULLUNB[i]/(CTTth[i] + DetNoise(i))) + np.log((yTTFULLUNB[i]/(CTTth[i] + DetNoise(i)))) -3)
        x.append(t)
    likli = (fs/2)*sum(x)
    return likli

    

