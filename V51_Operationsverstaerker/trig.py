
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import json
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import uarray
from uncertainties import unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,std_devs as stds)
from scipy.stats import sem

import matplotlib as mpl
mpl.use('pgf')
mpl.rcParams.update({
'font.family': 'serif',
'text.usetex': True,
'pgf.rcfonts': False,
'pgf.texsystem': 'lualatex',
'pgf.preamble': r'\usepackage{unicode-math}\usepackage{siunitx}',
})



data= np.genfromtxt('usb/Trigger.csv', delimiter=',', dtype=None)
f=data[:,0:1]
u1=data[:,1:2]
u2=data[:,2:3]

#print(f)

#plt.plot(f, u1 , label= "u1")
#plt.plot(f, u2 , label= "u2")
#plt.legend(loc="best")
#plt.savefig("build/trig.pdf")
#plt.clf()

print(np.amin(u2))
#########################
r1= 10e3
r2= 100e3
r3=1e3
c=1e-6

print(r2/4 / c / r1 /r3)

data= np.genfromtxt('usb/Signal.csv', delimiter=',', dtype=None)
f=data[:,0:1]
u1=data[:,1:2]
u2=data[:,2:3]

print(r1/r2)
print(np.amax(u1)/np.amax(u2))
print(np.amin(u1)/np.amin(u2))
