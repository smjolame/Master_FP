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

l, phi1, phi2 = np.genfromtxt("data/hochrein.txt", unpack= True)
l, phi1_136, phi2_136 = np.genfromtxt("data/136.txt", unpack= True)
l, phi1_1296, phi2_1296 = np.genfromtxt("data/hochrein.txt", unpack= True)

phi = np.abs(phi1 - phi2)/2
phi136 = np.abs(phi1_136 - phi2_136) / 2
phi1296 = np.abs(phi1_1296 - phi2_1296) / 2

plt.plot(l**2, phi)
plt.plot(l**2, phi136-phi)
plt.plot(l**2, phi1296-phi)
plt.savefig("build/plot.pdf")

##Curvefit
#def BeispielFunktion(x,a,b):
#    return a*x+b 
#params, cov = curve_fit(BeispielFunktion, x-Werte, y-Werte,sigma=fehler_der_y_werte,p0=[schätzwert_1,#schätzwert_2])
#a = ufloat(params[0],np.absolute(cov[0][0])**0.5)
#b = ufloat(params[1],np.absolute(cov[1][1])**0.5)
#
#
##Json
#Ergebnisse = json.load(open('data/Ergebnisse.json','r'))
#if not 'Name' in Ergebnisse:
#    Ergebnisse['Name'] = {}
#Ergebnisse['Name']['Name des Wertes']=Wert
#
#json.dump(Ergebnisse,open('data/Ergebnisse.json','w'),indent=4)
