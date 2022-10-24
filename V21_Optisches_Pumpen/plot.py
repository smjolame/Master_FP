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

from scipy.constants import mu_0 

def abw(exact,approx):
    return (exact-approx)*100/exact  #Abweichnung


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




#Daten imortieren:
f, B_s1, B_s2, B_h1, B_h2 = np.genfromtxt('data/B_felder.csv', delimiter=',', unpack=True)
U1, T1 = np.genfromtxt('data/Periodendauern1.csv', delimiter=',', unpack=True)
U2, T2 = np.genfromtxt('data/Periodendauern2.csv', delimiter=',', unpack=True)



bruchteile_period1 = np.array([1,3,4,5,6,9,10,10,11,10,13,15,16,18])
bruchteile_period2 = np.array([1,2,3,4,5,7,8,8,9,10,10,8,10,8])



T1 = T1/bruchteile_period1 
T2 = T2/bruchteile_period2 


#Test-Plots
plt.plot(f,B_s1*0.1+B_h1*0.3)
plt.plot(f,B_s2*0.1+B_h2*0.3)
plt.show()
plt.clf()
plt.plot(U1,T1)
plt.plot(U2,T2)
plt.show()



# Belder berechnen:

def B_helmholtz(N,R,I):
    return 8/(np.sqrt(125))*mu_0*N*I/R

