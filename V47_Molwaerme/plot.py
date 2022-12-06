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


def T(r):
    return  0.00134*r**2 + 2.296*r - 243.02

def alpha(T, a, b ):
    return a / T + b


def tab(name ,x, *args):
    i= 0
    s= 1
    for arg in args:
        s+=1
    d = " & " * (s - 1 )
    #name = input("Dateiname?")
    file= open(f'{name}.tex','w')
    while i <  len(x) :
        a = f"{x[i]}"
        a = a.replace(".", ",")
        a=a.replace("+/-", " \pm ")
        file.write(f" {a} " )
        if s==1:
         file.write(" \\\\ ")
        else:
            file.write(" & ")
        j=1
        for arg in args:
                    j+=1
                    b = f"{arg[i]}"
                    b = b.replace(".", ",")
                    b=b.replace("+/-", " \pm ")
                    file.write(f" {b} " )
                    if j != s:
                     file.write(" & ")
                    else:
                     file.write(" \\\\")
        file.write("\n")
        i+=1
    file.close()


############################
t, r_r, r_z, i , u = np.genfromtxt("data/data.txt", unpack= True) #in ohm
t =  t[1:] - t[:t.size -1] 
i= i /10**3 # von mA in A

T_r = T(r_r)
dT =  T_r[1:] - T_r[:T_r.size -1] 

########################
m = 342 # Masse in g
V_0 =  7.11 /10**6  # m**3 / mol
kappa= 140 # N /m**2
E = i[1:] * u[1:] *t 
M= 3.55 # in g / mol
c_p = E * M / m /  dT # /mol
##############
T_alpha, a = np.genfromtxt("data/alpha.txt", unpack = True)
p,cov = curve_fit(alpha, T_alpha ,a)
e = np.sqrt(np.diag(cov))

a_u= ufloat(p[0],e[0])
b_u= ufloat(p[1],e[1])

print(a_u)
print(b_u)


x = np.linspace(T_alpha[0], T_alpha[T_alpha.size  -1], 100)

plt.plot(T_alpha, a , "x", label = "Data")
plt.plot(x, alpha(x, p[0],p[1]) , label = "Fit")
plt.xlabel(r" $T$ in $K$")
plt.ylabel(r"$\alpha$ in $10^{-6} \si{\per\degreeCelsius}$")
plt.legend(loc="best")
plt.savefig("build/alpha.pdf")
plt.clf()
############################################
#alpha_T_u=alpha(T,a_u,b_u) * 10**(-6)
#c_v_u = c_p - 9 * alpha_T_u**2 * kappa * V_0 * T 

alpha_T=alpha(T_r[1:],p[0],p[1]) * 10**(-6)
c_v = c_p - 9 * alpha_T**2 * kappa * V_0 * T_r[1:]



plt.plot(T_r[1:], c_v , "x" )
plt.xlabel(r" $T$ in $\si{\degreeCelsius}$")
plt.ylabel(r"$c_v$ in $\si{\joule\per\mole\per\degreeCelsius}$")
plt.savefig("build/c_v.pdf")

#######################################
T = T_r[1:] + 273.15 
Cv= c_v[T<=170]
print(Cv)