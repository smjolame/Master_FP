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
m = 342 / 10**3# Masse in kg
V_0 =  7.11 /10**6  # m**3 / mol
rho= 8.92 * 10**6 /10**3 # g in m**3
kappa= 140 # N /m**2
E = i[1:] * u[1:] *t 
M= 63.55 /10**3# in kg / mol
c_p = E * M / m /  dT # /mol
##############
T_alpha, a = np.genfromtxt("data/alpha.txt", unpack = True)
p,cov = curve_fit(alpha, T_alpha ,a)
e = np.sqrt(np.diag(cov))

a_u= ufloat(p[0],e[0])
b_u= ufloat(p[1],e[1])

#print(a_u)
#print(b_u)


x = np.linspace(T_alpha[0], T_alpha[T_alpha.size  -1], 100)

plt.plot(T_alpha, a , "x", label = "Data")
plt.plot(x, alpha(x, p[0],p[1]) , label = "Fit")
plt.xlabel(r" $T$ in $K$")
plt.ylabel(r"$\alpha$ in $10^{-6} \si{\per\kelvin}$")
plt.legend(loc="best")
plt.savefig("build/alpha.pdf")
plt.clf()
############################################

T = T_r[1:] + 273.15 #in K 
alpha_T_u=alpha(T,a_u,b_u) * 10**(-6)
c_v_u = c_p - 9 * alpha_T_u**2 * kappa * V_0 * T

alpha_T=alpha(T,p[0],p[1]) * 10**(-6)
c_v = c_p - 9 * alpha_T**2 * kappa * V_0 * T

#tab("build/cv",np.round(T,1),c_p,  c_v_u)

plt.plot(T, c_v , "x", label= r"$c_v$" )
plt.axhline(y = 3 * 8.3143, color = 'r', linestyle = '--', label= r"$3 R $")
plt.xlabel(r" $T$ in $\si{\kelvin}$")
plt.ylabel(r"$c_v$ in $\si{\joule\per\mole\per\kelvin}$")
plt.legend(loc="best")
plt.savefig("build/c_v.pdf")
"""
#######################################
Cv= c_v[T<=170]
tab("build/Deb", T[T<=170], Cv)
#############################################
Q_T= np.array(
[4.6,
 4.7,
 4.9,
 3.5,
 3.6,
 2.0,
 2.1])

Q = Q_T * T[T<=170]

#tab("build/Deb", np.round(T[T<=170],1), np.round(Cv,2), Q_T, np.round(Q,1))

Q_u= ufloat(np.mean(Q), np.std(Q))
print(Q_u)
####################################################################################
v_l = 4.7 * 10**3 #m/s
v_t = 2.26*10**3 #m/s
Qd= const.hbar / const.k *np.power(18 * np.pi**2 * rho * const.N_A / M  / (1/v_l**3 + 2/v_t**3), 1/3  )
print(Qd)
print(1-Q_u /345)
print(1-Qd /345)
"""