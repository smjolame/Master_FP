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


def lin(x,a,):
    return a*x

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



l, phi1, phi2 = np.genfromtxt("data/hochrein.txt", unpack= True)
l, phi1_136, phi2_136 = np.genfromtxt("data/136.txt", unpack= True)
l, phi1_1296, phi2_1296 = np.genfromtxt("data/hochrein.txt", unpack= True)

dhoch= 5.11 /10**3
d136= 1.36 / 10 **3
d1296=1.296 /10**3


phi = np.abs(phi1 - phi2)/2 / dhoch
phi136 = np.abs(phi1_136 - phi2_136) / 2 /d136
phi1296 = np.abs(phi1_1296 - phi2_1296) / 2 /d1296

#tab("build/tab1", l, phi1, phi2 , np.round(phi))
#tab("build/tab2", l,  phi1_136, phi2_136, np.round(phi136 ))
#tab("build/tab3", l,  phi1_1296, phi2_1296, np.round(phi1296 ))


plt.plot(l**2, phi, "x", label="hochrein")
plt.plot(l**2, phi136,"x",  label= r"$N = 1.2 \cdot 10^{18} \si{\per\cubic\cm}$")
plt.plot(l**2, phi1296,"x", label= r"$N= 2.8 \cdot 10^{18} \si{\per\cubic\cm}$")
plt.xlabel(r"$\lambda^2$ in $ \si{\um \squared}")
plt.ylabel(r"$\frac{\Theta}{d}$")
plt.legend( loc="best")
plt.savefig("build/plot1.pdf")
plt.clf()

a136, cov136 = curve_fit(lin, (l)**2 ,phi136-phi)
a_e136 = np.sqrt(cov136[0][0])
a_u136= ufloat(a136, a_e136)
#print(a_u136)

a1296, cov1296 = curve_fit(lin, l**2 ,phi1296-phi)
a_e1296 = np.sqrt(cov1296[0][0])
a_u1296= ufloat(a1296, a_e1296)
#print(a_u1296)

x=np.linspace(0, l[l.size -1], 100)

plt.plot(l**2, phi136-phi, "rx", label= r"$N = 1.2 \cdot 10^{18} \si{\per\cubic\cm}$")
plt.plot(x**2, a136 * x**2, label=r"Ausgleichsgerade, $N = 1.2 \cdot 10^{18} \si{\per\cubic\cm}$" )
plt.plot(l**2, phi1296-phi, "gx", label= r"$N = 2.8 \cdot 10^{18} \si{\per\cubic\cm}$")
plt.plot(x**2, a1296 * x**2, label=r"Ausgleichsgerade,$N = 2.8 \cdot 10^{18} \si{\per\cubic\cm}$")
plt.legend(loc="best")
plt.xlabel(r"$\lambda^2$ in $ \si{\um \squared}")
plt.ylabel(r"$\frac{\Theta}{d} - \frac{\Theta_\text{hochrein}}{d} $")
plt.savefig("build/plot2.pdf")
plt.clf()


#################################################
a_u136=a_u136 * 10**12
a_u1296=a_u1296 * 10**12
N136 = 1.2 *10**(18) * 10**6
N1296 = 2.8 *10**(18) * 10**6
B = 406 /10*3

n = 3.57
m136= unp.sqrt(const.e**3 * N136 *B / 8 / np.pi**2 / const.epsilon_0 / const.c**3 /n/ a_u136)
#print(m136)
#print(m136/const.m_e)
m1296= unp.sqrt(const.e**3 * N1296 *B / 8 / np.pi**2 / const.epsilon_0 / const.c**3 /n/ a_u136)
#print(m1296)
#print(m1296/const.m_e)
print(1 - m136/const.m_e / 0.067)
print(1 - m1296/const.m_e / 0.067)
