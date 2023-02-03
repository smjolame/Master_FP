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

#f / Hz , U_A / V, -phi/Grad
f1, Ua1, phi1 = np.genfromtxt("data/M1.txt", unpack= True, delimiter= ",")
f2, Ua2, phi2 = np.genfromtxt("data/M2.txt", unpack= True, delimiter= ",")
f3, Ua3, phi3 = np.genfromtxt("data/M3.txt", unpack= True, delimiter= ",")

Ue= 0.139 #V
############
def exp(x, a, b):
    return a *x**b

###########
V1= Ua1/Ue
V2= Ua2/Ue
V3= Ua3/Ue
##########################################################################################
#R=100kohm

i1=6
V0_1,cov1 = np.polyfit(f1[:i1] , V1[:i1], 0, cov=True)
e1 = np.sqrt(np.diag(cov1))

s1, covv1= curve_fit(exp, f1[i1+1:], V1[i1+1:])
ee1 = np.sqrt(np.diag(covv1))

x1=np.linspace(f1[i1+1], f1[f1.size -1], 100)

a1= ufloat(s1[0], ee1[0])
b1= ufloat(s1[1], ee1[1])
V0_1u=ufloat(V0_1, e1[0])

fgrenz1_u= (V0_1u / np.sqrt(2) / a1)**(1/b1)
band1 = fgrenz1_u * V0_1u
##############

plt.plot(f1, V1, "x" , label=r"Messwerte")
plt.hlines(y=V0_1, xmin= f1[0], xmax= f1[i1],colors="orange", label = r"Ausgleichsgerade $V_0$")
plt.plot(x1, exp(x1, s1[0], s1[1]), label=r"Ausgleichskurve $V= a f^b$")
plt.plot((V0_1 / np.sqrt(2) / s1[0])**(1/s1[1]), V0_1 /np.sqrt(2), "x" , label= r"$f_\text{grenz}$")

plt.legend(loc="best")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$f$ in $\qty{}{\hertz}$")
plt.ylabel(r"$V$")
plt.title(r"$R_2=\qty{100}{\kilo\ohm}$")
plt.savefig("build/int100.pdf")
plt.clf()

##########################################################################################
#R=150kohm

i2=6
V0_2,cov2 = np.polyfit(f2[:i2] , V2[:i2], 0, cov=True)
e2 = np.sqrt(np.diag(cov2))

s2, covv2= curve_fit(exp, f2[i2+1:], V2[i2+1:])
ee2 = np.sqrt(np.diag(covv2))

x2=np.linspace(f2[i2+1], f2[f2.size -1], 100)

a2= ufloat(s2[0], ee2[0])
b2= ufloat(s2[1], ee2[1])
V0_2u=ufloat(V0_2, e2[0])

fgrenz2_u= (V0_2u / np.sqrt(2) / a2)**(1/b2)
band2 = fgrenz2_u * V0_2u
##############

plt.plot(f2, V2, "x" , label=r"Messwerte")
plt.hlines(y=V0_2, xmin= f2[0], xmax= f2[i2],colors="orange", label = r"Ausgleichsgerade $V_0$")
plt.plot(x2, exp(x2, s2[0], s2[1]), label=r"Ausgleichskurve $V= a f^b$")
plt.plot((V0_2 / np.sqrt(2) / s2[0])**(1/s2[1]), V0_2 /np.sqrt(2), "x" , label= r"$f_\text{grenz}$")

plt.legend(loc="best")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$f$ in $\qty{}{\hertz}$")
plt.ylabel(r"$V$")
plt.title(r"$R_2=\qty{150}{\kilo\ohm}$")
plt.savefig("build/int150.pdf")
plt.clf()
#################################################################################
#R=47kohm

i3=6
V0_3,cov3 = np.polyfit(f3[:i3] , V3[:i3], 0, cov=True)
e3 = np.sqrt(np.diag(cov3))

s3, covv3= curve_fit(exp, f3[i3+1:], V3[i3+1:])
ee3 = np.sqrt(np.diag(covv3))

x3=np.linspace(f3[i3+1], f3[f3.size -1], 100)

a3= ufloat(s3[0], ee3[0])
b3= ufloat(s3[1], ee3[1])
V0_3u=ufloat(V0_3, e3[0])

fgrenz3_u= (V0_3u / np.sqrt(2) / a3)**(1/b3)
band3 = fgrenz3_u * V0_3u
##############

plt.plot(f3, V3, "x" , label=r"Messwerte")
plt.hlines(y=V0_3, xmin= f3[0], xmax= f3[i3],colors="orange", label = r"Ausgleichsgerade $V_0$")
plt.plot(x3, exp(x3, s3[0], s3[1]), label=r"Ausgleichskurve $V= a f^b$")
plt.plot((V0_3 / np.sqrt(2) / s3[0])**(1/s3[1]), V0_3 /np.sqrt(2), "x" , label= r"$f_\text{grenz}$")

plt.legend(loc="best")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$f$ in $\qty{}{\hertz}$")
plt.ylabel(r"$V$")
plt.title(r"$R_2=\qty{47}{\kilo\ohm}$")
plt.savefig("build/int47.pdf")
plt.clf()

#print("Leerlaufverstärkung")
print(1-V0_1u/100)
print(1-V0_2u/150)
print(1-V0_3u/47)

#print("fitpararmeter")
#print("a1: ", a1, "b1: ", b1)
#print("a2: ", a2, "b2: ", b2)
#print("a3: ", a3, "b3: ", b3)
#
#print("grenfrequenz")
#print(fgrenz1_u)
#print(fgrenz2_u)
#print(fgrenz3_u)
#
#print("bandprodukt")
#print(band1)
#print(band2)
#print(band3)
################################################
#phasenabhängigkeit

plt.plot(f1, phi1 , "x")
plt.xscale("log")
plt.xlabel(r"$f$ in $\qty{}{\hertz}$")
plt.ylabel(r"$\phi$ in \qty{}{\degree}")
plt.title(r"$R_2=\qty{100}{\kilo\ohm}$")
plt.savefig("build/int_f100.pdf")
plt.clf()

plt.plot(f2, phi2 , "x")
plt.xscale("log")
plt.xlabel(r"$f$ in $\qty{}{\hertz}$")
plt.ylabel(r"$\phi$ in \qty{}{\degree}")
plt.title(r"$R_2=\qty{150}{\kilo\ohm}$")
plt.savefig("build/int_f150.pdf")
plt.clf()

plt.plot(f3, phi3 , "x")
plt.xscale("log")
plt.xlabel(r"$f$ in $\qty{}{\hertz}$")
plt.ylabel(r"$\phi$ in \qty{}{\degree}")
plt.title(r"$R_2=\qty{47}{\kilo\ohm}$")
plt.savefig("build/int_f47.pdf")
plt.clf()