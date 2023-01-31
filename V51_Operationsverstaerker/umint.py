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


def exp(x, a, b):
    return a *x**b


fi, Uei, Uai = np.genfromtxt("data/int_sin.txt", unpack= True, delimiter= ",")
Vi= Uai/ Uei

pi, covi= curve_fit(exp, fi, Vi )
ei = np.sqrt(np.diag(covi))

xi=np.linspace(fi[0], fi[fi.size -1], 100)

a = ufloat(pi[0], ei[0])
b = ufloat(pi[1], ei[1])
print("a:", a)
print("b:", b)

plt.plot(fi, Vi, "x" , label=r"Messwerte")
plt.plot(xi, exp(xi, pi[0], pi[1]), label=r"Ausgleichskurve $V= a f^b$")

plt.legend(loc="best")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$f$ in $\qty{}{\hertz}$")
plt.ylabel(r"$V$")
plt.savefig("build/umint.pdf")
plt.clf()