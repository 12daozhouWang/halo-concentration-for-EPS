######------ Calculation halo variance------######
import numpy as np
import sympy as sp
from scipy import integrate
from constants import h, omega_m, n, sigma_8

k, R = sp.symbols('k R')

# Define transfer function
def T(k):
    q = k / (omega_m * h)
    HS_value = sp.log(1.0 + 2.34*q) / (2.34*q) * (1.0 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**(-0.25)
    return HS_value

# Define the integral function
def s2_int(k, R, AA):
    pi = sp.pi
    inv2pi = 1.0 / (2.0 * pi * pi)
    trf = T(k)
    window_f = 3.0 * ((sp.sin(k * R) - k * R * sp.cos(k * R)) / ((k * R) ** 3))
    s2_int_value = inv2pi * k ** 2 * AA * k ** n * trf ** 2 * window_f ** 2
    return s2_int_value

# Optimize integration with scipy quad and sympy lambdify
s2_int_lambdified = sp.lambdify((k, R, 'AA'), s2_int(k, R, sp.symbols('AA')), 'numpy')

def sigma2(R, AA):
    result, _  = integrate.quad(s2_int_lambdified, 0, np.inf, args=(R, AA), epsabs=1.49e-8, epsrel=1.49e-8, limit=1000)
    return result

# Normalize P(k) to the given sigma_8:
def normalize_sigma8():
    AA = 1.0
    s2 = sigma2(8.0, AA)
    AA *= sigma_8 ** 2 / s2
    return AA
