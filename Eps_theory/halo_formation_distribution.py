######------ Calculation halo formation redshift distribution------######
import numpy as np
from scipy import integrate
from constants import omega_m, delta_c

# Calculate deltac_t(z) as per Carroll, Press & Turner (1992) and Lin et al. (2003)
def deltac_t(z):
    omega_z = omega_m * ((1.0 + z)**3 / (1 - omega_m + (1.0 + z)**3 * omega_m))
    g_omega0 = 2.5 * omega_m * (1 / 70 + (209 / 140) * omega_m - omega_m**2 / 140 + omega_m**(4 / 7))**(-1)
    g_omegaz = 2.5 * omega_z * (1 / 70 + (209 / 140) * omega_z - omega_z**2 / 140 + omega_z**(4 / 7))**(-1)
    deltac_t_value = (1.0 + z) * (g_omega0 / g_omegaz) * delta_c
    return deltac_t_value

# Define P(z1) and related functions
def M_S1(S1, S2, Sh, n):
    S_bar = (S1 - S2) / (Sh - S2)
    alpha = (n + 3) / 3
    return (1 + (2**alpha - 1) * S_bar)**(1 / alpha)

def f_S1(S1, z1, S2, Sh, n, z2):
    delta1 = deltac_t(z1)
    delta2 = deltac_t(z2)
    nu = (delta1 - delta2) / np.sqrt(S1 - S2)
    result = M_S1(S1, S2, Sh, n) * (1.0 / np.sqrt(2.0 * np.pi)) * nu * np.exp(-0.5 * nu**2) / (S1 - S2)
    return result

def P(S2, Sh, n, z2, z1_values):
    def P_z1(z1):
        f_z1, _ = integrate.quad(lambda S1: f_S1(S1, z1, S2, Sh, n, z2), S2, Sh)
        return f_z1

    return [P_z1(z1) for z1 in z1_values]
