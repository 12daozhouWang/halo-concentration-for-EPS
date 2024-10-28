######------ Constant module------######

import numpy as np
# Cosmology parameters from Planck Collaboration XIII 2016 (TNG simulations)
h = 0.6774       # Hubble constant in 100 km/s/Mpc
omega_m = 0.3089  # Total matter density
n = 0.9667        # Spectral index
sigma_8 = 0.816  # Variance of the linear density field
delta_c = 1.686   # Linear overdensity at virialization
# Calculate critical density and current density
rho_c = 3.0e4 / (8.0 * np.pi * 6.67430e-11) * 1.5513826e-2  # rho_c = 3*H^2/(8pi*G), in h^2 M_sun Mpc^(-3)
rho_0 = rho_c * omega_m # Matter density