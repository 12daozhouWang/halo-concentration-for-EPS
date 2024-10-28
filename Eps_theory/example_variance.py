#####----- Halo mass variance------#####
import numpy as np
from constants import rho_0
from growth_factor import linear_growth
from halo_mass_variance import sigma2, normalize_sigma8

# Redshift and growth factor
z = 0.0
dz, ddot = linear_growth(z)
# Load halo catalog
data = np.loadtxt('.../TNG100dark_99.dat')
M0 = data[:, 0]  # M_sun/h
# Normalize sigma8
AA = normalize_sigma8()

# Numerical integration results in a redshift distribution
Result_File = []
for i in range(0,10,1):
    # Calculate sigmaM for every M_200c
    R2 = (3.0 * M0[i] / (4.0 * np.pi * rho_0)) ** (1/3)
    sigmaM = np.sqrt(sigma2(R2, AA))
    sigmaM_z = sigmaM * dz  # Output results
    # Store results
    Result_File.append([M0[i], sigmaM_z])

# Convert to numpy array
data_array = np.array(Result_File)
# Save to file using scientific notation
output_file = '.../Eps_theory/TNG100dark_99.dat'
np.savetxt(output_file, data_array, fmt='%.6e %.16f', comments='')
print("Data saved to:", output_file)
