######------ Halo formation redshift distribution------######
import numpy as np
import matplotlib.pyplot as plt
from constants import n, rho_0
from halo_mass_variance import sigma2, normalize_sigma8
from halo_formation_distribution import P
from metropolis_hastings import metropolis_hastings
from calculate_peak_and_fwhm import calculate_peak_and_fwhm
from growth_factor import linear_growth

# Load halo catalog data
data = np.loadtxt('.../TNG100dark_99.dat')
M0 = data[:, 0]
M_S = M0 / 2
z2 = 0.0
dz, ddot = linear_growth(z2)

# Numerical integration results in a redshift distribution
AA = normalize_sigma8()
M2, Mh = [], []
output_results = []
for i in range(0, 1):  # Assuming you want to loop over the first 5 halos
    R2 = (3.0 * M0[i] / (4.0 * np.pi * rho_0))**(1 / 3)
    Rh = (3.0 * M_S[i] / (4.0 * np.pi * rho_0))**(1 / 3)
    M2.append(sigma2(R2, AA))
    Mh.append(sigma2(Rh, AA))
    S2 = M2[-1]
    Sh = Mh[-1]

    sigmaM = np.sqrt(sigma2(R2, AA))
    sigmaM_z = sigmaM * dz  # Output results
    
    # Calculate P(z1) distribution
    z1_values = np.linspace(0.01, 6, 2000)
    P_values = P(S2, Sh, n, z2, z1_values)
    dP_dz1 = -np.gradient(P_values, z1_values)

    z_peak, left_z, right_z = calculate_peak_and_fwhm(z1_values, dP_dz1)

    # Generate random samples using Metropolis-Hastings algorithm
    def proposal_distribution(current_value):
        return np.random.normal(loc=current_value, scale=z_peak * 0.2)
    
    initial_value = z_peak
    iterations = 50000
    samples = metropolis_hastings(lambda z: np.interp(z, z1_values, dP_dz1), proposal_distribution, initial_value, iterations)

    # Sort and select the samples to be plotted, 
    # Here I select 50 data for each halo, but you can select more.
    sorted_samples = np.array(sorted(samples))
    value_interval = (sorted_samples[-1] - sorted_samples[0]) / 99
    selected_indices = [0, len(sorted_samples) - 1]
    for j in range(1, 99):
        value = sorted_samples[0] + j * value_interval
        index = np.abs(sorted_samples - value).argmin()
        selected_indices.append(index)
    selected_z1_values = [sorted_samples[i] for i in selected_indices]
    selected_z1_values.sort()

    # Plot histogram of Metropolis-Hastings samples
    plt.figure(figsize=(6, 6))
    plt.hist(samples, bins="auto", density=True, histtype='step', color='blue', edgecolor='black', label="Samples")
    plt.plot(z1_values, dP_dz1, 'r-', label="True Distribution")
    plt.xlabel('$z_{f}$')
    plt.ylabel('$dP / dz_{f}$')
    plt.axis([-0.2, 3, 0, 2.0])
    plt.legend()

    # Save the plot with a unique filename
    plot_save_path = f'.../Eps_theory/plot_TNG100dark_formation_99_halo_{i}.png'
    plt.savefig(plot_save_path)
    print(f"Plot {i} saved to:", plot_save_path)
    
# Store results in the output_results list
    output_results.append([i, M0[i], sigmaM_z, z_peak, left_z, right_z, *selected_z1_values])
data_array = np.array(output_results)
# Determine the number of columns
num_columns = data_array.shape[1]
# Create the format string based on the number of columns
fmt = '%.0f ' +'%.6e ' + ' '.join(['%.16f'] * (num_columns - 2))
# Save to file using scientific notation
output_file = '.../Eps_theory/TNG100dark_formation_99.dat'
np.savetxt(output_file, data_array, fmt=fmt, comments='')
print("Data saved to:", output_file)

