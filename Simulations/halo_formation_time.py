##############################################################
#     Halo formation time for IllustrisTNG simulations       #
##############################################################
#      The formation time of the halo, defined by the        #
#      relation M200c(a=a_form)=M200c(a=1)/2. We also        #
#      use a similar method to calculate the formation       #
#      time of halo under <a.                                #

#      Note: 1. To calculate the halo formation time, you    #
#      need to obtain the TNG dark halo's merger trees,      #
#      Groupcat and Offsets.                                 #
#            2. For the core code, see                       #
#           https://www.tng-project.org/data/docs/scripts/   #
##############################################################

import illustris_python as ill
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

basePath = '.../TNG-100-dark/Groupcat'  # Data path
snapshot_redshift_file = '.../Shapshot_redshift.tab'# Snapshot to redshift mapping file
output_file = '...TNG100dark_formation_99.dat'  # Output file for results
snapshot = 99  # Snapshot number

GroupFirstSub = ill.groupcat.loadHalos(basePath, snapshot, fields=['GroupFirstSub'])
fields = ['SubhaloMass', 'SubfindID', 'SnapNum', 'Group_M_Crit200','GroupMass']
snap_to_redshift = np.loadtxt(snapshot_redshift_file)     # Load the snapshot to redshift mapping
snap_to_redshift_dict = {int(row[0]): row[1] for row in snap_to_redshift} # Convert to a dictionary for easy lookup

output_results = []  # Initialize a list to store results
for i in range(0, 1):
    tree = ill.sublink.loadTree(basePath, snapshot, GroupFirstSub[i], fields=fields, onlyMPB=True)
    # Convert SnapNum to redshifts
    redshifts = np.array([snap_to_redshift_dict.get(snap, np.nan) for snap in tree['SnapNum']])
    # Interpolation for Group_M_Crit200
    interp_M200c = interp1d(np.log10(1 / (1 + redshifts)), np.log10(tree['Group_M_Crit200']), kind='linear')
    # Generate a continuous range of scale factors for interpolation
    x_continuous = np.linspace(np.log10(1 / (1 + redshifts.min())), np.log10(1 / (1 + redshifts.max())), 5000)
    y_M200c = interp_M200c(x_continuous)

    initial_mass = 10**y_M200c[0]
    M_200c_half = np.log10(initial_mass / 2)     # Compute half of the current Group_M_Crit200
    condition = y_M200c >= M_200c_half           # Condition where the mass is greater than or equal to half
    if np.any(condition):
        first_index = np.where(condition)[0][-1] # Get the last occurrence that meets the condition
        a_form = x_continuous[first_index]       # formation time a_form
        mass_value = tree['Group_M_Crit200'][0]  # Current mass
    redshift =  1/10**a_form-1                   # Calculate the corresponding redshift
    
    ####--------- Plotting the results--------####
    plt.figure(figsize=(8, 6))
    plt.plot(np.log10(1/(1+redshifts)), np.log10(tree['Group_M_Crit200']), 'k.', label="Group_M_Crit200")
    plt.plot(np.log10(1/(1+redshifts)), np.log10(tree['SubhaloMass']), 'b-', label="SubhaloMass")
    plt.plot(x_continuous, y_M200c, 'r-', label='Interpolated Group_M_Crit200')
    # Indicate the formation time on the plot
    plt.axhline(y=M_200c_half, color='r', linestyle='--')
    plt.axvline(x=a_form, color='r', linestyle='--', label='Formation Redshift')
    plt.xlabel(r'$\log_{10} a_{form}$')
    plt.ylabel(r'$\log_{10}$ Group_M_Crit200 [code units]')
    plt.legend(loc='best')
    # Save the plot with a unique filename
    plot_save_path = f'.../TNG100dark_formation_99_halo_{i}.png'
    plt.savefig(plot_save_path)
    plt.close()  # Close the plot to free up memory
    print(f"Plot {i} saved to:", plot_save_path)
    
    ####--------- Store results in the output_results list--------####
    M = tree['Group_M_Crit200'][0] #M_200/10^10 Msun/h
    a = 10**a_form
    z_f = redshift
    output_results.append([M, a, z_f])
# Save the output results to a file
output_array = np.array(output_results)
np.savetxt(output_file, output_array, fmt='%.6e %.16f %.16f', comments='')
print("Data saved to:", output_file)
