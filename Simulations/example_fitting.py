##############################################################
#          NFW and gNFW fit to halo density profile          #
##############################################################
#      we construct 50 equal-mass log radial bins between    #
#       3h < r < R_200c, h is the force softening scale      #
#       of the simulations                                   #
##############################################################

import numpy as np
import matplotlib.pyplot as plt
from Fitting_functions import best_fit_NFW, best_fit_gNFW, rho_function, errfunc_nfw, sigma_function

file = '.../Group_profile_99.dat'
output_file = '.../TNG100dark_99.dat'
data = np.loadtxt(file)
DM_ID2 = data[0:1]
Mvir2 = data[1:2].T
M500 = data[2:3].T
Rvir2 = data[3:4].T
stellar_r2 = data[4:54].T   #entire data point and remove 前尾部虚假数据
stellar_dis2 = data[54:104].T 
dm_r2 = data[104:154].T   #entire data point and remove 前尾部虚假数据
dm_dis2 = data[154:204].T  #entire data point and remove 前尾部虚假数据


Result_File = []
for i in range(10,11,1):
    rvir = Rvir2[i]
    M200 = Mvir2[i]
    mvir=(M200*10**10)
    v = dm_r2[i]                                                        
    w = dm_dis2[i]*10**10
    p1=best_fit_NFW(v,w,rvir,mvir) #fitting output parameters

###..........compute parameter error, goodness-of-fit and draw figure.............####
    dm_nfw=rho_function(v,p1,'nfw')

    # #.........compute residual distribution ..........###
    para_err=errfunc_nfw(p1,v,np.log10(w),rvir,mvir,w=1) #model-simuation values
    # The rms of the fitting residuals between the simulated log(y) and the model is used to evaluate the relative goodness of fit
    rms_NFW =sigma_function(np.log10(w),np.log10(dm_nfw))
    
    # # #.........figure..........###
    plt.figure(figsize=(8, 6))
    plt.loglog(v,w,'k.',markersize=8,label='DM profile') #draw DM profile
    plt.loglog(v,dm_nfw,'r-',lw=1.5,label='NFW') #fitting models profiles
    plt.xlabel(r'$r\ (kpc)$')
    plt.ylabel(r'DM mass density ($M_{\odot}kpc^{-3}$)')
    plt.legend(frameon=False,loc='upper right')
 
    # Save profiles
    plot_save_path = f'.../halo_profile_99_{i}.png'
    plt.savefig(plot_save_path)
    print(f"Plot {i} saved to:", plot_save_path)
    
    # Store results
    c_NFW = p1[0] #halo concentration
    r_200 = p1[1] #halo virial radius
    M_200 = p1[2] #halo virial mass
    Result_File.append([M_200, r_200, c_NFW])
# Convert to numpy array
data_array = np.array(Result_File)
np.savetxt(output_file, data_array, fmt='%.6e %.3f %.16f', comments='')
print("Data saved to:", output_file)

