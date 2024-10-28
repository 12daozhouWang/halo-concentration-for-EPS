# Halo Concentrations of Extended Press–Schechter

This Python script calculates the formation redshift distribution ($z_{f}$) of dark matter halos based on the conditional mass function of the excursion set theory, the mass variance ($\sigma_{M}$) of the linear density fluctuation field, and the halo concentration parameter by fitting the radial mass density distribution of TNG simulations using the NFW model. Additionally, it provides halo formation times based on the IllustrisTNG simulation merger tree. The concentration-mass relation of EPS theory is reproduced using these halo properties and their relationships. The distribution of the $z_{f}$ and $\sigma_{M}$ for the sample halos is calculated in detail within the EPS framework, following the physical processes outlined in Ling et al. (2003) and Lukic et al. (2007).

## Files Contained in the Package

- **EPS_theory/constants.py**: Mainly contains the cosmological parameters in Planck Collaboration XIII (2016).
- **EPS_theory/halo_mass_variance.py**: Calculates the variance of the linear density fluctuation field.
- **EPS_theory/halo_formation_distribution.py**: Calculates the formation redshift distribution ($z_{f}$) of dark matter halos based on the conditional mass function of EPS theory.
- **EPS_theory/growth_factor.py**: Calculates the growth factor for a given redshift.
- **EPS_theory/metroplis_hastings.py**: A Markov Chain Monte Carlo (MCMC) method to obtain a random sample sequence of halo formation redshifts.
- **EPS_theory/example_.py**: An example of calculating $z_{f}$ and $\sigma_{M}$ under the EPS theory.
- **Simulations/Fitting_functions.py**: The NFW profile fits the radial mass density profile of the halo.
- **Simulations/example_fitting.py**: An example of fitting a halo with the NFW profile.
- **Simulations/halo_formation_time.py**: Calculates halo formation time for IllustrisTNG simulations based on the merger tree.

## Reference

If you use this code in your work, we kindly request you to cite the following paper:
