import numpy as np
from scipy.integrate import quad
from scipy import optimize

def rho_function(r, params, model='nfw'):
    # NFW PROFILE
    if model == 'nfw':
        c, Rvir, Mvir = params
        rs = Rvir / c
        x = r / rs
        qb = Mvir * c**3 / (4. * np.pi * Rvir**3 * (np.log(1. + c) - c / (1. + c)))
        return qb / (x * (1. + x)**2)
    
    # GENERALIZED NFW PROFILE (gNFW)
    elif model == 'gnfw':
        c, a, Rvir, Mvir = params
        rs = Rvir / c
        x = r / rs
        rhoc = Mvir / quad(lambda r: 4 * np.pi * r**2 / ((r / rs)**a * (1. + r / rs)**(3 - a)),
                           1.e-8 * Rvir, Rvir, epsabs=0., epsrel=1.e-5)[0]
        return rhoc / (x**a * (1. + x)**(3 - a))

def errfunc_gnfw(p, r, logrho, Rvir, Mvir, w=1):
    c, a = p 
    fun = np.log10(rho_function(r, (c, a, Rvir, Mvir), model='gnfw'))
    return w * (fun - logrho)

def errfunc_nfw(p, r, logrho, Rvir, Mvir, w=1):
    c = p[0]
    fun = np.log10(rho_function(r, (c, Rvir, Mvir), model='nfw')) 
    return w * (fun - logrho)

def sigma_function(y, ymodel):
    residuals = y - ymodel
    N = np.size(y)
    return np.sqrt(np.sum(residuals**2) / N)
    
# NFW rho fit
def best_fit_NFW(r, rho, Rvir, Mvir):
    pfit, perr = optimize.leastsq(errfunc_nfw, (10.,), args=(r, np.log10(rho), Rvir, Mvir, 1.), full_output=0)
    c_nfw = pfit[0]
    params_nfw = (c_nfw, Rvir, Mvir)
    return params_nfw

# gNFW rho fit
def best_fit_gNFW(r, rho, Rvir, Mvir):
    pfit, perr = optimize.leastsq(errfunc_gnfw, (10., 1.), args=(r, np.log10(rho), Rvir, Mvir, 1.), full_output=0)
    c_gnfw, a_gnfw = pfit
    params_gnfw = (c_gnfw, a_gnfw, Rvir, Mvir)
    return params_gnfw
