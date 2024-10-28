######------ Calculation growth factor------######
import numpy as np
from scipy.integrate import solve_ivp
from constants import h, omega_m

def Hz(a, omega_m, h):
    H0 = h * 100.0 / 3.08567758e19  # Hubble constant in 1/s
    return H0 * np.sqrt(omega_m / a**3 + (1.0 - omega_m))

def fcn(a, y):
    y1, y2 = y
    Hz_value = Hz(a, omega_m, h)
    dy1_da = y2 / (a * Hz_value)
    dy2_da = -2.0 * (y2 / a) + 1.5 * omega_m * (h * 100.0 / 3.08567758e19)**2 * y1 / (Hz_value * a**4)
    return [dy1_da, dy2_da]

def linear_growth(z, tol=1e-6):
    a_start = 1.0e-5
    a_final = 1.0 / (1.0 + z)
    
    # Initial conditions
    y0 = [a_start, 0.0]
    
    # Solve ODE from a_start to 1 (normalization)
    sol_norm = solve_ivp(fcn, [a_start, 1.0], y0, rtol=tol, atol=tol)
    norm = sol_norm.y[0, -1]
    
    # Solve ODE from a_start to a_final
    sol = solve_ivp(fcn, [a_start, a_final], y0, rtol=tol, atol=tol)
    
    delta = sol.y[0, -1] / norm
    ddot = sol.y[1, -1] / norm
    
    return delta, ddot
