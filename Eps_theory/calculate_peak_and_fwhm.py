import numpy as np
from scipy.signal import find_peaks

def calculate_peak_and_fwhm(z_values, dP_dz):
    """
    Calculate the peak position and the points where the curve crosses the half-maximum.

    Parameters:
    z_values (numpy array): Independent variable values (x-axis)
    dP_dz (numpy array): Target curve (or its derivative)

    Returns:
    z_peak_value (float): The position of the peak (x-coordinate)
    left_z (float): The x-coordinate where the curve crosses half-maximum on the left
    right_z (float): The x-coordinate where the curve crosses half-maximum on the right
    """
    # Find peaks in the curve (we assume we're looking for the highest peak)
    peaks, _ = find_peaks(dP_dz, height=0.2)
    if len(peaks) == 0:
        raise ValueError("No peaks found, please check the input data.")

    # Find the index of the highest peak
    max_peak_index = peaks[np.argmax(dP_dz[peaks])]
    z_peak = z_values[max_peak_index]

    # Calculate the half-maximum value (i.e., the value at half the peak height)
    half_max_value = dP_dz[max_peak_index] / 2.0

    # Find the indices where the curve crosses the half-maximum on both sides of the peak
    left_indices = np.where(dP_dz[:max_peak_index] < half_max_value)[0]
    right_indices = np.where(dP_dz[max_peak_index:] < half_max_value)[0] + max_peak_index

    if len(left_indices) == 0 or len(right_indices) == 0:
        raise ValueError("Unable to find the positions where the curve crosses the half-maximum, please check the peak shape.")

    # Find the corresponding z-values for the left and right half-maximum crossings
    left_z = z_values[left_indices[-1]]
    right_z = z_values[right_indices[0]]
    # FWHM = right_z - left_z
    return z_peak, left_z, right_z
