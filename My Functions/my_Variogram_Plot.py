import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rasterio as rio
import skgstat as skg
from scipy.interpolate import interp1d
from scipy.integrate import quad

# Create a cubic model to plot curve
def crea_Cubic_Model(fitted_model, bin, interp1d_kind = "cubic"):
    def fx(x):
        return fitted_model(x)
    cubic_model = interp1d(np.concatenate((np.array([0]),bin)), fx(np.concatenate((np.array([0]),bin)),), kind = interp1d_kind)
    return cubic_model

# Create linspace continuous x coordinates from bin
def crea_Linspace_X(bin):
    return np.linspace(0,bin[-1])

# Single plot
def plot_Single_Variogram(distance, bin, fitted_model, range, sill, nugget, figsize = (20,10)):
    plt.figure(figsize = figsize)
    temp_x = crea_Linspace_X(bin)
    plt.plot(temp_x, crea_Cubic_Model(temp_x), label = f'30m Theoretical; a = {range:.2F}; c = {sill:.2F}; c0 = {nugget:.2E}', color = 'Grey', linestyle = '-')