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

# Default Plot
def plot_Variogram(list_Distance, list_Bin, list_Fitted_Model, list_Range, list_Sill, list_Nugget, list_Empiral_X,list_Empiral_Y, path_SaveFig, bool_SaveFig = False, bool_CloseFig = False, site_Name_Graph = "Site", figsize = (20,10), my_Colors = ["Grey","Red","Black","Orange","Green","Blue"], my_Markers = ["2","o","v","H","s","*"]):
    plt.figure(figsize = figsize)
    for i in range(len(list_Distance)):
        temp_x = crea_Linspace_X(list_Bin[i])
        if len(list_Bin[i]) >= 4:
            temp_cubic_model = crea_Cubic_Model(list_Fitted_Model[i], list_Bin[i])
            plt.plot(temp_x, temp_cubic_model(temp_x), label = f'{list_Distance[i]} Theoretical; a = {list_Range[i]:.2F}; c = {list_Sill[i]:.2F}; c0 = {list_Nugget[i]:.2E}', color = my_Colors[i], linestyle = '-')
        else:
            plt.plot(temp_x, list_Fitted_Model[i](temp_x), label = f'{list_Distance[i]} Theoretical; a = {list_Range[i]:.2F}; c = {list_Sill[i]:.2F}; c0 = {list_Nugget[i]:.2E}', color = my_Colors[i], linestyle = '-')
    for i in range(len(list_Distance)):
        plt.scatter(list_Empiral_X[i],list_Empiral_Y[i], color = my_Colors[i], label = f'{list_Distance[i]} Experimental', marker = my_Markers[i])
    plt.title("Variogram of " + site_Name_Graph)
    plt.xlabel("Distance or Lag")
    plt.ylabel("Semivariance")
    plt.legend()
    plt.grid()
    if bool_SaveFig:
        plt.savefig(path_SaveFig)
    if bool_CloseFig:
        plt.close()


# # Single plot
# def plot_Single_Variogram(distance, bin, fitted_model, range, sill, nugget, figsize = (20,10)):
#     plt.figure(figsize = figsize)
#     temp_x = crea_Linspace_X(bin)
#     plt.plot(temp_x, crea_Cubic_Model(temp_x), label = f'30m Theoretical; a = {range:.2F}; c = {sill:.2F}; c0 = {nugget:.2E}', color = 'Grey', linestyle = '-')