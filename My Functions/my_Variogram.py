# A python script file used to output 

import numpy as np
import skgstat as skg

# The function to retrieve all the geographical coordinates from a rasterio img input
def get_Coordinates_from_Rio(img):
    temp_list_Coordinates = []
    for x in range(img.shape[0]):
        for y in range(img.shape[1]):
            temp_x, temp_y = img.xy(x, y)
            temp_list = [temp_x, temp_y]
            temp_list_Coordinates.append(temp_list)
    temp_list_Coordinates = np.array(temp_list_Coordinates)
    return(temp_list_Coordinates)

# Calculate the bin required to calculate variogram and then roman metrics
def calc_Bin(distance, resolution = 10, step = 10):
    if step < resolution:
        print("Error: The step must be a multiple of the resolution! ")
        return
    else:
        temp_Diagonal = distance * (2 ** 0.5) / 2
        temp_Bin = np.arange(resolution, temp_Diagonal, step)  
        return(temp_Bin)

# Calculate variogram and return its parameters, bins and fitted models, which will be used
def calc_Variogram_Parameters(band, coordinates, distance, resolution = 10, step = 10, nugget_Bool = True):
    if len(band.shape) != 1:
        print("Please reshape the band to one dimension! ")
        return
    else: 
        temp_Band_Reshape = band
        temp_Coorindates = coordinates
        temp_Bin_Func = calc_Bin(distance, resolution, step)
        temp_Variogram = skg.Variogram(temp_Coorindates, temp_Band_Reshape, use_nugget=nugget_Bool, bin_func=temp_Bin_Func, maxlag = temp_Bin_Func[-1])
        temp_Range = temp_Variogram.parameters[0]
        temp_Nugget = temp_Variogram.parameters[2]
        temp_Sill = temp_Variogram.parameters[1] + temp_Nugget
        temp_Bins = temp_Variogram.bins
        temp_Experimental = temp_Variogram.experimental
        temp_FittedModel = temp_Variogram.fitted_model
        return temp_Range, temp_Nugget, temp_Sill, temp_Bins, temp_Experimental, temp_FittedModel

# Get the semivariance of the distance equal to range in the experimental model
def get_Semivar_Exp(distance, bins, experimental_model):
    for i in range(len(bins)):
        if i + 1 == len(bins):
            index = i
        elif distance >= bins[i] and distance <= bins[i+1]:
            index = i + 1
            break
    temp_Semivar_Exp = experimental_model[index]
    return temp_Semivar_Exp

def calc_STD(band):
    return np.std(band)

def calc_CV(band):
    return np.std(band) / np.mean(band)