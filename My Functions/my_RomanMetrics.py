import math
from scipy.interpolate import interp1d
from scipy.integrate import quad

# Rcv metric
def calc_Rcv(CV_Ext, CV_Int):
    return (CV_Ext - CV_Int) / CV_Int

# RAW Score
def calc_RAWScore(Rcv):
    return (abs(2 * Rcv)) ** (-1)

# Rse metric
def calc_Rse(g, range_Ext, range_Int):
    return math.exp(-((g / range_Int) ** 2 + (g / range_Ext) ** 2) ** (1/2))

# ST
def calc_ST(semivar_range_Exp, nugget):
    return (semivar_range_Exp - nugget) / semivar_range_Exp

# Rst metric
def calc_Rst(ST_Ext, ST_Int):
    return (ST_Ext - ST_Int) / ST_Int

# SV
def calc_SV(fitted_model, range_, nugget, sill):
    def f_SV(x):
        return (fitted_model(x) - nugget) / sill
    return quad(f_SV, 0, range_)[0]

# Rsv metric
def calc_Rsv(SV_Ext, SV_Int):
    return (SV_Ext - SV_Int) / SV_Int

# ST Score
def calc_STScore(Rcv, Rse, Rst, Rsv):
    return (abs(Rcv) / 3 + abs(Rst) / 3 + abs(Rsv) / 3 +Rse) ** (-1)