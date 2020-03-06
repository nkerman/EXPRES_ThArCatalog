import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io.idl as idl
from scipy.io import readsav
from astropy.io import fits
from scipy.interpolate import LSQUnivariateSpline
import pandas as pd
import pickle
import glob
import time
from scipy.special import erfinv as erfinv
# %%
def Chauv(N_):
     #Implemenetation of Chauvinet's Principle
     sig_cutoff = np.sqrt(2) * erfinv((1 - (1/(2*N_))))
     return( sig_cutoff) #this should be the sigma cuttoff from the median

# %%

def removeOutliers(df, line, verbose = True):
    lineNum = (list(df.columns.values)).index(line)
    array_ = np.array(df[line][0:701])
    array_orig = array_.copy()
    numExist, curr_Chauv_cutoff = countExist[lineNum-3], Chauv(countExist[lineNum-3])
    medVal = np.nanmedian(array_)
    deltaVal = (np.nanstd(array_) * curr_Chauv_cutoff)
    minAccept, maxAccept = medVal - deltaVal, medVal + deltaVal
    for i, val in enumerate(array_):
        if np.abs(np.subtract(val, medVal)) > deltaVal:
            if verbose == True:
                print(i, val, val - medVal)
            array_[i] = np.NaN
            dfm[line].iat[i] = np.NaN
            dfm[line].iat[i+701] = np.NaN
    return(numExist, curr_Chauv_cutoff, medVal, deltaVal,minAccept, maxAccept, array_, array_orig)
# %%
def countN(df, all_lines):
    countNaN = []
    for i, line in enumerate(all_lines):
        numNans = 0
        col = df[line]
        for j, val in enumerate(col[:701]):
            if np.isnan(val):
                numNans += 1
        countNaN.append(numNans)
    countExist = np.subtract(701, countNaN)
    rtN = np.sqrt(countExist)
# %%
