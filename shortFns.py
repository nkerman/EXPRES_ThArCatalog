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
num_thids = 701

# %%
def Chauv(N_):
     #Implemenetation of Chauvinet's Principle
     sig_cutoff = np.sqrt(2) * erfinv((1 - (1/(2*N_))))
     return( sig_cutoff) #this should be the sigma cuttoff from the median

# %%
def countReal(array_):
    try:
        return np.count_nonzero(~np.isnan(array_))
    except:
        return -1
# %%
def removeOutliersChauv(df, line, verbose = True): #Removing outliers via a Statistical Chauvinet Criterion and self-defined sigma for a line
    lineNum = (list(df.columns.values)).index(line)
    array_ = np.array(df[line][0:num_thids])
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
            dfm[line].iat[i+num_thids] = np.NaN
    return(numExist, curr_Chauv_cutoff, medVal, deltaVal,minAccept, maxAccept, array_, array_orig)
# %%
def removeOutliers_setSigma(df_, line,  setSigma, centerMethod = 'median', verbose = False, oldDf =None, numSigma = 3): #Removing outliers via a number of a specific numerical value
    unfilteredLine = (df_[line][:num_thids][df_[line][:num_thids].notnull()])
    init_NExist = (np.size(df_[line][:num_thids][df_[line][:num_thids].notnull()]))
    # unfilteredLine

    if centerMethod == 'median':
        center = unfilteredLine.median()
    if centerMethod == 'mean':
        center = unfilteredLine.mean()
    if centerMethod not in ['mean','median']:
        print("Unacceptable method")
        center = 1
    condition = abs(unfilteredLine - center) > numSigma*setSigma
    invcondition = abs(unfilteredLine - center) < numSigma*setSigma
    df_[line][:num_thids].mask(condition, inplace = True)
    filteredLine = (df_[line][:num_thids][df_[line][:num_thids].notnull()])
    end_NExist = (np.size(df_[line][:num_thids][df_[line][:num_thids].notnull()]))
    percentPassed = 100* end_NExist/init_NExist
    # print(percentPassed
    if verbose == True:
        print(centerMethod)
        a, b, c = plt.hist(oldDf[line][:num_thids],color = 'r', alpha = 1,bins = 30 ,label = "Before $\sigma$-cut") #plt.cm.rainbow(0.15) looks nice
        plt.hist(df_[line][:num_thids], bins = b, color= 'k' , alpha = 0.9 , label = f"After {numSigma} $\sigma$-cut")
        plt.legend()
        plt.title(f"Set $\sigma$-clipping\nLine = {line} $\AA$ | $\sigma$ = {setSigma} | $N_{{\sigma}}$ = {numSigma}")

    return unfilteredLine, filteredLine, percentPassed
# %%
t0 = time.time()
def tester():
    return time.time() - t0
tester()
# %%
def countN(df, all_lines):
    countNaN = []
    for i, line in enumerate(all_lines):
        numNans = 0
        col = df[line]
        for j, val in enumerate(col[:num_thids]):
            if np.isnan(val):
                numNans += 1
        countNaN.append(numNans)
    countExist = np.subtract(num_thids, countNaN)
    rtN = np.sqrt(countExist)
    return rtN, countExist, countNaN
# %%
def readCat(pathToFile, catalog_file = 'thid_data_II.xdr', plot_ = False):
        # Gets the ThAr catalog as th_wvac, th_inten, th_wnum
    thidDataFile = idl.readsav(pathToFile + catalog_file)
    th_inten = thidDataFile['th_inten']
    th_wvac = thidDataFile['th_wvac']
    th_wnum = thidDataFile['th_wnum']
    if plot_ == True:
        exec("plt.plot(%s)"%key)
        plt.title(key)
        plt.xlabel("DataPoint Number")
        plt.ylabel(key)
        plt.show() #Plot catalog stuff if you want
    return th_wvac, th_inten, th_wnum

# # %%
def loadLFC(checkpointfile = '/Users/nkerman/OneDrive - Yale University/Spring_2020/Thesis/Compare_LFC_ThAr/singleData/checkpoint/LFC_191115.1144.npy', returnAll = False):
    checkpointdata = np.load(checkpointfile, allow_pickle = True)[()]
    lfc_params= checkpointdata["params"]
    lfc_fwhm= checkpointdata["fwhm"]
    lfc_wvln= checkpointdata["wvln"]
    lfc_cov= checkpointdata["cov"]
    lfc_quality= checkpointdata["quality"]
    lfc_pixels_allOrders = dict.fromkeys(range(0,160))
    for relOrder in range(np.size(lfc_params)):
        if np.size(lfc_params[relOrder]) > 1:
            lfc_pixels = lfc_params[relOrder][:,1]
            lfc_pixels_allOrders[relOrder] = lfc_pixels
    if returnAll == True:
        return lfc_pixels_allOrders, lfc_wvln, lfc_fwhm, lfc_params, lfc_cov,lfc_quality
    else:
        return lfc_pixels_allOrders, lfc_wvln, lfc_fwhm
# %%

# %%
