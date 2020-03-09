import shortFns as sf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.io import readsav
from astropy.io import fits
# from scipy.interpolate import LSQUnivariateSpline
import pickle
# %%#These directories are the most likely things needing changing
inputsDir = '/Users/nkerman/OneDrive - Yale University/Spring_2020/Thesis/BuildCatalog/'
plotsDir = '/Users/nkerman/OneDrive - Yale University/Spring_2020/Thesis/BuildCatalog/plots/'
thidDatXdrDir ='/Users/nkerman/OneDrive - Yale University/Spring_2020/Thesis/ThAr_Solutions/data/cat/'
lineHistDir = plotsDir+ 'lineHist/'
errorbarDir = plotsDir+ 'errorBars_I/'

# %% ###PARAMS you need to check!
totalNumThids = 701

# %%
df = pd.read_pickle(inputsDir +'total_dataFrame_AsOf200213.pkl')
df.reset_index(inplace=True)
df.replace(to_replace=[None], value=np.nan, inplace=True)
df_ = df.copy()
# %%
all_lines = df.columns.values[3:]
rtN, countExist, countNaN = sf.countN(df, all_lines)
th_wvac, th_inten, th_wnum = sf.readCat(thidDatXdrDir)

# %%
# def singleResid(thidFile, LFCfiles, catFile)
lfc_pixels_allOrders, lfc_wvln, lfc_fwhm = sf.loadLFC() # This line loads the useful data for a single LFC file
