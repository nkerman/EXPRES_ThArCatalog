import shortFns as sf
import pandas as pd
import numpy as np
# %%
inputsDir = '/Users/nkerman/OneDrive - Yale University/Spring_2020/Thesis/BuildCatalog/'
plotsDir = '/Users/nkerman/OneDrive - Yale University/Spring_2020/Thesis/BuildCatalog/plots/'
lineHistDir = plotsDir+ 'lineHist/'
errorbarDir = plotsDir+ 'errorBars_I/'

# %%
reloaded_total_df = pd.read_pickle(inputsDir +'total_dataFrame_AsOf200213.pkl')
df = reloaded_total_df.copy()
del reloaded_total_df
# all_lines =list(df.columns.values)[2:]
df.reset_index(inplace=True)
df.replace(to_replace=[None], value=np.nan, inplace=True)
df_ = df.copy()
# %%
all_lines = df.columns.values[3:]
sf.countN(df, all_lines)
# %%
df
all_lines
