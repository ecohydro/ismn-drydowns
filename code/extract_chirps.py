#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Thu 18 Apr 24 18:22:14'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           extract_chirps.py
Compatibility:  Python 3.12.0

Description:   Extract CHIRPS precipitation data for ISMN stations.

URL:            https://

Requires:       list of libraries required

Dev ToDo:       None

AUTHOR:         Bryn Morgan
ORGANIZATION:   University of California, Santa Barbara
Contact:        bryn.morgan@geog.ucsb.edu
Copyright:      (c) Bryn Morgan 2024


"""


# IMPORTS
import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import time
import multiprocessing as mp
import ast

#import drydowns
from drydowns.config import Config


import utils

# VARIABLES

config = Config('soil/config.ini')
# cfg = config[config.get('RUN','type')]
proj_dir = config.config.get('ISMN', 'project_dir')
data_dir = os.path.join(proj_dir, 'data')


chirps_dir = '/home/chc-data-out/products/CHIRPS-2.0/global_daily/netcdf/'

# grid_dir = os.path.join(data_dir, 'grids')

# chirps_out_dir = os.path.join(proj_dir, 'outputs', 'chirps')
out_dir = os.path.join(data_dir, 'ISMN', 'ancillary', 'CHIRPS')

anc_file = os.path.join(data_dir, 'ISMN', 'ismn_ancillary.csv')

res = '05'
# Use .25° for testing; will use .05° for final analysis
fp_chirps = os.path.join(chirps_dir, f'p{res}')
files_yr = sorted([f for f in glob.glob(fp_chirps + '/*.nc') if int(f[-16:-12]) < 2024])



# FUNCTIONS

def extract_xr_values(
    df, ds_in, 
    subset_args = {'dim':'location', 'lat_col':'latitude', 'lon_col':'longitude'},
    rolling_args = {'window':60, 'min_periods' : 60, 'var':{'total60d':'precip'}},
    return_ds = True
):
    # dim='location', lat_col='latitude', lon_col='longitude', 
    # ):
    # Subset the dataset
    ds = utils.subset_xr(df, ds_in, **subset_args)
    # Calculate rolling sum
    ds = calculate_rolling_sum(ds, **rolling_args)
    # Extract values to dataframe
    print(f"Extracting values from {len(df)} stations")
    vals = ds.to_dataframe().reset_index()

    print(f"Extracted {len(vals)/len(df)} values each from {len(df)} stations")
    if return_ds:
        return vals, ds
    else:
        return vals

# def subset_xr(df, ds, dim='location', lat_col='latitude', lon_col='longitude'):
#     print(f"Subsetting dataset with {len(df)} stations")
#     lats = xr.DataArray(df[lat_col].values, dims=dim, coords={dim: df[dim]})
#     lons = xr.DataArray(df[lon_col].values, dims=dim, coords={dim: df[dim]})
    
#     coords = dict(zip(utils.get_coord_names(ds), [lats, lons]))
    
#     ds_sub = ds.sel(coords, method='nearest')

#     return ds_sub

def calculate_rolling_sum(ds, window=60, min_periods=1, var={'total60d':'precip'}):
    print(f"Calculating rolling sum over {window} days")
    roll = ds.rolling(time=window, min_periods=min_periods).sum()
    var_map = {k : roll[v] for k,v in var.items()}
    ds = ds.assign(**var_map)

    return ds

# def expand_col(df, col='location', cols=['network','station']):
#     in_cols = list(df.columns[2:])
#     df[col] = df[col].apply(ast.literal_eval)
#     df[cols] = pd.DataFrame(df[col].to_list(), index=df.index)
#     # Reorder columns
#     df = df[['time']+cols+in_cols]
#     return df


# def agg_file(
#     file, col='precip', roll_args={'window':60, 'min_periods':1}, 
#     names=['total60d','total_month']
# ):
#     df = pd.read_csv(file)
#     if all([n in df.columns for n in names]):
#         print(f"File {os.path.basename(file)} already aggregated. Skipping.")
#         # return df
#         return

#     df = utils.expand_col(df)

#     try:
#         df = utils.agg_data(df, col=col, roll_args=roll_args, names=names)
#     except Exception as e:
#         print(f"Error aggregating file {file}: {e}")
#     df.to_csv(file, index=False)
#     print(f"Saved {os.path.basename(file)}"
# )



#%MAIN
def main():
    start = time.perf_counter()
    # Load rainfall data
    ds_all = xr.open_mfdataset(files_yr)

    if os.path.exists(anc_file):
        ismn = pd.read_csv(anc_file)
    else:
        from soil.station import META
        ismn = META[['network','station','latitude','longitude']].drop_duplicates()

    ismn['location'] = ismn.apply(lambda x: tuple(x[['network','station']]), axis=1)

    # Extract values
    vals,ds = extract_xr_values(ismn, ds_all)
    
    m1 = time.perf_counter()
    print(f"Finished extracting in {m1-start:.6f} seconds.")
    # Make sure out_dir exists
    utils.create_dir(out_dir)
    # Save to netcdf
    # nc_out = f'ismn_chirps{res}.nc'
    # print(f'Saving dataset to {nc_out}')
    # ds.to_netcdf(os.path.join(out_dir, nc_out))

    m2 = time.perf_counter()
    # print(f"Finished saving netcdf in {m2-m1:.6f} seconds.")

    ds_all.close()


    # Save to file
    utils.separate_df(
        df=vals, by=ismn,#['network','station'], 
        suff=f'_chirps{res}', out_dir=out_dir,
        sub_cols = ['location']
    )

    end = time.perf_counter()
    print(f"Finished separation in {end-m2:.6f} seconds.")
    print(f"Finished processing in {end-start:.6f} seconds.")

    # NOTE: code below was for running alone to add columns to existing files.
    # It needs to be incorporated above if running again.
    # TODO: Don't overwrite existing files.
    
    # files = sorted([f for f in glob.glob(os.path.join(out_dir, '*.csv'))])
    # roll_args = {'window': 60, 'min_periods': 60}
    # names = ['total60d', 'total_month']
    # recalc = True
    # expand_kwargs = {'col':'location', 'cols':['network','station']}

    # with mp.Pool(50) as pool:
    #     pool.starmap(
    #         utils.agg_file, 
    #         [(f, 'precip', roll_args, names, recalc, expand_kwargs) for f in files]
    #     )
    # pool.close()
    # pool.join()

#%%
if __name__ == "__main__":
    main()