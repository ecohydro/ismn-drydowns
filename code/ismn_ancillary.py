#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Wed 10 Apr 24 17:17:04'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           sand_ismn_rain.py
Compatibility:  Python 3.12.0

Description:    Create + update ISMN ancillary data.

URL:            https://

Requires:       list of libraries required

Dev ToDo:       None

AUTHOR:         Bryn Morgan
ORGANIZATION:   University of California, Santa Barbara
Contact:        bryn.morgan@geog.ucsb.edu
Copyright:      (c) Bryn Morgan 2024


"""


#%% IMPORTS

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr


import drydowns
from drydowns.config import Config

# import soil.station as station
from soil.station import META

import utils

#%% FUNCTIONS

config = Config('soil/config.ini')
# cfg = config[config.get('RUN','type')]
proj_dir = config.config.get('ISMN', 'project_dir')
data_dir = os.path.join(proj_dir, 'data')
out_dir = os.path.join(proj_dir, 'outputs')

chirps_dir = os.path.join(out_dir, 'chirps')

grid_dir = os.path.join(data_dir, 'grids')

anc_file = os.path.join(data_dir, 'ISMN', 'ismn_ancillary.csv')

if os.path.exists(anc_file):
    ismn = pd.read_csv(anc_file)
else:
    ismn = META[['network','station','latitude','longitude']].drop_duplicates()

ismn['location'] = ismn.apply(lambda x: tuple(x[['network','station']]), axis=1)



#%% TIMERANGE

def get_timerange(df, start_col='timerange_from', end_col='timerange_to'):
    # Start
    start = META.groupby(
        by=['network','station']
    )[start_col].min().dt.year.reset_index()
    # End
    end = META.groupby(
        by=['network','station']
    )[end_col].max().dt.year.reset_index()

    df = df.merge(start, on=['network','station'], how='left')
    df = df.merge(end, on=['network','station'], how='left')

    df.rename(
        columns={'timerange_from': 'start_year', 'timerange_to': 'end_year'}, 
        inplace=True
    )

    return df

if 'start_year' not in ismn.columns:
    ismn = get_timerange(ismn)



#%% GRIDS

# Add grid coords to this dict for each dataset

coord_file_dict = {
    'smap' : ['coord_info_unique_row.csv', 'coord_info_unique_column.csv'],
    'chirps25' : ['chirps25_coord_info_row.csv', 'chirps25_coord_info_column.csv'],
    'chirps05' : ['chirps05_coord_info_row.csv', 'chirps05_coord_info_column.csv'],
}

# to_add = [name for name in coord_file_dict.keys() if f'{name}_lat' not in ismn.columns]

def get_data_grid(df, grid_dir, name='smap', lat_file=None, lon_file=None):
    if not f'{name}_lat' in df.columns:
        grid, lats, lons = get_data_coords(
            df, grid_dir, name=name, lat_file=lat_file, lon_file=lon_file
        )

    if not f'{name}_row_index' and not 'EASE_row_index' in df.columns:
        df = get_grid_indices(df, grid_dir, name=name, lats=lats, lons=lons)

    return df

def get_data_coords(df, grid_dir, name='smap', lat_file=None, lon_file=None):
    if not lat_file:
        lat_file, lon_file = coord_file_dict[name]

    grid, lats, lons = utils.create_grid(
        grid_dir, 
        lat_file=lat_file, 
        lon_file=lon_file,
        return_dfs=True
    )
    # Get the grid coordinates + add to df
    df[[f'{name}_lat',f'{name}_lon']] = utils.get_xr_values(
        df, grid, variables=['latitude','longitude'], 
        lat_col='latitude', lon_col='longitude', 
        func=utils.extract_xr_value
    ).rename(columns={'latitude':f'{name}_lat','longitude':f'{name}_lon'})
    
    print(f"Added {name} coordinates to df")

    return grid, lats, lons

def get_grid_indices(df, name='smap', lats=None, lons=None):
    if not lats or not lons:
        grid, lats, lons = get_data_coords(
            df, grid_dir, name=name
        )
    # Get the grid row and column indices
    for dfi, var in zip([lats, lons], ['latitude','longitude']):
        dfi.rename(
            columns={col:f'{name}_{col}' for col in ['row_index','column_index']}, 
            inplace=True
        )
        df = df.merge(
            dfi, left_on=f'{name}_{var[:3]}', right_on=var, how='left', suffixes=('','_y')
        )
        del df[var+'_y']
    
    print(f"Added {name} coordinates to df")
    
    return df


for name in coord_file_dict.keys():
    ismn = get_data_grid(ismn, grid_dir, name=name)



#%% MAP

# TODO: Move this function to somewhere if it's generalizable
def get_data(df, file, col_map={'precip' : 'MAP', 'pet':'MEAN_PET'}):
    ds = xr.open_dataset(file)
        # Extract values + add to df
    df[list(ds.data_vars)] = utils.get_xr_values(
        df, ds, variables=list(ds.data_vars), 
        lat_col='latitude', lon_col='longitude',
        func = utils.extract_xr_value
    )
    # Rename columns
    df.rename(columns=col_map, inplace=True)
    return df

def get_map(df, res='05', fdate_map='12mar'):
    map_file = utils.get_filename('map', f'p{res}', fdate_map)
    # Open file
    map_ds = utils.open_file(
        map_file, chirps_dir, resample=False
    )
    # Extract values + add to df
    df[list(map_ds.data_vars)] = utils.get_xr_values(
        df, map_ds, variables=list(map_ds.data_vars), 
        lat_col='latitude', lon_col='longitude',
        func = utils.extract_xr_value
    )
    # df.rename(columns={'precip':f'MAP_{res}'}, inplace=True)
    df.rename(columns={'precip': 'MAP'}, inplace=True)
    return df

res = '05' 
fdate_map = '12mar'

# if f'MAP_{res}' not in ismn.columns:
if 'MAP' not in ismn.columns:
    # ismn = get_map(ismn, res=res, fdate_map=fdate_map)
    map_file = utils.get_filename('map', f'p{res}', fdate_map)
    ismn = get_data(ismn, map_file, col_map={'precip':'MAP'})



#%% MEAN PET  + AI


# def get_pet(df, file):
#     # Open file
#     ds_in = xr.open_dataset(file)
    # # Subset the dataset
    # ds = utils.subset_xr(
    #     df, ds_in, dim='location', lat_col='latitude', lon_col='longitude'
    # )
    # # Extract values to dataframe
    # vals = ds.to_dataframe().reset_index()
    # # Add network + station cols and remove location
    # vals = utils.add_network_station(vals)
    # Add to df
    # df = df.merge(vals, on=['network','station'], how='left')
    # return df

if 'MEAN_PET' not in ismn.columns:
    pet_file = os.path.join(out_dir, 'pet', 'mean_annual_pet.nc')
    ismn = get_data(ismn, pet_file, col_map={'pet':'MEAN_PET'})

# AI
if 'AI' not in ismn.columns:
    ismn['AI'] = ismn['MAP'] / ismn['MEAN_PET']




#%%

# ADD STUFF HERE AS IT IS USEFUL




#%% SAVE

save = True

if save:
    ismn.to_csv(
        os.path.join(data_dir, 'ISMN', 'ismn_ancillary.csv'), index=False
    )


#%%

# Chunks of ISMN coords for APPEEARS requests.
# for i in range(6):
#     n0 = i*515
#     n = (i+1)*515
#     if i == 5:
#         n += 5
    
#     df = ismn_coords.iloc[n0:n]
#     df.to_csv(
#         os.path.join(config.config.get('ISMN','data_dir'),'..',f'ismn_coords_0{i}.csv'), index=False
#     )
