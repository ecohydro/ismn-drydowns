#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Thu 04 Apr 24 11:15:18'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           sand_soil.py
Compatibility:  Python 3.7.0
Description:    Description of what program does

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
import numpy as np
import pandas as pd
import xarray as xr

import cartopy.crs as ccrs

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib import rcParams
from mpl_toolkits import mplot3d

%matplotlib qt

import seaborn as sns

import fluxtower
import utils
from drydowns.config import Config

#%% FUNCTIONS

proj_dir = '~/EcohydroRegimes'
data_dir = os.path.join(proj_dir, 'data')
gldas_dir = os.path.join(data_dir, 'GLDAS')

fn_frac = os.path.join(gldas_dir, 'GLDASp4_soilfraction_025d.nc4')
fn_porosity = os.path.join(gldas_dir, 'GLDASp5_porosity_025d.nc4')
fn_texture = os.path.join(gldas_dir, 'GLDASp5_soiltexture_025d.nc4')

#%%

texture = xr.open_dataset(fn_texture)
porosity = xr.open_dataset(fn_porosity)
soil_frac = xr.open_dataset(fn_frac)

soil = xr.Dataset({
    'texture_code': texture['GLDAS_soiltex'],
    'porosity': porosity['GLDAS_porosity'],
    'clay_fraction' : soil_frac['GLDAS_soilfraction_clay'],
    'sand_fraction' : soil_frac['GLDAS_soilfraction_sand'],
    'silt_fraction' : soil_frac['GLDAS_soilfraction_silt'],
})

gldas_codes = {
    1: 'sand',
    2: 'loamy sand',
    3: 'sandy loam',
    4: 'silt loam',
    5: 'silt',
    6: 'loam',
    7: 'sandy clay loam',
    8: 'silty clay loam',
    9: 'clay loam',
    10: 'sandy clay',
    11: 'silty clay',
    12: 'clay',
    13: 'organic materials',
    14: 'water',
    15: 'bedrock',
    16: 'other'
}


# %%

meta = fluxtower.flx_tower.META
# badm = fluxtower.flx_tower.FLX_BADM_DICT['DD'] 

# dup_sites = meta[meta.VARIABLE == 'LOCATION_LAT'][
#     meta[meta.VARIABLE == 'LOCATION_LAT'
# ].duplicated(subset='SITE_ID', keep=False)].SITE_ID.unique()


# flx_coords = meta[
#         (meta.VARIABLE.isin(['LOCATION_LAT','LOCATION_LONG']))
#         & (~meta.SITE_ID.isin(dup_sites)
#     ].pivot(
#             index='SITE_ID', columns='VARIABLE', values='DATAVALUE'
#         ).reset_index().rename(columns={'LOCATION_LAT':'LAT', 'LOCATION_LONG':'LON'})

# dup_sites = pd.DataFrame({
#     'SITE_ID' : meta[(meta.SITE_ID.isin(dup_sites)) & (meta.VARIABLE == 'LOCATION_LAT')].SITE_ID.reset_index(drop=True),
#     'LAT' : meta[(meta.SITE_ID.isin(dup_sites)) & (meta.VARIABLE == 'LOCATION_LAT')].DATAVALUE.reset_index(drop=True),
#     'LON' : meta[(meta.SITE_ID.isin(dup_sites)) & (meta.VARIABLE == 'LOCATION_LONG')].DATAVALUE.reset_index(drop=True),
# })


flx_coords = meta[
        (meta.VARIABLE.isin(['LOCATION_LAT','LOCATION_LONG']))
        # & (~meta.SITE_ID.isin(dup_sites)
    ].pivot_table(
            index='SITE_ID', columns='VARIABLE', values='DATAVALUE', aggfunc='first'
        ).reset_index().rename(columns={'LOCATION_LAT':'LAT', 'LOCATION_LONG':'LON'})

# def extract_xr_values(ds, lat, lon, variables=['texture', 'porosity']):
#     """
#     Get soil values from a dataset at a given lat, lon
#     """
#     values = {var : ds[var].sel(lat=lat, lon=lon, method='nearest').item() for var in variables}
#     return values

# def get_soil_values(df, ds, variables=['texture', 'porosity'], lat_col='LAT', lon_col='LON'):
#     """
#     Get soil values from a dataset at a given lat, lon
#     """
#     values = df.apply(
#         lambda row: pd.Series(extract_xr_values(ds, row[lat_col], row[lon_col], variables)), axis=1
#     )
#     return values

# Extract texture and porosity values into separate columns
# flx_coords[['texture', 'porosity']] = flx_coords.apply(
#     lambda row: pd.Series(get_soil_values(soil, row['LAT'], row['LON'])), axis=1
# )
flx_coords[list(soil.data_vars)] = utils.get_soil_values(flx_coords, soil, variables=list(soil.data_vars))
flx_coords.texture_code = flx_coords.texture_code.round(0)
flx_coords.insert(4, 'texture', flx_coords.texture_code.map(gldas_codes))

# flx_coords.to_csv(os.path.join(data_dir, 'fluxnet_soilproperties.csv'), index=False)
# flx_coords.to_csv('/Users/brynmorgan/dev/ecflux/src/fluxtower/metadata/soil_properties.csv', index=False)

#%% SMAP GRID

config = Config('soil/config.ini')
# cfg = config[config.get('RUN','type')]
proj_dir = config.config.get('ISMN', 'project_dir')
data_dir = os.path.join(proj_dir, 'data')
out_dir = os.path.join(proj_dir, 'outputs')

chirps_dir = os.path.join(out_dir, 'chirps')
grid_dir = os.path.join(data_dir, 'grids')

grid, lats, lons = utils.create_grid(
    grid_dir, 
    lat_file='coord_info_unique_column.csv', 
    lon_file='coord_info_unique_row.csv',
    return_dfs=True
)


gldas_smap = utils.resample_res(soil, grid)

gldas_smap.flatten()