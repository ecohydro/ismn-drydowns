#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Mon 15 Apr 24 00:20:07'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           extract_pet.py
Compatibility:  Python 3.12.0

Description:    Extract dPET data for ISMN stations.

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


# import drydowns
from drydowns.config import Config

# import soil.station as station
# from soil.station import Station, META, ismn_data

import utils


#% DIRECTORIES

config = Config('soil/config.ini')
# cfg = config[config.get('RUN','type')]
proj_dir = config.config.get('ISMN', 'project_dir')
data_dir = os.path.join(proj_dir, 'data', 'ISMN')
# out_dir = os.path.join(proj_dir, 'outputs')
out_dir = os.path.join(data_dir, 'dPET', 'byyear')

chirps_dir = os.path.join(out_dir, 'chirps')
grid_dir = os.path.join(data_dir, 'grids')

pet_dir = os.path.join(config.config.get('SMAP','project_dir'), 'data', 'PET')

# pet_files = [f for f in glob.glob(os.path.join(pet_dir, '*_daily_pet.nc'))]
# pet_files.sort()


# ismn_coords = META[['network','station','latitude','longitude']].drop_duplicates()
ismn_coords = META.groupby(
    by=['network','station','latitude','longitude']
).timerange_from.min().dt.year.reset_index()
ismn_coords.rename(columns={'timerange_from':'start_year'}, inplace=True)

def process_coords(ds, coords, file, log=False):
    
    if log:
        print(
            f"Extracting values at {coords['latitude']}, {coords['longitude']}",
            f"for file {file}"
        )
    pet_vals = utils.extract_xr_values(
        ds, coords['latitude'], coords['longitude']
    )
    return pet_vals


def process_file(file, df):
    pet = xr.open_dataset(file)

    # pet_dfs = []
    pet_dfs = {}

    for n,(i, row) in enumerate(df.iterrows()):
        stat = row.to_dict()
        pet_vals = process_coords(pet, stat, file)
        pet_vals.insert(0, 'network', stat['network'])
        pet_vals.insert(1, 'station', stat['station'])
        pet_dfs[(stat['network'],stat['station'])] = pet_vals
        # pet_dfs.append(pet_vals)
        # Print update every 200 stations
        if (n+1) % 10 == 0:
            print(f"Processed {n+1}/{len(df)} stations")

    print(f"Finished processing {file}")
    
    pet.close()

    return pet_dfs

        # print(la)
        # break


# Separate out skipped dfs into individual station dfs
def separate_yr_df(file):
    print(f"Reading in {file}")
    df = pd.read_csv(os.path.join(out_dir, file))#, index_col=0)
    df['time'] = pd.to_datetime(df['time'])
    df.set_index('time', inplace=True)

    # Get unique stations
    ustations = df[['network','station']].drop_duplicates()
    pet_dfs = {}
    # Iterate through unique stations + extract station data
    for i,row in ustations.iterrows():
        stat = row.to_dict()
        df_yr = df[(df.network == stat['network']) & (df.station == stat['station'])]
        pet_dfs[(stat['network'],stat['station'])] = df_yr
    return pet_dfs

# Groupby start year and get list of coords for each year
# yr_dict = {}
# for yr in range(1981,2024):
#     # if yr == 1981:
#     #     df = ismn_coords[ismn_coords.start_year <= yr]
#     # else:
#     #     df = ismn_coords[ismn_coords.start_year == yr]
#     df = ismn_coords[ismn_coords.start_year <= yr]
#     yr_dict[yr] = df
# yrs = ismn_coords.groupby('start_year').apply(
#     lambda grp : grp.to_dict(orient='records'),
#     include_groups=False
# )


    # utils.extract_xr_values(pet, la)

def fix_2007(yr, file):

    df = ismn_coords[ismn_coords.start_year <= yr]
    df_yr = pd.read_csv(os.path.join(out_dir, 'ismn_pet_2007_01.csv'))
    ns_all = df[['network','station']].apply(tuple, axis=1)
    ns_done = df_yr[['network','station']].apply(tuple, axis=1).unique()
    i_todo = ns_all[~ns_all.isin(ns_done)].index
    df = df.loc[i_todo]
    print(f"Processing {file}")
    try:
        pet_dfs = process_file(file, df)
    except Exception as e:
        print(f"Error processing {file}: {e}")
    df_yr2 = pd.concat(pet_dfs.values())
    df_yr = pd.concat([df_yr, df_yr2])

    return df_yr

# MAIN
def main():
    pet_files = [f for f in glob.glob(os.path.join(pet_dir, '*_daily_pet.nc'))]
    pet_files.sort()

    pet_list = []
    # yr_list = []
    skip_list = []

    # for file,yr in zip(pet_files,range(1981,2024)):
    # for file in pet_files:
    for yr in range(1981,2024):
    # for yr in [2019, 2018, 2017, 2016, 2015, 2014, 2013, 2012, 2011, 2010]:
    # for yr in range(1999, 2002):
        file = [f for f in pet_files if str(yr) in f][0]
        # yr = int(pet_files[0].split('_')[-2])
        out_file = f'ismn_pet_{yr}.csv'

        # Check if file exists (has already been processed)
        if os.path.exists(os.path.join(out_dir, out_file)):
            print(f"{out_file} exists. Skipping for now; will read in later")
            skip_list.append(out_file)
            continue

        if yr == 2007:
            df_yr = fix_2007(yr, file)
            print(f"Saving {out_file}")
            df_yr.to_csv(os.path.join(out_dir, out_file))

            pet_list.append(pet_dfs)
            continue
        
        print(f"Processing {file}")
        df = ismn_coords[ismn_coords.start_year <= yr]
        try:
            pet_dfs = process_file(file, df)
        except Exception as e:
            print(f"Error processing {file}: {e}")
            continue
        df_yr = pd.concat(pet_dfs.values())

        # filename = f'ismn_pet_{yr}.csv'

        print(f"Saving {out_file}")
        df_yr.to_csv(os.path.join(out_dir, out_file))

        pet_list.append(pet_dfs)
        # yr_list.append(yr)
    
    # Read in skipped files
    for file in skip_list:
        pet_dfs = separate_yr_df(file)
        pet_list.append(pet_dfs)
    

    for i,row in ismn_coords.iterrows():
        stat = row.to_dict()
        pets = [ dic[(stat['network'], stat['station'])] for dic in pet_list ]

        stat_df = pd.concat(pets.values())
        stat_df.sort_index(inplace=True)
        filename = f"{stat['network']}_{stat['station']}_pet.csv"
        print(f"Saving {filename}")
        stat_df.to_csv(os.path.join(
            #data_dir, 'dPET', filename
            out_dir, filename
        ))


#%%
if __name__ == "__main__":
    main()