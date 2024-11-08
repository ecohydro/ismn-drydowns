#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Thu 18 Apr 24 16:23:17'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           separate_pet.py
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


# IMPORTS
import os
import glob
import re
# import numpy as np
import pandas as pd
# import xarray as xr

import multiprocessing as mp
# import drydowns
from drydowns.config import Config

import utils

# VARIABLES
config = Config('soil/config.ini')
# cfg = config[config.get('RUN','type')]
proj_dir = config.config.get('ISMN', 'project_dir')
data_dir = os.path.join(proj_dir, 'data')
# anc_dir = config.config.get('ISMN', 'anc_dir')
pet_dir = config.config.get('ISMN', 'pet_dir')

anc_file = os.path.join(data_dir, 'ISMN', 'ismn_ancillary.csv')

# FUNCTIONS
def create_dir(out_dir):
    # out_dir = os.path.join(parent_dir, subdir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"Created directory: {out_dir}")
    else:
        print(f"Directory already exists: {out_dir}")

def process_subset(row, df, suff, out_dir=None, parent_dir=None):
    info = row.to_dict()
    
    out_file = f"{info['network']}_{info['station']}{suff}.csv"
    if not out_dir:
        out_dir = os.path.join(parent_dir, info['network'])
        create_dir(out_dir)

    if os.path.exists(os.path.join(out_dir, out_file)):
        print(f"File {out_file} already exists. Skipping.")
    else:
        df_ns = df[
            (df['network'] == info['network']) 
            & (df['station'] == info['station'])
        ]
        df_ns.reset_index(drop=True, inplace=True)
        df_ns.to_csv(os.path.join(out_dir, out_file), index=False)
        print(f"Saved {out_file}")

def separate_df(df, by=['network', 'station'], threads=20, suff=None, out_dir=None, parent_dir=None):
    if isinstance(by, list):
        stations = df[by].drop_duplicates()
    elif isinstance(by, pd.DataFrame):
        stations = by
    
    if threads:
        with mp.Pool(threads) as pool:
            pool.starmap(
                process_subset, 
                [(row, df, suff, out_dir, parent_dir) for i,row in stations.iterrows()]
            )
        pool.close()
        pool.join()
    else:
        for i,row in stations.iterrows():
            process_subset(row, df, suff, out_dir, parent_dir)


def process_station(row, out_dir=pet_dir, parent_dir=os.path.join(pet_dir,'tmp')):
    if isinstance(row, tuple):
        row = row[1]
    network, station = row.to_dict().values()
    net_dir = os.path.join(parent_dir, network)
    file_pre = f"{network}_{station}"
    
    out_file = f"{file_pre}_dPET.csv"
    if os.path.exists(os.path.join(out_dir, out_file)):
        print(f"File {out_file} already exists. Skipping.")
        return
    
    files = sorted([f for f in glob.glob(os.path.join(net_dir, f"{file_pre}*.csv"))])
    df = pd.concat([pd.read_csv(f) for f in files])
    df = utils.agg_data(df)
    df.to_csv(os.path.join(out_dir, out_file), index=False)
    print(f"Saved {out_file}")


# def agg_data(
#     df, col='pet', roll_args={'window':60, 'min_periods':1}, 
#     names=['total60d','total_month']
# ):
#     df.time = pd.to_datetime(df.time, format='mixed') #format='%Y-%m-%d')
#     df.set_index('time', inplace=True)
#     df.sort_index(inplace=True)

#     if not names[0] in df.columns:
#         # Calculate rolling 60-day sum
#         df[names[0]] = df[col].rolling(**roll_args).sum()
#     if not names[1] in df.columns:
#         # Calculate monthly sum
#         df[names[1]] = df[col].resample('MS').sum()
#         df[names[1]] = df[names[1]].ffill()

#     return df


# def agg_file(
#     file, col='pet', roll_args={'window':60, 'min_periods':1}, 
#     names=['total60d','total_month']
# ):
#     df = pd.read_csv(file)
#     if all([n in df.columns for n in names]):
#         print(f"File {os.path.basename(file)} already aggregated. Skipping.")
#         # return df
#         return
#     try:
#         df = agg_data(df, col=col, roll_args=roll_args, names=names)
#     except Exception as e:
#         print(f"Error aggregating {os.path.basename(file)}: {e}")
#         return

#     df.to_csv(file)
#     print(f"Saved {os.path.basename(file)}")

    # return df


# MAIN
def main():
    # Get list of files
    pet_files = sorted([f for f in glob.glob(os.path.join(pet_dir, 'byyear', '*.csv'))])
    # pet_files.sort()

    # Loop through files
    # for file in pet_files:
    for yr in range(1980, 2024, 10):
        print(f"Processing {yr}s")
        files = sorted([f for f in pet_files if str(yr)[:-1] in f])
        # Load file
        # df = pd.read_csv(file)
        print(f"Reading {len(files)} files")
        df = pd.concat([pd.read_csv(f) for f in files])
        # Get year
        # yr = int(re.search(r'\d{4}', file).group())
        with mp.Pool(20) as pool:
            pool.imap(separate_df, df.iterrows())
        # Separate by network and station
        separate_df(
            df, parent_dir=os.path.join(pet_dir, 'tmp'), suff=f'_{yr}s'
        )
        print(f"Finished processing {yr}s")
    
    if os.path.exists(anc_file):
        ismn = pd.read_csv(anc_file)
    else:
        from soil.station import META
        ismn = META[['network','station']].drop_duplicates()

    with mp.Pool(20) as pool:
        pool.starmap(
            process_station, 
            [(row, pet_dir, os.path.join(pet_dir, 'tmp')) for i, row in ismn.iterrows()]
        )
    pool.close()
    pool.join()

    # NOTE: code below was for running alone to add columns to existing files.
    # This is now incorporated in process_station, so shouldn't be necessary.

    # pet_files = sorted([f for f in glob.glob(os.path.join(pet_dir, '*.csv'))])
    # roll_args = {'window': 60, 'min_periods': 60}
    # names = ['total60d', 'total_month']

    # with mp.Pool(50) as pool:
    #     pool.starmap(
    #         utils.agg_file, 
    #         [(f, 'pet', roll_args, names, True) for f in pet_files]
    #     )
    # pool.close()
    # pool.join()

#
if __name__ == "__main__":
    main()