#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Sun 21 Apr 24 13:06:57'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           extract_rainfall.py
Compatibility:  Python 3.12.0

Description:    Code for extracting data from rainfall stats.

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

out_dir = os.path.join(proj_dir, 'outputs')

chirps_dir = os.path.join(out_dir, 'chirps')

# grid_dir = os.path.join(data_dir, 'grids')

# chirps_out_dir = os.path.join(proj_dir, 'outputs', 'chirps')

anc_file = os.path.join(data_dir, 'ISMN', 'ismn_ancillary.csv')

fdate_dict = {
    'p05': '27mar',
    'p25': '20mar',
    'map' : '12mar'
}

# Use .25° for testing; will use .05° for final analysis

# FUNCTIONS
def extract_xr_values(
    df, ds_in, 
    subset_args = {'dim':'location', 'lat_col':'latitude', 'lon_col':'longitude'},
    return_ds = False
):
    # dim='location', lat_col='latitude', lon_col='longitude', 
    # ):
    # Subset the dataset
    ds = utils.subset_xr(df, ds_in, **subset_args)
    # Extract values to dataframe
    print(f"Extracting values from {len(df)} stations")
    vals = ds.to_dataframe().reset_index()
    # Add network + station cols and remove location
    vals = utils.add_network_station(vals)
    print(f"Extracted {len(vals)/len(df)} values each from {len(df)} stations")
    if return_ds:
        return vals, ds
    else:
        return vals

# def add_network_station(df):
#     cols = list(df.columns)
#     df['network'], df['station'] = zip(*df['location'])
#     return df[['network','station'] + cols[1:]]



# MAIN
def main():
    
    res = 'p05'

    files = [
        os.path.join(
            chirps_dir, utils.get_filename(pre, res, fdate_dict.get(res))
        ) for pre in ['alpha', 'lambda', 'cv']
    ]

    ds = xr.open_mfdataset(files)

    if os.path.exists(anc_file):
        ismn = pd.read_csv(anc_file)
    else:
        from soil.station import META
        ismn = META[['network','station','latitude','longitude']].drop_duplicates()

    ismn['location'] = ismn.apply(lambda x: tuple(x[['network','station']]), axis=1)

    vals = extract_xr_values(ismn, ds)

    vals.to_csv(
        os.path.join(out_dir, f'ismn_rain_stats_{res}.csv'), index=False
    )

    ds.close()
    
if __name__ == "__main__":
    main()