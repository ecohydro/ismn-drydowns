#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Sun 21 Apr 24 17:39:35'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           calc_mean_pet.py
Compatibility:  Python 3.7.0

Description:    Calculate mean annual potential evapotranspiration (PET) from
                daily PET (dPET) data.

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
import multiprocessing as mp


# import drydowns
from drydowns.config import Config

# import soil.station as station
# from soil.station import Station, META, ismn_data

# import utils


#% DIRECTORIES

config = Config('soil/config.ini')
# cfg = config[config.get('RUN','type')]
proj_dir = config.config.get('ISMN', 'project_dir')
data_dir = os.path.join(proj_dir, 'data', 'ISMN')
# out_dir = os.path.join(proj_dir, 'outputs')
# out_dir = os.path.join(data_dir, 'dPET', 'byyear')
out_dir = os.path.join(proj_dir, 'outputs', 'pet')

chirps_dir = os.path.join(out_dir, 'chirps')
grid_dir = os.path.join(data_dir, 'grids')

pet_dir = os.path.join(config.config.get('SMAP','project_dir'), 'data', 'PET')

# pet_files = [f for f in glob.glob(os.path.join(pet_dir, '*_daily_pet.nc'))]
# pet_files.sort()

# FUNCTIONS
def calc_mean(file):

    print(f"Opening {file}.")
    ds = xr.open_dataset(file)
    year = ds['time.year'][0].item()
    print(f"Masking nans for {year}.")
    mask = ~np.isnan(ds.where(ds.pet < 100.)).pet.isel(time=0).drop_vars('time')
    print(f"Calculating total PET for {year}.")
    yr = ds.where(ds.pet < 100.).where(ds.pet > 0.).sum(dim='time').where(mask)

    yr = yr.assign_coords(year=year)
    yr = yr.expand_dims('year')

    ds.close()

    print(f"Finished calculating total PET for {year}.")

    return yr


# MAIN
def main():

    pet_files = sorted([f for f in glob.glob(os.path.join(pet_dir, '*_daily_pet.nc'))])

    with mp.Pool(5) as pool:

        yr_list = list(pool.map(calc_mean, pet_files))
    pool.close()
    pool.join()
    
    yrs = xr.concat(yr_list, dim='year')
    print(f"Finished calculating total PET for all years. Calculating mean annual PET.")
    mask = ~np.isnan(yr_list[0].pet.isel(year=0)).drop_vars('year')
    mean_pet = yrs.mean(dim='year').where(mask)
    print(f"Saving mean annual PET to {out_dir}.")
    # Save
    mean_pet.to_netcdf(os.path.join(out_dir, 'mean_annual_pet.nc'))



if __name__ == '__main__':
    main()
