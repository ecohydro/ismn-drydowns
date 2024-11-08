#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Tue 12 Mar 24 12:44:28'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           calc_map.py
Compatibility:  Python 3.12.0

Description:    Calculate annual rainfall totals and mean annual precipitation (MAP)
                from CHIRPS data.

Updated:        12 Mar 2024

URL:            https://

Requires:       list of libraries required

Dev ToDo:       None

AUTHOR:         Bryn Morgan
ORGANIZATION:   University of California, Santa Barbara
Contact:        bryn.morgan@geog.ucsb.edu
Copyright:      (c) Bryn Morgan 2024


"""

#-------------------------------------------------------------------------------
# IMPORTS
#-------------------------------------------------------------------------------

#Packages
import os
import glob
import datetime

import rainfall

#-------------------------------------------------------------------------------
# READ DATA
#-------------------------------------------------------------------------------
chirps_dir = '/home/chc-data-out/products/CHIRPS-2.0/global_daily/netcdf/'
res = 'p05'

# Use .25° for testing; will use .05° for final analysis
fp_chirps = os.path.join(chirps_dir, res)
fp_chirps_monthly = os.path.join(chirps_dir, res, 'by_month')
# '/home/chc-data-out/products/CHIRPS-2.0/global_daily/netcdf/p25/by_month/'

# Get all files
# files = glob.glob(fp_chirps + '*.nc')
# Only get files pre-2024
files_yr = [f for f in glob.glob(fp_chirps + '/*.nc') if int(f[-16:-12]) < 2024]
files_mth = [f for f in glob.glob(fp_chirps_monthly + '/*.nc') if int(f[-19:-15]) < 2024]
files_yr.sort()
files_mth.sort()


#-------------------------------------------------------------------------------
# RUN CODE 
#-------------------------------------------------------------------------------
def main():

    yrly_tot, ds_map = rainfall.calc_map(files_yr, mask=None)
    
    # SAVE TO NETCDF

    out_dir = '/home/waves/projects/ecohydro-regimes/outputs/chirps/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"Created directory: {out_dir}")
    else:
        print(f"Directory already exists: {out_dir}")

    date = datetime.datetime.now().strftime('%d%b').lower()

    # suff = f"_monthly_{res}_7mar"
    suff = f"_{res}_{date}"
    # alpha.to_netcdf(os.path.join(out_dir, f'alpha{id}.nc'))
    # freq.to_netcdf(os.path.join(out_dir, f'lambda{id}.nc'))

    rainfall.save_netcdf(yrly_tot, out_dir, f'total_annual{suff}')
    rainfall.save_netcdf(ds_map, out_dir, f'map{suff}')


if __name__ == '__main__':
    main()

