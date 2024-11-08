#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Thu 07 Mar 24 17:46:15'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           calc_rainfall_monthly.py
Compatibility:  Python 3.12.0

Description:    Calculate rainfall statistics (α, λ) for each month across the 
                CHIRPS record.

Updated:        07 Mar 2024

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
# import calendar

import rainfall

# READ DATA
chirps_dir = '/home/chc-data-out/products/CHIRPS-2.0/global_daily/netcdf/'
res = 'p05'
thresh = 0.1

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

def main():

    freq, alpha_mn, alpha = rainfall.process_months(
        files_mth, thresh, mask=None, parallel=False, calc_lagged=False,
    )

    cv = rainfall.calc_cv(alpha.intensity, freq.frequency)

    date = datetime.datetime.now().strftime('%d%b').lower()
    # date = datetime.datetime.now().strftime('%Y_%m_%d')

    # out_dir = f'/home/waves/projects/ecohydro-regimes/outputs/chirps/{date}'
    out_dir = f'/home/waves/projects/ecohydro-regimes/outputs/chirps/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        # os.makedirs(os.path.join(out_dir, 'stats'))
        print(f"Created directory: {out_dir}")
    else:
        print(f"Directory already exists: {out_dir}")

    # names = ['lambda', 'mean_alpha', 'alpha', 'cv']
    # names = {
    #     'lambda' : freq.frequency,
    #     'mean_alpha' : alpha_mn.intensity,
    #     'alpha' : alpha.intensity,
    #     'cv'  : cv
    # }
    names = {
        'lambda' : freq.frequency,
        'mean_alpha' : alpha_mn.intensity,
        'alpha' : alpha.intensity,
        'cv'  : cv,
        'alpha_stats' : alpha[['ks_stat','p_value']],
        # 'lambda_stats' : freq.lag1_probability
    }
    # suff = f'_monthly_{res}'
    suff = f'_monthly_{res}_{date}'

    # for i, arr in enumerate([freq, alpha_mn, alpha, cv]):
    for name,ds in names.items():
        rainfall.save_netcdf(ds, out_dir, f'{name}{suff}')


if __name__ == '__main__':
    main()
