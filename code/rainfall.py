#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Tue 12 Mar 24 12:47:18'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           rainfall.py
Compatibility:  Python 3.712.0

Description:    Utility functions for calculating rainfall statistics.

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
import os
import numpy as np
import xarray as xr

import calendar

from scipy import stats
# from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor

#-------------------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------------------


def calc_lambda(ds, thresh, mask=None):
    """
    Calculate the frequency of rainfall for a given threshold.

    Parameters
    ----------
    ds : xarray.Dataset
        CHIRPS dataset.
    thresh : float
        Rainfall threshold for defining a rainy day. The default is 0.1 mm.
    mask : xarray.DataArray
        Mask for removing NaN values. The default is None.
    
    Returns
    -------
    freq : xarray.DataArray
        Rainfall frequency, lambda [days-1].
    n_rainy : xarray.DataArray
        Number of rainy days.
    """
    # Get nan mask
    if mask is None:
        mask = ~np.isnan(ds['precip'].isel(time=0)).drop_vars('time')
    # Calculate the number of rainy days
    n_rainy = (ds['precip'] > thresh).sum(dim='time').where(mask)
    # Number of days
    n_days = ds.sizes['time']
    # Calculate frequency of rainfall.
    freq = n_rainy / n_days
    return freq, n_rainy

def calc_lag_lambda(ds, thresh, mask=None, lag=1, groupby='time.month'):
    """
    Calculate the frequency of rainfall for a given threshold.

    Parameters
    ----------
    ds : xarray.Dataset
        CHIRPS dataset.
    thresh : float
        Rainfall threshold for defining a rainy day. The default is 0.1 mm.
    mask : xarray.DataArray
        Mask for removing NaN values. The default is None.
    lag : int
        Number of days to lag the rainfall. The default is 1.
    
    Returns
    -------
    freq : xarray.DataArray
        Rainfall frequency, lambda [days-1].
    n_rainy : xarray.DataArray
        Number of rainy days.
    """
    # Get nan mask
    if mask is None:
        mask = ~np.isnan(ds['precip'].isel(time=0)).drop_vars('time')
    # Shift by lag
    yesterday = ds['precip'].shift(time=lag)
    # Masks
    if groupby:
        n_rainy_yest = (yesterday > thresh).groupby(groupby).sum(dim='time').where(mask)
        n_rainy_both = ((ds['precip'] > thresh) & (yesterday > thresh)).groupby(groupby).sum(dim='time').where(mask)
    else:
        # Calculate the number of days with rain the day before
        n_rainy_yest = (yesterday > thresh).sum(dim='time').where(mask)
        # Calculate the number of days with rain today and yesterday
        n_rainy_both = ((ds['precip'] > thresh) & (yesterday > thresh)).sum(dim='time').where(mask)
    # Probability of rain today if it rained yesterday
    freq_both = n_rainy_both / n_rainy_yest
    
    return freq_both



def calc_mean_alpha(ds, n_rainy, mask=None):
    """
    Calculate the mean rainfall intensity.

    Parameters
    ----------
    ds : xarray.Dataset
        CHIRPS dataset.
    n_rainy : xarray.DataArray
        Number of rainy days.
    mask : xarray.DataArray
        Mask for removing NaN values. The default is None.

    Returns
    -------
    alpha : xarray.DataArray
        Rainfall intensity, alpha [mm].
    """
    # Get nan mask
    if mask is None:
        mask = ~np.isnan(ds['precip'].isel(time=0)).drop_vars('time')
    # Calculate monthly alpha
    alpha = xr.where(
        n_rainy > 0,        # If there are rainy days
        ds['precip'].sum(dim='time').where(mask) / n_rainy,     # Calculate mean intensity
        n_rainy            # else, return the number of rainy days (0)
    )
    return alpha

def fit_exponential(arr, var='precip'):#, return_stats=True):
    """
    Fit an exponential distribution to the data for a given latitude and longitude.

    Parameters
    ----------
    ds : xarray.Dataset
        CHIRPS dataset.
    lat : float
        Latitude.
    lon : float
        Longitude.
    var : str
        Variable to fit the distribution to. The default is 'precip'.

    Returns
    -------
    scale : np.array
        Scale parameter (mean) for the exponential distribution, Kolmorogov-Smirnov
        test statistic, and p-value.
    """
    # Extract the data for the given latitude and longitude along time dimension
    # arr = ds[var].sel(latitude=lat, longitude=lon).values.flatten()
    # Fit the exponential distribution to the data
    try:
        loc, scale = stats.expon.fit(arr[arr > 0.], floc=0)
        ks, pval = stats.kstest(arr[arr > 0.], 'expon', args=(0, scale))
    except ValueError:
        scale = np.nan
        ks = np.nan
        pval = np.nan
    # if return_stats:
        # if isinstance(arr, xr.DataArray):
    # if hasattr(arr, 'coords'):
    #     coords = {coord : arr[coord] for coord in arr.coords if coord not in arr.dims}
        # scale = xr.DataArray(scale, coords=coords, dims=coords.keys())#dims=arr.dims[1:])
        # scale = xr.DataArray(scale, coords=coords, dims=arr.dims[1:])
    # print(f"Alpha for lat: {lat}, lon: {lon} is {scale}")
    return np.array([scale, ks, pval])


def calc_alpha(ds, var='precip', parallel=False):
    """
    Calculate the mean rainfall intensity across the CHIRPS record.

    Parameters
    ----------
    ds : xarray.Dataset
        CHIRPS dataset.
    pixels : list
        List of tuples containing latitude and longitude.
    var : str
        Variable to fit the distribution to. The default is 'precip'.
    parallel : bool
        Whether to use parallel processing. The default is False.

    Returns
    -------
    scales_xr : xarray.Dataset
        Contains 3 variables: 'alpha,' the mean rainfall intensity, 'ks_stat,' the
        Kolmogorov-Smirnov test statistic, and 'p_value,' the p-value from the
        Kolmogorov-Smirnov test.
    """
    # Get list of pixels
    # pixels = [(lat, lon) for lat in ds['latitude'].values for lon in ds['longitude'].values]
    # if parallel:
    #     # Create a multiprocessing Pool
    #     # with Pool(processes=None) as pool:  # Use None to automatically detect number of CPU cores
    #     #     # Use map to apply the function to each pixel in parallel
    #     #     scales = pool.starmap(fit_exponential, [(ds, lat, lon, var) for lat, lon in pixels])
    #     args = [
    #         {
    #             'arr': ds[var].sel(latitude=lat, longitude=lon), 
    #             'lat': lat, 'lon': lon, 'var': 'precip'
    #         } for lat in ds['latitude'].values for lon in ds['longitude'].values
    #     ]
    #     with ThreadPoolExecutor() as executor:
    #         scales = list(executor.map(lambda x: fit_exponential(**x), args))
    # # else:
    # #     # Use a for loop to apply the function to each pixel in series
    # #     scales = [
    # #         fit_exponential(ds[var].sel(latitude=lat, longitude=lon), lat, lon, var) for lat, lon in pixels
    # #     ]
    # # # Reshape the scale parameters array to match the shape of the latitude and longitude arrays
    # # scales = np.array(scales).reshape(len(ds['latitude']), len(ds['longitude']))
    # # # Convert the scale parameters to a DataArray
    # # scales_xr = xr.DataArray(
    # #     scales, coords=[ds['latitude'], ds['longitude']], dims=['latitude', 'longitude']
    # # )

    scales_xr = xr.apply_ufunc(
        fit_exponential, ds[var], input_core_dims=[['time']], 
        output_core_dims=[['params']], output_sizes={'params': 3},
        dask='allowed', vectorize=True,
    )

    scales_xr = scales_xr.to_dataset(dim='params').rename_vars(
        {0: 'intensity', 1: 'ks_stat', 2: 'p_value'}
    )

    return scales_xr

#-------------------------------------------------------------------------------
# MONTHLY STUFF
#-------------------------------------------------------------------------------

def get_files_by_month(files, month=1):
    return [f for f in files if int(f[-14:-12]) == month]

def calc_month(file):
    """
    Calculate total rainfall and number of rainy days for a single month.

    Parameters
    ----------
    file : str
        Path to the CHIRPS file.

    Returns
    -------
    mth_rain : xarray.DataArray
        Monthly rainfall sums.
    n_rainy : xarray.DataArray
        Number of rainy days.
    """
    # Read in the file
    ds = xr.open_dataset(file)
    # ds = xr.open_dataset(mth_files[0])
    # Calculate the monthly sums
    # monthly_sums = ds['precip'].groupby('time.month').sum(dim='time')
    mth_rain = ds.sum(dim='time') #.to_dataarray(dim='precip')
    # n_rainy = (ds['precip'] > 0).groupby('time.month').sum(dim='time')
    n_rainy = (ds['precip'] > 0).sum(dim='time')

    ds.close()

    return mth_rain, n_rainy



def calc_monthly(files, month=1, thresh=0.1, mask=None):
    """
    Calculate rainfall frequency and intensity for a given month across the CHIRPS record.

    Parameters
    ----------
    files : list
        List of CHIRPS files.
    month : int
        Month to calculate.
    thresh : float
        Rainfall threshold for defining a rainy day.
    mask : xarray.DataArray
        Mask for removing NaN values.

    Returns 
    -------
    alpha : xarray.DataArray
        Rainfall intensity, alpha [mm].
    freq : xarray.DataArray
        Rainfall frequency, lambda [days-1].
    """
    # Get files for the month
    mth_files = get_files_by_month(files, month)
    # Read in the file
    ds = xr.open_mfdataset(mth_files)
    # Get nan mask
    if mask is None:
        mask = ~np.isnan(ds['precip'].isel(time=0)).drop_vars('time')
    # Calculate number of rainy days
    n_rainy = (ds['precip'] > thresh).sum(dim='time').where(mask)
    # Number of days
    n_days = ds.sizes['time']
    # Calculate frequency
    freq = n_rainy / n_days
    # Calculate the monthly sums
    # tot_rain = ds['precip'].sum(dim='time').where(mask)
    # Calculate monthly alpha
    # alpha = ds['precip'].sum(dim='time').where(mask) / n_rainy

    alpha = xr.where(n_rainy > 0, ds['precip'].sum(dim='time').where(mask) / n_rainy, n_rainy)

    ds.close()

    # return tot_rain, n_rainy
    return alpha, freq


def calc_all_months(files_mth, thresh=0.1, mask=None):
    """
    Calculate rainfall statistics for all months across the CHIRPS record.

    Parameters
    ----------
    thresh : float
        Rainfall threshold for defining a rainy day.
    mask : xarray.DataArray
        Mask for removing NaN values.

    Returns
    -------
    alpha_all : xarray.Dataset
        Rainfall intensity, alpha [mm]. The dataset contains 3 dimensions: 
        latitude, longitude, month.

    freq_all : xarray.Dataset
        Rainfall frequency, lambda [days-1]. The dataset contains 3 dimensions:
        latitude, longitude, month.
    """
    for month in range(1, 13):
        alpha, freq = calc_monthly(files_mth, month, thresh, mask)
        if month == 1:
            alpha_all = alpha
            freq_all = freq
        else:
            alpha_all = xr.concat([alpha_all, alpha], dim='month')
            freq_all = xr.concat([freq_all, freq], dim='month')
    alpha_all = alpha_all.assign_coords(month=np.arange(1,13)).to_dataset().rename_vars({'precip': 'intensity'})

    freq_all = freq_all.assign_coords(month=np.arange(1,13)).to_dataset().rename_vars({'precip': 'frequency'})

    # return tot_rain, n_rainy
    return alpha_all, freq_all


def calc_lag_lambda_monthly(files_all, thresh=0.1, mask=None):
    ds_yr = xr.open_mfdataset(files_all)
    freq_both = calc_lag_lambda(ds_yr, thresh=0.1, mask=None, lag=1, groupby='time.month')
    ds_yr.close()
    freq_both = freq_both.to_dataset().rename_vars({'precip': 'lag1_probability'})
    return freq_both


def process_month(files_mth, month, thresh, mask=None, parallel=False):
    # 1. Read in data for the month
    # Get files for the month
    mth_files = get_files_by_month(files_mth, month)
    print(f"Found {len(mth_files)} files for month {month}")
    # Read in the file
    ds = xr.open_mfdataset(mth_files)
    print(f"Reading files for month {month}")
    ds.load()
    # 2. Calculate the frequency of rainfall + number of rainy days
    print(f"Calculating rainfall frequency for month {month}")
    freq, n_rainy = calc_lambda(ds, thresh, mask)
    # 3. Calculate the average intensity of rainfall
    print(f"Calculating mean intensity for month {month}")
    alpha_mn = calc_mean_alpha(ds, n_rainy, mask)
    # 4. Calculate the intensity of rainfall for each pixel by fitting an exponential distribution
    # Get list of pixels
    pixels = [(lat, lon) for lat in ds['latitude'].values for lon in ds['longitude'].values]
    # 5. Calculate the scale parameter for each pixel
    print(f"Calculating rainfall intensity for month {month}")
    alpha = calc_alpha(ds, parallel=parallel)

    print(f"Finished processing month {month}")# \n freq = {freq} days^-1 \n alpha_mn = {alpha_mn} mm \n alpha = {alpha} mm")
    # 6. Close the dataset
    ds.close()

    return freq, alpha_mn, alpha


def process_months(files_mth, thresh=0.1, mask=None, parallel=False, calc_lagged=True):
    # Iterate through months
    for month in range(1, 13):
        # Calculate frequency, mean intensity, and intensity (exponential fit) for each month
        freq, alpha_mn, alpha = process_month(files_mth, month, thresh, mask, parallel)
        # If it's the first month, create the datasets
        if month == 1:
            alpha_all = alpha
            freq_all = freq
            alpha_mn_all = alpha_mn
        # Otherwise, concatenate
        else:
            alpha_all = xr.concat([alpha_all, alpha], dim='month')
            freq_all = xr.concat([freq_all, freq], dim='month')
            alpha_mn_all = xr.concat([alpha_mn_all, alpha_mn], dim='month')
    # Assign coordinates and rename variables
    alpha_all = alpha_all.assign_coords(month=np.arange(1,13))
        # month=np.arange(1,13)).to_dataset().rename_vars({'precip': 'intensity'})
    freq_all = freq_all.assign_coords(
        month=np.arange(1,13)).to_dataset().rename_vars({'precip': 'frequency'})
    alpha_mn_all = alpha_mn_all.assign_coords(
        month=np.arange(1,13)).to_dataset().rename_vars({'precip': 'intensity'})
    if calc_lagged:
        print(f"Calculating lagged rainfall probability by month")
        freq_both = calc_lag_lambda_monthly(files_mth, thresh=0.1, mask=None)
        freq_all = xr.merge([freq_all,freq_both])
        # fits_all = xr.merge([alpha_all.ks_stat, alpha_all.p_value, freq_both])

    return freq_all, alpha_mn_all, alpha_all


#-------------------------------------------------------------------------------
# COEFFICIENT OF VARIATION
#-------------------------------------------------------------------------------
def get_n_days(month):
    if month != 2:
        return calendar.mdays[month]
    else:
        return (33*28 + 10*29)/43



def calc_mu(alpha, freq):
    """
    Calculate the mean event size, mu.

    Parameters
    ----------
    alpha : xarray.DataArray
        Rainfall intensity, alpha [mm].
    freq : xarray.DataArray
        Rainfall frequency, lambda [days-1].

    Returns
    -------
    mu : xarray.Dataset
        Total rainfall [mm].
    """
    n_days = xr.DataArray(
        [get_n_days(month) for month in range(1,13)], 
        dims=('month'), coords={'month': range(1, 13)}
    )

    mu = alpha * freq * n_days
    
    return mu

def calc_cv(alpha, freq, mu=None):
    """
    Calculate the coefficient of variation (CV) for a given month.

    Parameters
    ----------
    alpha : xarray.DataArray
        Rainfall intensity, alpha [mm].
    freq : xarray.DataArray
        Rainfall frequency, lambda [days-1].

    Returns
    -------
    cv : xarray.DataArray
        Coefficient of variation.
    """
    if mu is None:
        mu = calc_mu(alpha, freq)

    cv = np.sqrt(2) * (alpha/mu)**(1/2) #(mu / alpha)**(-1/2)

    cv = cv.to_dataset(name='coeff_var')

    return cv

#-------------------------------------------------------------------------------
# YEARLY RAINFALL TOTALS + MEAN ANNUAL PRECIPITATION
#-------------------------------------------------------------------------------

# def calc_yearly_rain(file, month=1, mask=None):
#     """
#     Calculate yearly rainfall statistics across the CHIRPS record.

#     Parameters
#     ----------
#     files : list
#         List of CHIRPS files.
#     thresh : float
#         Rainfall threshold for defining a rainy day.
#     mask : xarray.DataArray
#         Mask for removing NaN values.

#     Returns
#     -------
#     alpha_all : xarray.Dataset
#         Rainfall intensity, alpha [mm]. The dataset contains 3 dimensions: 
#         latitude, longitude, year.

#     freq_all : xarray.Dataset
#         Rainfall frequency, lambda [days-1]. The dataset contains 3 dimensions:
#         latitude, longitude, year.
#     """
#     # Open file
#     ds = xr.open_dataset(file)
#     # Get nan mask
#     if mask is None:
#         mask = ~np.isnan(ds['precip'].isel(time=0)).drop_vars('time')
#     # Calculate total rainfall for the year

#     # return tot_rain, n_rainy
#     return alpha_all, freq_all


def calc_map(files_yr, mask=None):
    """
    Calculate the mean annual precipitation across the CHIRPS record.

    Parameters
    ----------
    files_yr : list
        List of CHIRPS files.
    mask : xarray.DataArray
        Mask for removing NaN values.

    Returns
    -------
    ds_map : xarray.Dataset
        Mean annual precipitation [mm/yr].
    """
    # Calculate yearly rainfall totals
    ds_all = xr.open_mfdataset(files_yr)
    # Get mask
    if mask is None:
        mask = ~np.isnan(ds_all['precip'].isel(time=0)).drop_vars('time')
    # Calculate sums
    yrly_tot = ds_all.groupby('time.year').sum(dim='time').where(mask)

    yrly_tot['precip'].attrs['units'] = 'mm/yr' #r"mm yr$^{-1}$"

    # Calculate MAP
    ds_map = yrly_tot.mean(dim='year').where(mask)

    ds_all.close()

    return yrly_tot, ds_map



#-------------------------------------------------------------------------------
# SAVE
#-------------------------------------------------------------------------------
def save_netcdf(ds, out_dir, id):
    filename = f'{id}.nc'
    print(f"Saving {filename} to {out_dir}")
    ds.to_netcdf(os.path.join(out_dir, f'{id}.nc'))
    print(f"Saved {filename} to {out_dir}")

