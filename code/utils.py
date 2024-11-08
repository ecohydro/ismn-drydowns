#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Wed 13 Mar 24 15:37:15'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           utils.py
Compatibility:  Python 3.12.0

Description:    Utility functions for preprocessing.

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
import pandas as pd
import xarray as xr
import multiprocessing as mp
import ast
#-------------------------------------------------------------------------------
# VARIABLES
#-------------------------------------------------------------------------------


var_dict = {
    'alpha' : {
        # 'ds' : alpha,
        'var' : 'intensity',
        'label' : r"$\alpha$ (mm)",
        'kwargs' : {
            'vmin' : 0,
            'vmax' : 15,
            'cmap' : 'YlGnBu',
        }
    },
    'lambda' : {
        # 'ds' : freq,
        'var' : 'frequency',
        'label' : r"$\lambda$ (days$^{-1}$)",
        'kwargs' : {
            'vmin' : 0,
            'vmax' : 1,
            'cmap' : 'viridis',
        },
    },
    'cv' : {
        # 'ds' : cv,
        'var' : 'coeff_var',
        'label' : r"Coefficient of variation (CV)",
        'kwargs' : {
            'vmin' : 0,
            'vmax' : 3,
            'cmap' : 'YlGnBu',
        },
    },
    'mu' : {
        # 'ds' : mu,
        'var' : 'precip',
        'label' : r"Mean monthly precipitation (mm)",
        'kwargs' : {
            'vmin' : 0,
            'vmax' : 300,
            'cmap' : 'YlGnBu', 
        },
    },
    'map' : {
        # 'ds' : map25,
        'var' : 'precip',
        # 'label' : r"MAP (mm)", # (mm yr$^{-1}$)",
        'label' : r"Mean annual precipitation (mm)", # (mm yr$^{-1}$)",
        'kwargs' : {
            'vmin' : 0,
            'vmax' : 4000,
            'cmap' : 'YlGnBu',
        },
    },
    'tap' : {
        # 'ds' : tap25,
        'var' : 'precip',
        'label' : r"Total annual precipitation (mm)",
        'kwargs' : {
            'vmin' : 0,
            'vmax' : 4000,
            'cmap' : 'YlGnBu',
        },
    },
    'theta' : {
        'symbol' : r"$\theta$",
        'label' : r"$\theta$ (m$^3$ m$^{-3}$)",
        'units' : r"m$^3$ m$^{-3}$",
    },
    'dtheta' : {
        'symbol' : r"$\frac{d\theta}{dt}$",
        'label' : r"$\frac{d\theta}{dt}$ (m$^3$ m$^{-3}$ day$^{-1}$)",
        'units' : r"m$^3$ m$^{-3}$ day$^{-1}$",
    },
    'dtheta_mm' : {
        'symbol' : r"$\frac{d\theta}{dt}$",
        'label' : r"$\frac{d\theta}{dt}$ (mm day$^{-1}$)",
        'units' : r"mm day$^{-1}$",
    },
    't' : {
        'symbol' : r"$t$",
        'label' : r"$t$ (days)",
        'units' : r"days",
    },
    'et' : {
        'symbol' : r"$ET$",
        'label' : r"$ET$ (mm day$^{-1}$)",
        'units' : r"mm day$^{-1}$",
    },
    'LAI' : {
        'symbol' : r"LAI",
        'label' : r"LAI (m$^2$ m$^{-2}$)",
        'units' : r"m$^2$ m$^{-2}$",
    },
    # "theta": {
    #     "col": "sm",
    #     "symbol": r"$\theta$",
    #     # "label": r"Soil moisture",
    #     "label": r"Soil moisture, $\theta$",
    #     "unit": r"(m$^3$ m$^{-3}$)",
    #     "lim": [0, 0.50],
    # },
    # "dtheta": {
    #     "col": "",
    #     "symbol": r"$-\frac{d\theta}{dt}$",
    #     "label": r"Change in soil moisture",
    #     # "label": r"Change in soil moisture, $-\frac{d\theta}{dt}$",
    #     "unit": r"(m$^3$ m$^{-3}$ day$^{-1}$)",
    #     "lim": [-0.10, 0],
    # },
    "theta_w": {
        "col": "min_sm",
        "symbol": r"$\theta_w$",
        "label": r"Estimated $\theta_w$",
        "unit": r"(m$^3$ m$^{-3}$)",
        "lim": [0.1, 0.45],
    },
    "q": {
        "col": "q_q",
        "symbol": r"$q$",
        # "label": r"Nonlinear parameter $q$",
        "label": r"$q$",
        # "unit": "[-]",
        "unit": "",
        "lim": [0.5, 4.0],
    },
    "ET_max": {
        "col": "q_ETmax",
        "symbol": r"$ET_{\mathrm{max}}$",
        "label": r"Estimated $ET_{max}$",
        "unit": r"(mm day$^{-1}$)",
        "lim": [0, 17.5],
    },
    "theta_star": {
        "col": "max_sm",
        "symbol": r"$\theta_*$",
        "label": r"Estimated $\theta_*$",
        "unit": r"(m$^3$ m$^{-3}$)",
        "lim": [0.1, 0.45],
    },
    "sand_bins": {
        "col": "sand_bins",
        "symbol": r"",
        "label": r"Sand fraction",
        "unit": "[-]",
        "lim": [0.0, 1.0],
    },
    "ai_bins": {
        "col": "ai_bins",
        "symbol": r"AI",
        "label": r"Aridity Index",
        "unit": "(MAP/MAE)",
        "lim": [0.0, 2.0],
    },
    "veg_class": {
        "col": "name",
        "symbol": r"",
        "label": r"IGBP Landcover Class",
        "unit": "",
        "lim": [0, 1],
    },
    "ai": {
        "col": "AI",
        "symbol": r"AI",
        "label": r"Aridity Index",
        "unit": "(MAP/MAE)",
        "lim": [0.0, 1.1],
    },
    "diff_R2": {
        "col": "diff_R2",
        "symbol": r"$R^2$",
        "label": r"$R^2$ (Nonlinear - linear)",
        "unit": "[-]",
        "lim": [-0.02, 0.02],
    },
}


#-------------------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------------------
def get_filename(var, res, fdate):
    file_dict = {
        'alpha' : f"alpha_monthly_{res}_{fdate}.nc",
        'alpha_mn' : f"mean_alpha_monthly_{res}_{fdate}.nc",
        'lambda' : f"lambda_monthly_{res}_{fdate}.nc",
        'cv' : f"cv_monthly_{res}_{fdate}.nc",
        'map' : f"map_{res}_{fdate}.nc",
        'tap' : f"total_annual_{res}_{fdate}.nc",
    }
    return file_dict.get(var)

def open_file(file, path, resample=False, grid=None):
    print(f"Opening {file}")
    ds = xr.open_dataset(
        os.path.join(path, file)
    )
    if resample:
        print(f"Resampling {file}")
        ds.rio.set_crs("EPSG:4326")
        ds = resample_res(ds, grid)
        ds.rio.set_crs("EPSG:4326")
    
    return ds


def resample_res(in_array, proj_array):
    """
    Resample the resolution of a DataArray.

    Parameters
    ----------
    in_array : xarray.DataArray
        Array to be resampled
    proj_arr : xarray.DataArray
        DataArray to which to project the input array.

    Returns
    -------
    arr_resamp : xarray.DataArray
        Resampled array.

    TODO: Consider using rioxarray methods instead:
        arr_resamp = in_array.rio.reproject_match(proj_array)
        arr_resamp = arr_resamp.assign_coords({'x': in_array.x, 'y': in_array.y})
    This updates the resolution (etc.) and maintains handling of raster stuff.
    The downside is that it doesn't handle NoData well, but not an issue if 
    subsequently performing operations with an array with proper NoData.

    """
    in_array_coords = get_coord_names(in_array)
    proj_array_coords = get_coord_names(proj_array)

    arr_resamp = in_array.interp(
        **{
            in_array_coords[0]: proj_array[proj_array_coords[0]], 
            in_array_coords[1]: proj_array[proj_array_coords[1]]
        }
    )


    # arr_resamp = in_array.interp(
    #     longitude=proj_array.longitude, latitude=proj_array.latitude
    # )

    return arr_resamp

def get_val(df, var, stat, by=['EASE_row_index', 'EASE_column_index']):

    df_stat = df.groupby(
        by=by
    )[var_dict[var].get('col',var)].agg(stat) #.reindex(new_index, fill_value=np.nan)
    df_stat = df_stat.reset_index()

    return df_stat

def add_coord_info(df, var, coord_info,):
    new_index = pd.MultiIndex.from_tuples(
        zip(coord_info["EASE_row_index"], coord_info["EASE_column_index"]),
        names=["EASE_row_index", "EASE_column_index"],
    )
    df_pad = df.reindex(new_index, fill_value=np.nan)

    merged = (
        df_pad.reset_index().merge(
            coord_info[["EASE_row_index", "EASE_column_index", "latitude", "longitude"]], 
            on=["EASE_row_index", "EASE_column_index"],
            how="left",
        ).set_index(["EASE_row_index", "EASE_column_index"])
    )
    merged.rename(columns={var_dict[var].get('col',var) : var}, inplace=True)
    # merged = merged[merged.latitude > -60]  # Exclude antarctica in the map (no data)
    return merged


def create_pivot(df, var, stat, coord_info):
    # Get single value for variable
    df_stat = get_val(df, var, stat)
    # Add coordinate information
    merged = add_coord_info(df_stat, var, coord_info)
    # Create pivot array
    piv = merged.pivot(
        index="latitude", columns="longitude", values=var
    )
    # pivot_array[pivot_array.index > -60]  # Exclude antarctica in the map (no data)
    return piv


# def create_xarray(df, var, stat, coord_info):
#     # Get single value for variable
#     df_stat = get_val(df, var, stat)
#     # Add coordinate information
#     merged = add_coord_info(df_stat, var, coord_info)#.rename(columns={var: var_dict[var].get('col',var)})
#     # Create xarray
#     df_xr = merged.set_index(["latitude", "longitude"])[var].to_xarray()
#     return df_xr

def create_xarray(df : pd.DataFrame, grid : xr.DataArray, dims : list = ['latitude', 'longitude']):
    """
    Create xarray from dataframe. First pads the dataframe with all pixels in the grid.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe to convert to xarray.
    grid : xarray.DataArray
        Grid to use for the xarray lat/long coordinates.

    Returns
    -------
    ds : xarray.Dataset
        Dataset with lat/long coordinates and variable(s) from the dataframe.
    """
    # 1. Set index to lat/long
    # if 'latitude' not in df.index.names:
    #     df = df.set_index(['latitude', 'longitude'])
    # # 1. Pad with all pixels
    # index = grid.to_dataframe(name='precip').index
    # df_pad = df.reindex(index, fill_value=np.nan)

    # 1. Get lat/long from grid
    coords = grid.to_dataframe(name='precip').reset_index()
    del coords['precip']
    # 2. Merge with dataframe
    df_pad = df.merge(coords, on=["latitude", "longitude"], how="right")
    # 3. Set index
    df_pad.set_index(dims, inplace=True)
    # 2. Convert to xarray
    ds = df_pad.to_xarray()
    return ds


def create_grid(path : str, lat_file : str, lon_file : str, return_dfs=False):
    """
    Create xarray grid of latitudes and longitudes.

    Parameters
    ----------
    path : str
        Directory containing latitude and longitude files.
    lat_file : str
        File with latitude information.
    lon_file : str
        File with longitude information.

    Returns
    -------
    grid : xarray.DataArray
        Empty array of latitudes and longitudes.
    """
    # Read in lat/long files
    lats = pd.read_csv(os.path.join(path, lat_file))
    lons = pd.read_csv(os.path.join(path, lon_file))
    # Create xarray grid
    grid = xr.DataArray(
        coords={
            'latitude': sorted(lats.latitude.values), 
            'longitude': sorted(lons.longitude.values)
        },
        dims=['latitude', 'longitude']
    )
    if return_dfs:
        return grid, lats, lons
    else:
        return grid

def get_coord_names(ds):
    [lat, lon] = sorted([coord for coord in ds.coords if 'lat' in coord or 'lon' in coord])
    return lat, lon


def add_network_station(df):
    cols = list(df.columns)
    df['network'], df['station'] = zip(*df['location'])
    return df[['network','station'] + cols[1:]]

def subset_xr(df, ds, dim='location', lat_col='latitude', lon_col='longitude'):
    print(f"Subsetting dataset with {len(df)} stations")
    lats = xr.DataArray(df[lat_col].values, dims=dim, coords={dim: df[dim]})
    lons = xr.DataArray(df[lon_col].values, dims=dim, coords={dim: df[dim]})
    
    coords = dict(zip(get_coord_names(ds), [lats, lons]))
    
    ds_sub = ds.sel(coords, method='nearest')

    return ds_sub


def extract_da_value(da, lat, lon, variables=None):
    # TODO: Fix this. Rough for now to work with get_xr_values()
    """
    Get values from a DataArray at a given lat, lon
    """
    coords = dict(zip(get_coord_names(da), [lat, lon]))
    df = da.sel(
        **coords, method='nearest'
    ).to_dataframe()
    return df

def extract_xr_value(ds, lat, lon, variables=['texture', 'porosity']):
    """
    Get values from a dataset at a given lat, lon
    """
    coords = dict(zip(get_coord_names(ds), [lat, lon]))
    values = {var : ds.sel(
        **coords, method='nearest'
    )[var].item() for var in variables}
    return values

def extract_xr_values(ds, lat, lon, variables=['pet']):
    coords = dict(zip(get_coord_names(ds), [lat, lon]))
    values = ds.sel(
        **coords, method='nearest'
    ).compute().to_dataframe()#.reset_index()
    return values


def get_xr_values(
        df, ds, variables=['texture', 'porosity'], 
        lat_col='LAT', lon_col='LON', func=extract_xr_values
    ):
    """
    Get soil values from a dataset at a given lat, lon
    """
    values = df.apply(
        lambda row: pd.Series(
            func(ds, row[lat_col], row[lon_col], variables)
        ), axis=1
    )
    return values

def create_dir(out_dir):
    # out_dir = os.path.join(parent_dir, subdir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"Created directory: {out_dir}")
    else:
        print(f"Directory already exists: {out_dir}")

# def get_out_filename(row, suff):
#     info = row.to_dict()
#     try:
#         out_file =  f"{row['network']}_{row['station']}{suff}.csv"
#     except:



def process_subset(row, df, suff, out_dir, sub_cols=['network','station'], overwrite=False):
    info = row.to_dict()

    # out_file = get_out_filename(row, suff) 
    out_file = f"{info['network']}_{info['station']}{suff}.csv"

    if os.path.exists(os.path.join(out_dir, out_file)) and not overwrite:
        print(f"File {out_file} already exists. Skipping.")
    else:
        mask = (df[sub_cols] == row[sub_cols]).all(axis=1)
        df_ns = df[mask]
        if not df_ns.empty:
            df_ns.reset_index(drop=True, inplace=True)
            df_ns.to_csv(os.path.join(out_dir, out_file), index=False)
            print(f"Saved {out_file}")


def separate_df(
    df, by=['network', 'station'], threads=20, 
    suff=None, out_dir=None, sub_cols=['network','station'], overwrite=False #parent_dir=None
):
    if isinstance(by, list):
        stations = df[by].drop_duplicates()
    elif isinstance(by, pd.DataFrame):
        stations = by
    
    if threads:
        with mp.Pool(threads) as pool:
            pool.starmap(
                process_subset, 
                [(row, df, suff, out_dir, sub_cols, overwrite) for i,row in stations.iterrows()]
            )
        pool.close()
        pool.join()
    else:
        for i,row in stations.iterrows():
            process_subset(row, df, suff, out_dir, sub_cols) #parent_dir)

def expand_col(df, col='location', cols=['network','station'], drop=True):
    in_cols = list(df.columns)
    df[col] = df[col].apply(ast.literal_eval)
    df[cols] = pd.DataFrame(df[col].to_list(), index=df.index)
    # Reorder columns
    for i, c in enumerate(cols):
        in_cols.insert(i+1, c)

    df = df[in_cols]
    if drop:
        # df.drop(col, axis=1, inplace=True)
        del df[col]
    return df

def agg_data(
    df, col='pet', roll_args={'window':60, 'min_periods':60}, 
    names=['total60d','total_month'], mask=100., recalc=False
):
    df.time = pd.to_datetime(df.time, format='mixed') #format='%Y-%m-%d')
    df.set_index('time', inplace=True)
    df.sort_index(inplace=True)

    if mask:
        df.insert(df.columns.get_loc(col), col+'_raw', df[col].copy())
        df[col] = df[col].where(df[col] < mask)
        #df.insert(df.columns.get_loc(col)+1, col, df[col].where(df[col] < mask))

    if not names[0] in df.columns or recalc:
        # Calculate rolling 60-day sum
        df[names[0]] = df[col].rolling(**roll_args).sum()
    if not names[1] in df.columns or recalc:
        # Calculate monthly sum
        df[names[1]] = df[col].resample('MS').sum()
        df[names[1]] = df[names[1]].ffill()

    return df


def agg_file(
    file, col='pet', roll_args={'window':60, 'min_periods':1}, 
    names=['total60d','total_month'], recalc=False, expand_kwargs=None
):
    df = pd.read_csv(file)
    if not recalc and all([n in df.columns for n in names]):
        print(f"File {os.path.basename(file)} already aggregated. Skipping.")
        # return df
        return

    if expand_kwargs:
        df = expand_col(df, **expand_kwargs)

    try:
        df = agg_data(df, col=col, roll_args=roll_args, names=names, recalc=recalc)
    except Exception as e:
        print(f"Error aggregating {os.path.basename(file)}: {e}")
        return

    df.to_csv(file)
    print(f"Saved {os.path.basename(file)}")