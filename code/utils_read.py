#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Wed 25 Sep 24 12:08:56'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           utils_read.py
Compatibility:  Python 3.12.0

Description:    Code for reading and processing outputs from drydown model.

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
# import sys
import re
import numpy as np
import pandas as pd

from scipy import stats


import seaborn as sns

# from utils_figs import igbp_dict, var_dict, get_continuous_cmap
from utils_figs import bin_dict #, igbp_cats, cmap


from soil.station import config, META #, ismn_data
# import drydowns
from drydowns.config import Config

# import utils

# ANCILLARY DATA FUNCTIONS

run_dict = {
    'ISMN' : {
        'dir' : 'ISMN',
        'prefix' : 'ismn',
        # 'id_cols' : ['network', 'station'],
    },
    'FLUXNET' : {
        'dir' : 'FluxData',
        'prefix' : 'flx',
        # 'id_cols' : ['site_id'],
    }
}

anc_dict = {
    'LC' : 'MCD12Q1_061_LC',
    'VCF' : 'MOD44B_061_VCF',
}

lc_dict = {
    10: 'Cropland, rainfed',
    11: 'Cropland, rainfed / Herbaceous cover',
    12: 'Cropland, rainfed / Tree or shrub cover',
    20: 'Cropland, irrigated or post-flooding',
    30: 'Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous',
    40: 'Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)',
    50: 'Tree cover, broadleaved, evergreen, Closed to open (>15%)',
    60: 'Tree cover, broadleaved, deciduous, Closed to open (>15%)',
    61: 'Tree cover, broadleaved, deciduous, Closed (>40%)',
    62: 'Tree cover, broadleaved, deciduous, Open (15-40%)',
    70: 'Tree cover, needleleaved, evergreen, closed to open (>15%)',
    71: 'Tree cover, needleleaved, evergreen, closed (>40%)',
    72: 'Tree cover, needleleaved, evergreen, open (15-40%)',
    80: 'Tree cover, needleleaved, deciduous, closed to open (>15%)',
    81: 'Tree cover, needleleaved, deciduous, closed (>40%)',
    82: 'Tree cover, needleleaved, deciduous, open (15-40%)',
    90: 'Tree cover, mixed leaf type (broadleaved and needleleaved)',
    100: 'Mosaic tree and shrub (>50%) / herbaceous cover (<50%)',
    110: 'Mosaic herbaceous cover (>50%) / tree and shrub (<50%)',
    120: 'Shrubland',
    121: 'Shrubland / Evergreen Shrubland',
    122: 'Shrubland / Deciduous Shrubland',
    130: 'Grassland',
    140: 'Lichens and mosses',
    150: 'Sparse vegetation (tree, shrub, herbaceous cover) (<15%)',
    152: 'Sparse vegetation (tree, shrub, herbaceous cover) (<15%) / Sparse shrub (<15%)',
    153: 'Sparse vegetation (tree, shrub, herbaceous cover) (<15%) / Sparse herbaceous cover (<15%)',
    160: 'Tree cover, flooded, fresh or brakish water',
    170: 'Tree cover, flooded, saline water',
    180: 'Shrub or herbaceous cover, flooded, fresh/saline/brakish water',
    190: 'Urban areas',
    200: 'Bare areas',
    201: 'Consolidated bare areas',
    202: 'Unconsolidated bare areas',
    210: 'Water',
    220: 'Permanent snow and ice'
 }


def read_lc(run_type, anc_type='LC'):
    # Read file
    lc = pd.read_csv(os.path.join(
        data_dir, run_dict[run_type].get('dir', run_type), 'ancillary', anc_type, 
        f'{run_dict[run_type]["prefix"]}_{anc_dict[anc_type]}.csv'
    ))
    # Add date info
    lc['Date'] = pd.to_datetime(lc.Date)
    lc['year'] = lc.Date.dt.year
    return lc


# DATA
# run_type = 'ISMN'

config = Config('soil/config.ini')
proj_dir = config.config['ISMN']['project_dir']    
# data_dir = os.path.join(proj_dir, 'data')
data_dir = os.path.join('..','data')
res_dir = os.path.join(proj_dir, 'outputs')

# chirps_dir = os.path.join(res_dir, 'chirps')

# Ancillary data
# IGBP/CCI Landcover codes
cci = pd.read_csv(os.path.join(data_dir,'cci_landcover.csv'))    # CCI Landcover
# mod_igbp = pd.read_csv(os.path.join(data_dir,'modis_landcover','IGBP.csv')) # IGBP Landcover
mod_igbp = pd.read_csv(os.path.join(data_dir,'IGBP.csv')) # IGBP Landcover

# Ancillary info for ISMN stations
ismn_anc = pd.read_csv(os.path.join(data_dir, 'ISMN', 'ismn_ancillary.csv')) # Rainfall for ISMN stations

# lc_dict = {v:k for k,v in ismn_data.landcover.items()}

ismn = META.drop_duplicates(subset=['network','station']).copy()
ismn = ismn.merge(cci, left_on='lc_2010', right_on='CCI_Code')


# for run in ['ISMN', 'FLUXNET']:
for run in ['ISMN']:
    run_dict[run].update({
        'lc' : read_lc(run, 'LC'),
        'vcf' : read_lc(run, 'VCF'),
    })

# ismn_lc = run_dict['ISMN']['lc']

# FUNCTIONS (IMPORT)

def read_output(
        file, 
        run_type = 'ISMN',
        array_cols=['theta', 'time', 'exp_theta_opt', 'q_theta_opt', 'et', 'gpp'], 
        dt_cols=['event_start', 'event_end'],
        col_dict=None,
        **kwargs
    ):
    # Read file
    df = pd.read_csv(file, **kwargs)
    # Rename columns
    if col_dict is not None:
        df.rename(columns=col_dict, inplace=True)
    # Convert timeseries to arrays
    try:
        print(f'Converting arrays in {os.path.basename(file)}')
        df = convert_arrays(df, array_cols=array_cols)
    except:
        pass
    # Convert datetime columns
    for col in dt_cols:
        if col in df.columns:
            df[col] = pd.to_datetime(df[col])
    # Add theta_w if available
    if 'end_dsdt_mm' in df.columns:
        get_theta_w(df)
    # Add ancillary info
    # if 'LC' in df.columns:
    df = add_info(df, run_type=run_type)
    add_season(df)
    print(f'Read {df.shape[0]} events from {os.path.basename(file)}')
    return df


def convert_arrays(df, array_cols=['theta', 'time', 'exp_theta_opt', 'q_theta_opt', 'gpp']):
    for col in array_cols:
        if col in df.columns:
            df.loc[:,col] = df.loc[:,col].apply(
                # lambda x: np.array(x.strip('[]').replace(',', '').split(), dtype=float)
                str_to_array
            )
    calc_ranges(df)
    return df

def str_to_array(x):
    try:
        arr = np.array(x.strip('[]').replace(',', '').split(), dtype=float)
        return arr
    except:
        return x

def calc_ranges(df):
    df['event_min'] = df.theta.apply(lambda x: x.min() if isinstance(x, np.ndarray) else np.nan)
    df['event_max'] = df.theta.apply(lambda x: x.max() if isinstance(x, np.ndarray) else np.nan)
    df['event_range'] = df.event_max - df.event_min
    df['frac_range'] = df.event_range/(df.max_sm-df.min_sm)
    # return df

def get_theta_w(df):
    df['num_theta_w'] = np.nan
    df.loc[df.end_dsdt_mm < 0., 'num_theta_w'] = df.event_min[df.end_dsdt_mm < 0.]


def add_info(df, run_type='ISMN'):
    # Rename columns
    df.rename(columns={'SITE_ID':'station', 'Sensor_ID':'sensor'}, inplace=True)

    # Add metadata (for ISMN stations only, for now)
    if run_type == 'ISMN':
        df = df.merge(
            META[['station','sensor', 'network','climate_KG','elevation']], 
            # left_on=['SITE_ID','Sensor_ID'], right_on=['station','sensor']
        )
        print(f'Adding ancillary data to {df.shape[0]} events')
        df = df.merge(ismn_anc)#[['network','station','MAP']])
        # Add CCI landcover codes
        df = df.merge(cci, left_on='LC', right_on='CCI_Code')
    
    if 'network' not in df.columns:
        df['network'] = run_type
        df['location'] = df.station.copy()

    # Add landcover data
    print(f'Adding landcover data to {df.shape[0]} events')
    # IGBP
    df = add_modis_lc(df, lc=run_dict[run_type].get('lc', None))
    # Vegetation cover fraction
    df = add_modis_vcf(df, vcf=run_dict[run_type].get('vcf', None))

    # Add rainfall data
    print(f'Adding rainfall data to {df.shape[0]} events')
    df = add_rainfall(df, run_type=run_type, how='GPP_mean')

    # Calculate AI + bin data
    if 'AI_60d' in df:
        df['AI_rel'] = df.AI_60d/df.AI
        df['AIi'] = 1/df.AI
        df['AIi_60d'] = 1/df.AI_60d
        df['AIi_rel'] = df.AIi_60d/df.AIi
        # Calculate z-scores
        df['z_LAI'] = stats.zscore(df.LAI, nan_policy='omit')
        df['z_AI_60d'] = stats.zscore(df.AI_60d, nan_policy='omit')
        print(f'Binning data')
        bin_data(df)
    
    # if 'GPP' in df.columns:
    gpp_col = next((col for col in ['GPP', 'total_gpp'] if col in df.columns), None)
    if gpp_col is not None:
        if run_type == 'FLUXNET':
            df[gpp_col] = df[gpp_col] / 1000.     # TODO: Change this elsewhere. This is temporary to match MODIS units.
        df['GPP_avg'] = df[gpp_col]/df.duration

    return df

forr_codes = {
    'ENF' : 1, 'EBF' : 2, 'DNF' : 3, 'DBF' : 4, 'MF' : 5,
}

def add_modis_lc(df, lc=None):
    if lc is None:
        print('No landcover data found. Skipping...')
        return df
    df['year'] = df.event_start.dt.year
    df = df.merge(
        lc[[
            'network', 'station', 'year', 
            'LC_IGBP', 'LC_UMD', 'LC_LAI', 'LC_BGC', 'LC_PFT',
            'FAO_LCCS1_LC', 'FAO_LCCS2_LU'
        ]], 
    )
    # Calculate min in category
    df['LC_IGBP_min'] = df.groupby('location')['LC_IGBP'].transform(lambda x: x.min())
    # Calculate mode in category
    df['LC_IGBP_mode'] = df.groupby('location')['LC_IGBP'].transform(lambda x: x.mode()[0])
    
    # Columns with codes/classes
    codes = ['LC_IGBP', 'LC_IGBP_min', 'LC_IGBP_mode']
    cols = ['MOD_IGBP', 'iIGBP', 'mIGBP']

    # Edit
    if 'CCI_IGBP_ed' in df.columns:
        df['LC_IGBP_ed'] = df.LC_IGBP.copy()
        mask = (df.LC_IGBP == 8.) & (df.CCI_IGBP_ed == 'FOR')
        df.loc[mask, 'LC_IGBP_ed'] = df.loc[mask, 'CCI_IGBP'].map(forr_codes)

        df['LC_IGBP_ed_min'] = df.groupby('location')['LC_IGBP_ed'].transform(lambda x: x.min())

        codes += ['LC_IGBP_ed', 'LC_IGBP_ed_min']
        cols += ['MOD_IGBP_ed', 'iIGBP_ed']
    
    # Add IGBP classes
    for code,col in zip(codes, cols):
        df = get_igbp_classes(df, code=code, col=col)

    return df

def add_modis_vcf(df, vcf=None):
    if vcf is None:
        print('No vegetation cover fraction data found. Skipping...')
        return df
    df = df.merge(
        vcf[[
            'network', 'station', 'year', 'perc_nonwoodyveg', 
            'perc_woodyveg', 'perc_nonveg', 'perc_woodyveg_sd', 'perc_nonveg_sd'
        ]],
    )
    df.loc[df.perc_woodyveg > 100., 'perc_woodyveg'] = np.nan
    df.loc[df.perc_nonwoodyveg > 100., 'perc_nonwoodyveg'] = np.nan
    df.loc[df.perc_nonveg > 100., 'perc_nonveg'] = np.nan
    df.perc_woodyveg_sd = df.perc_woodyveg_sd/100.
    df.perc_nonveg_sd = df.perc_nonveg_sd/100.

    df['perc_woody'] = (df.perc_woodyveg / (100-df.perc_nonveg)) * 100
    df['iperc_woodyveg'] = df.groupby('location')['perc_woodyveg'].transform('max')
    df['iperc_woodyveg_avg'] = df.groupby('location')['perc_woodyveg'].transform('mean')

    return df

def get_igbp_classes(df, code='LC_IGBP', col='MOD_IGBP'):
    match = re.search(f"(?<={re.escape('LC_IGBP')}).*", code)
    suff = match.group(0) if match else '' #'_og'
    df = df.merge(
        mod_igbp[['IGBP_Code','MOD_IGBP']],
        left_on=code, right_on='IGBP_Code',
        suffixes=('',suff)
    )
    # if code == 'LC_IGBP':
    df.rename(columns={f'MOD_IGBP{suff}': col}, inplace=True)
        # suff = '_og'
    # Reclassify into combined categories
    reclass_all(df, col)
    return df

def reclass_all(df, col):
    reclass_igbp(df, new=f'{col}_5', col=col, sav=True, blf=True)   # combine savannas and BLFs
    reclass_igbp(df, new=f'{col}_6s', col=col, sav=False, blf=True) # combine BLFs
    reclass_igbp(df, new=f'{col}_6b', col=col, sav=True, blf=False) # combine savannas
    reclass_igbp(df, new=f'{col}_7', col=col, sav=False, blf=False) # leave separate
    reclass_igbp(df, new=f'{col}_2', col=col, sav=False, blf=False, wood=True) # combine woody vs. herbaceous
    reclass_igbp(df, new=f'{col}_3', col=col, sav=True, blf=False, wood=True) # combine woody vs. herbaceous

def reclass_igbp(df, new, col, sav=True, blf=False, wood=False):
    df[new] = df[col].copy()
    # 1. Combine shrublands
    df.loc[(df[new].isin(['OSH','CSH'])), new] = 'SHR'
    # 2. 
    if wood:
        if sav:
            df.loc[df[new].isin(['GRA']), new] = 'HERB'
            df.loc[df[new].isin(['SAV']), new] = 'SAV'
        else:
            df.loc[df[new].isin(['GRA','SHR','SAV']), new] = 'HERB'
        df.loc[df[new].isin(['WSA','EBF','DBF','ENF','MF','BLF']), new] = 'WOOD'

    # 2. Optionally combine savannas
    elif sav:
        df.loc[df[new] == 'WSA', new] = 'SAV'
    # 3. Optionally combine broadleaf forests
    if blf:
        df.loc[(df[new].isin(['EBF','DBF','MF'])), new] = 'BLF'
    # 4. Alternatively, combine into woody vs. herbaceous


def add_rainfall(df, run_type='ISMN', how='GPP_mean'):
    rain_file = os.path.join(res_dir, 'rain', f'{run_type}_rain_stats_p05.csv')
    seas_file = os.path.join(res_dir, 'rain', f'{run_type}_rain_season_{how}.csv')

    if os.path.exists(rain_file):
        df['month'] = df.event_start.dt.month
        rain = pd.read_csv(rain_file)
        rain.drop(columns=['latitude', 'longitude'], inplace=True) 
        df = df.merge(rain, on=['network','station','month'])

    if os.path.exists(seas_file):
        seas = pd.read_csv(seas_file)
        df = df.merge(seas, on=['network','station'], how='left', suffixes=('_monthly',''))
        df = convert_arrays(df, array_cols=['months'])
    
    return df

def bin_data(df):
    for key in bin_dict.keys():
        col = re.findall(r'(.*?)_bin', key)[0]
        df[key] = pd.cut(
            df[col], 
            bins=bin_dict.get(key)['bin'],
            labels=bin_dict.get(key)['bin_labels'],
            include_lowest=True, right=False
        )
        if key == 'MAP_bin':
            df[key+'_log'] = pd.cut(
                df[col], 
                bins=bin_dict.get(key)['bin_log'],
                labels=bin_dict.get(key)['bin_log_labels'],
                include_lowest=True, right=False
            )
        elif key in ['AI_60d_bin', 'AI_bin', 'AIi_60d_bin', 'AIi_bin', 'LAI_bin', 'perc_woodyveg_bin', 'perc_woody_bin']:
            df[col+'_bin2'] = pd.cut(
                df[col], 
                bins=bin_dict.get(key)['bin2'],
                labels=bin_dict.get(key)['bin2_labels'],
                include_lowest=True, right=False
            )
        elif col in ['frequency', 'intensity', 'coeff_var']:
            df[key+'_mth'] = pd.cut(
                df[col+'_monthly'], 
                bins=bin_dict.get(key)['bin_mth'],
                labels=bin_dict.get(key)['bin_mth_labels'],
                include_lowest=True, right=False
            )
        if key in ['perc_woodyveg_bin']:
            df[f'i{col}_bin'] = pd.cut(
                df[f'i{col}'], 
                bins=bin_dict.get(key)['bin'],
                labels=bin_dict.get(key)['bin_labels'],
                include_lowest=True, right=False
            )
            df[f'i{col}_avg_bin'] = pd.cut(
                df[f'i{col}_avg'], 
                bins=bin_dict.get(key)['bin'],
                labels=bin_dict.get(key)['bin_labels'],
                include_lowest=True, right=False
            )

    df['precip_bin'] = pd.cut(
        df.precip_total60d,
        bins=[0, 20, 50, 80, 100, 150, 250, 500, 3000],
        labels=['0-20','20-50','50-80','80-100','100-150','150-250','250-500','>500'],
        include_lowest=True, right=False
    )  

def add_season(df):
    if 'months' in df.columns:
        seas_mask = df.apply(lambda x : x.month in x.months if isinstance(x.months, np.ndarray) else False, axis=1)
        df['season'] = 'dry'
        df.loc[seas_mask, 'season'] = 'wet'


def filter_data(df, r2_thresh=0.7, duration : tuple=(0., 100.), rfrac : float = 0.2):
    filt = df[
        (df.q_r_squared >= r2_thresh)
        & (df.duration >= duration[0])
        & (df.duration <= duration[1])
        & (df.q_q > 0.01)
        & (df.frac_range > rfrac)
        & (df.z_m <= 0.5)
        & (df.q_theta_star - df.min_sm >= 0.005)
        # & (df.MAP <= 1200.)
    ]
    if 'q_theta_w' in df.columns:
        filt = filt[filt.q_theta_star - filt.q_theta_w >= 0.005]
    if 'dz' in df.columns:
        filt = filt[filt.dz <= 0.5]
    return filt


def clean_data(df, r2_thresh, duration : tuple=(0., 100.), rfrac : float = 0.2):
    df = filter_data(df, r2_thresh, duration, rfrac=rfrac)
    return df