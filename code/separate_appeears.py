#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Thu 18 Apr 24 13:57:00'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           separate_appeears.py
Compatibility:  Python 3.7.0
Description:    Extracts individual station data from the APPEEARS MOD15A2H output
                + saves to individual csv files.

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
import re
import glob
import numpy as np
import pandas as pd
import multiprocessing as mp
import time
import utils

from drydowns.config import Config

# VARIABLES
config = Config('soil/config.ini')
run_type = config.config.get('RUN', 'type')     # Can also edit manually ('ISMN' or 'FLUXNET')
# cfg = config[config.get('RUN','type')]
proj_dir = config.config.get(run_type, 'project_dir')
data_dir = os.path.join(proj_dir, 'data')
# anc_dir = config.config.get('ISMN', 'anc_dir')

anc_file = os.path.join(data_dir, 'ISMN', 'ismn_ancillary.csv')

if run_type == 'ISMN':
    coord_dir = os.path.join(data_dir, 'ISMN', 'coordinates')
else:
    coord_dir = os.path.join(data_dir, 'FluxData')


# VARIABLES
pre_dict = {
    'ISMN' : 'ismn',
    'FLUXNET' : 'flx',
    'SMAP' : 'smap',
}

suff_dict = {
    'LAI' : 'MOD15A2H',
    'GPP' : 'MOD17A2HGF',
    'LC' : 'MCD12Q1',
    'VCF' : 'MOD44B',
}

band_dict = {
    'MCD12Q1_061_LC_Prop1' : 'FAO_LCCS1_LC',
    'MCD12Q1_061_LC_Prop2' : 'FAO_LCCS2_LU',
    'MCD12Q1_061_LC_Prop3' : 'FAO_LCCS3_HYD',
    'MCD12Q1_061_LC_Prop1_Assessment' : 'FAO_LCCS1_confidence',
    'MCD12Q1_061_LC_Prop2_Assessment' : 'FAO_LCCS2_confidence',
    'MCD12Q1_061_LC_Prop3_Assessment' : 'FAO_LCCS3_confidence',
    'MCD12Q1_061_LC_Type1' : 'LC_IGBP',
    'MCD12Q1_061_LC_Type2' : 'LC_UMD',
    'MCD12Q1_061_LC_Type3' : 'LC_LAI',
    'MCD12Q1_061_LC_Type4' : 'LC_BGC',
    'MCD12Q1_061_LC_Type5' : 'LC_PFT',
    'MOD44B_061_Percent_NonTree_Vegetation' : 'perc_nonwoodyveg',
    'MOD44B_061_Percent_NonVegetated' : 'perc_nonveg',
    'MOD44B_061_Percent_Tree_Cover' : 'perc_woodyveg',
    'MOD44B_061_Percent_NonVegetated_SD' : 'perc_nonveg_sd',
    'MOD44B_061_Percent_Tree_Cover_SD' : 'perc_woodyveg_sd',
}


# FUNCTIONS
def rename_cols(
    df, 
    # id_cols=['station', 'network', 'latitude', 'longitude'], 
    id_cols=['station', 'network'], 
    exclude=r'MOD15A2H_061_'
):
    # Raw column names
    if 'network' not in df.columns:
        df.insert(1, 'network', run_type)
    raw_cols = df.columns
    # New column names
    cols = id_cols + ['latitude','longitude'] + [
        re.sub(re.compile(exclude), '', s) for s in raw_cols[len(id_cols)+2:]
    ]
    # Rename the columns
    df.rename(columns=dict(zip(raw_cols, cols)), inplace=True)
    # Reorder columns
    # df = df[['network','station'] + cols[2:]]
    df = df[id_cols[::-1] + cols[len(id_cols):]]
    return df


def fix_labels(df, coord_file):
    coords = pd.read_csv(coord_file)
    # Mask Antarctica site...
    coords = coords[coords.latitude > -60.].reset_index(drop=True)

    # Add network
    if 'network' not in coords.columns:
        coords = rename_cols(coords)

    df.sort_values(by=['network','station'], inplace=True)
    # df.reset_index(drop=True, inplace=True)

    dfi = df[['network','station','latitude','longitude']].drop_duplicates().reset_index(drop=True)
    

    df_map = coords.join(dfi, how='right', rsuffix='_y')
    # df_map.drop(columns=['latitude_y','longitude_y'], inplace=True)

    df = df_map.merge(
        df, left_on=['network_y','station_y', 'latitude_y','longitude_y'], 
        right_on=['network','station', 'latitude','longitude'], 
        suffixes=['','_y']
    )
    df.drop(columns=['network_y','station_y','latitude_y','longitude_y'], inplace=True)

    return df


# def separate_appeears_df(df, out_dir):
#     # Get the unique stations
#     stations = df[['network','station']].drop_duplicates()
#     # Separate out the dfs
#     for i,row in stations.iterrows():
#         label = row.to_dict()
#         out_file = f"{label['network']}_{label['station']}_MOD15A2H.csv"
#         if os.path.exists(os.path.join(out_dir, out_file)):
#             print(f"File {out_file} already exists. Skipping.")
#             continue
#         df_ns = df[
#             (df['network'] == label['network']) 
#             & (df['station'] == label['station'])
#         ]
#         df_ns.reset_index(drop=True, inplace=True)
#         df_ns.to_csv(os.path.join(out_dir, out_file), index=False)
#         print(f"Saved {out_file}")


# def process_appeears(file, coord_file, out_dir):
def process_appeears(file, coord_file):
    # Read in the data
    df = pd.read_csv(file)
    # Rename the columns
    df = rename_cols(df)
    # Fix labels
    df = fix_labels(df, coord_file)
    # # Separate dfs
    # separate_appeears_df(df, out_dir)
    print(f"Finished prepping {file}")
    return df


def separate_appeears(app_dir, coord_dir, out_dir, suff=f'_MOD15A2H', threads=40, overwrite=False):
    # Get the APPEEARS files
    files = sorted([f for f in glob.glob(os.path.join(app_dir, '*.csv'))])
    # Get the coordinate files
    coord_files = sorted([f for f in glob.glob(os.path.join(coord_dir, 'ismn_coords_*.csv'))])
    # Iterate through APPEEARS results files.
    for file, coord_file in zip(files, coord_files):
        start = time.perf_counter()

        df = process_appeears(file, coord_file)
        # Save to file
        utils.separate_df(
            df=df, by=['network','station'], 
            suff=suff, out_dir=out_dir,
            sub_cols = ['network', 'station'],
            overwrite=overwrite,
            threads = threads,
        )

        end = time.perf_counter()
        print(f"Finished processing in {end-start:.6f} seconds.")


def combine_appeears(var):
    # Get files
    out_dir = config.config.get(run_type, f'{var.lower()}_dir')

    # Read dataframes + concat
    if run_type == 'ISMN':
        files = sorted([f for f in glob.glob(os.path.join(out_dir, '*.csv'))])
        dfs = [pd.read_csv(f) for f in files]
        df = pd.concat(dfs, ignore_index=True)
    # Process if FLUXNET
    if run_type == 'FLUXNET' and var in ['LC', 'VCF']:
        app_dir = os.path.join(out_dir, 'APPEEARS')
        # Get file
        file = glob.glob(os.path.join(app_dir, '*.csv'))[0]
        # Process
        df = process_appeears(file, os.path.join(coord_dir, 'coords.csv'))

    df.rename(columns=band_dict, inplace=True)

    out_file = f'{pre_dict[run_type]}_{suff_dict[var]}_061_{var}.csv'
    print(f"Saving {var} data to {out_file}...")
    df.to_csv(os.path.join(os.path.dirname(out_dir), out_file), index=False)



# MAIN
def main():
    config = Config('soil/config.ini')

    var = 'LC'

    out_dir = config.config.get(run_type, f'{var.lower()}_dir')
    suff = f'_{suff_dict[var]}'


    # lai_dir = config.config.get('ISMN', 'lai_dir')
    # gpp_dir = config.config.get('ISMN', 'gpp_dir')
    # lc_dir = config.config.get('ISMN', 'lc_dir')
    # vcf_dir = config.config.get('ISMN', 'vcf_dir')

    # # out_dir = lai_dir
    # # suff = '_MOD15A2H'

    # # out_dir = gpp_dir
    # # suff = '_MOD17A2HGF'

    # # out_dir = lc_dir
    # # suff = '_MCD12Q1'

    # out_dir = vcf_dir
    # suff = '_MOD44B'

    threads = 40
    overwrite = True

    separate_appeears(
        app_dir=os.path.join(out_dir,'APPEEARS'), 
        coord_dir=coord_dir, 
        out_dir=out_dir,
        suff=suff, 
        threads=threads, 
        overwrite=overwrite
    )

    # files = sorted([f for f in glob.glob(os.path.join(lai_dir,'APPEEARS', '*.csv'))])
    # coord_files = sorted([f for f in glob.glob(os.path.join(coord_dir, 'ismn_coords_*.csv'))])
    # # files.sort()
    
    # for file,coord_file in zip(files, coord_files):
    #     start = time.perf_counter()

    #     df = process_appeears(file, coord_file)

    #     # Save to file
    #     utils.separate_df(
    #         df=df, by=['network','station'], 
    #         suff=f'_MOD15A2H', out_dir=out_dir,
    #         sub_cols = ['network', 'station'],
    #         overwrite=True,
    #         threads = 40,
    #     )

    #     end = time.perf_counter()
    #     print(f"Finished processing in {end-start:.6f} seconds.")



#
if __name__ == "__main__":
    main()

# For LC + VCF, after separating, combine the files into one.

var = 'LC'

combine_appeears(var)



# band_dict = {
#     'MCD12Q1_061_LC_Prop1' : 'FAO_LCCS1_LC',
#     'MCD12Q1_061_LC_Prop2' : 'FAO_LCCS2_LU',
#     'MCD12Q1_061_LC_Prop3' : 'FAO_LCCS3_HYD',
#     'MCD12Q1_061_LC_Prop1_Assessment' : 'FAO_LCCS1_confidence',
#     'MCD12Q1_061_LC_Prop2_Assessment' : 'FAO_LCCS2_confidence',
#     'MCD12Q1_061_LC_Prop3_Assessment' : 'FAO_LCCS3_confidence',
#     'MCD12Q1_061_LC_Type1' : 'LC_IGBP',
#     'MCD12Q1_061_LC_Type2' : 'LC_UMD',
#     'MCD12Q1_061_LC_Type3' : 'LC_LAI',
#     'MCD12Q1_061_LC_Type4' : 'LC_BGC',
#     'MCD12Q1_061_LC_Type5' : 'LC_PFT',
# }

# files = sorted([f for f in glob.glob(os.path.join(lc_dir, '*.csv'))])

# dfs = [pd.read_csv(f) for f in files]

# df = pd.concat(dfs, ignore_index=True)

# df.rename(columns=band_dict, inplace=True)

# df.to_csv(os.path.join(lc_dir, '..', 'ismn_MCD12Q1_061_LC.csv'), index=False)


# band_dict = {
#     'MOD44B_061_Percent_NonTree_Vegetation' : 'perc_nonwoodyveg',
#     'MOD44B_061_Percent_NonVegetated' : 'perc_nonveg',
#     'MOD44B_061_Percent_Tree_Cover' : 'perc_woodyveg',
#     'MOD44B_061_Percent_NonVegetated_SD' : 'perc_nonveg_sd',
#     'MOD44B_061_Percent_Tree_Cover_SD' : 'perc_woodyveg_sd',
# }


# files = sorted([f for f in glob.glob(os.path.join(vcf_dir, '*.csv'))])
# dfs = [pd.read_csv(f) for f in files]
# df = pd.concat(dfs, ignore_index=True)
# df.rename(columns=band_dict, inplace=True)

# df.to_csv(os.path.join(vcf_dir, '..', 'ismn_MOD44B_061_VCF.csv'), index=False)