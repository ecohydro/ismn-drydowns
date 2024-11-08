#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Sun 07 Apr 24 21:51:54'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           calc_ismn_daily.py
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


#-------------------------------------------------------------------------------
# IMPORTS
#-------------------------------------------------------------------------------
import os
import datetime
# from drydowns.config import Config

import station
from station import config, META, Station
#-------------------------------------------------------------------------------
# VARIABLES
#-------------------------------------------------------------------------------
# cfg_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','flx_config.ini')
# config = Config(cfg_file)

proj_dir = station.proj_dir
ismn_dir = station.ismn_dir #config.config.get('PATHS','ismn_dir')
net_dir = os.path.join(ismn_dir, 'Networks')


date = datetime.datetime.now().strftime('%d%b').lower()
daily_dir = os.path.join(ismn_dir, '..', 'daily')


#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------
def main():

    df = META[['network','station']].drop_duplicates()

    
    print(f"Initializing {len(df)} stations...")
    df['object'] = df.apply(
        lambda row: Station(config=config.config['ISMN'], network=row.network, station=row.station), 
        axis=1
    )
    print(f"Initialized {len(df)} stations.")


    for i,stat in enumerate(df.object):
        if os.path.exists(os.path.join(stat.config.get('daily_dir'), stat._daily_file())):
            print(f"Daily file for station {stat.name} already exists. Moving to next file")
            continue
        stat.save(stat.daily, fid='DD_all', out_dir=daily_dir)
        print(f'Finished processing {i+1}/{len(df)} stations.')



if __name__ == "__main__":
    main()