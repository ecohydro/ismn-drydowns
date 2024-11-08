#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Sun 07 Apr 24 21:28:59'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           station.py
Compatibility:  Python 3.10+
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
import re

import numpy as np
import pandas as pd


import timezonefinder as tzf



import ismn
import drydowns
from drydowns.config import Config 

from collections import OrderedDict
# from ismn.components import Station, Sensor
from ismn import components
from ismn.interface import ISMN_Interface

#-------------------------------------------------------------------------------
# VARIABLES
#-------------------------------------------------------------------------------
cfg_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'config.ini')
config = Config(cfg_file)

proj_dir = config.config.get('ISMN','project_dir')
# ismn_dir = config.config.get('ISMN','ismn_dir')
ismn_dir = config.config.get('ISMN','data_dir')
# net_dir = os.path.join(ismn_dir, 'Networks')
anc_dir = config.config.get('ISMN', 'anc_dir')

# TODO: First, change Sensor() to not need ismn_data objects.
# TODO: Remove init of ISMN_Interface (it takes too long). (Currently initing a sensor
# has a dependency on this, so first need to remove that).
# ismn_data = ISMN_Interface(ismn_dir, parallel=False) # for some reason, parallel isn't working

def load_meta(how='file', dt_cols=['timerange_to','timerange_from'], td_cols=['avg_data_freq']):
    # Load metadata from file
    meta_file = os.path.join(ismn_dir, '..', 'metadata.csv')
    if os.path.exists(meta_file) and how == 'file':
        META = load_meta_from_file(meta_file, dt_cols=dt_cols, td_cols=td_cols)
    elif how == 'file':
        meta_file = os.path.join('..', 'data', 'ISMN', 'metadata.csv')
        META = load_meta_from_file(meta_file, dt_cols=dt_cols, td_cols=td_cols)
    else:
        ismn_data = ISMN_Interface(ismn_dir, parallel=False) # for some reason, parallel isn't working
        META = ismn_data.metadata.copy()
        # Restructure metadata to get rid of unnecessary depth columns
        META[[('depth_from','val'),('depth_to','val')]] = META[
            [('variable','depth_from'),('variable','depth_to')]
        ]
        META[[('text_depth_from','val'),('text_depth_to','val')]] = META[
            [('clay_fraction','depth_from'),('clay_fraction','depth_to')]
        ]
        META = META.loc[:,META.columns.get_level_values(1)=='val']
        META.columns = META.columns.droplevel(1)
        META['extent'] = META.apply(lambda x : (x.depth_from, x.depth_to), axis=1)
        META['sensor'] = [s[2].name for s in ismn_data.collection.iter_sensors()]
        META['soil_texture'] = drydowns.soil.get_soil_texture(
            sand = META.sand_fraction, silt = META.silt_fraction, clay = META.clay_fraction
        )
    return META

def load_meta_from_file(file, dt_cols=['timerange_to','timerange_from'], td_cols=['avg_data_freq']):
    META = pd.read_csv(file)
    for col in dt_cols:
        META[col] = pd.to_datetime(META[col])
    for td_cols in td_cols:
        META[td_cols] = pd.to_timedelta(META[td_cols])
    return META

META = load_meta(how='file', dt_cols=['timerange_to','timerange_from'])

# META['lower'] = np.nan
# META['dz'] = np.nan

# If you add or remove folders from the collection, delete the python_metadata folder + regenerate it.
# ismn_data.networks
# Station variables: climate_KG, elevation, latitude, lc_2000, lc_2005, lc_2010, longitude, network, station, timerange_from (?), timerange_to (?)
# Sensor variables: instrument, clay_fraction, sand_fraction, silt_fraction, organic_carbon, depth_from, depth_to, timerange_from, timerange_to, name
# Instrument variables: instrument


# Variables:
# ['soil_temperature', 'soil_moisture', 'precipitation', 'air_temperature', 
#  'snow_depth', 'surface_temperature', 'snow_water_equivalent', 'soil_suction']


var_dict = {
    'soil_temperature': 'T_soil',
    'soil_moisture': 'SWC',
    'precipitation': 'P',
    'air_temperature': 'T_a',
    'snow_depth': 'Snow',
    'surface_temperature': 'T_surf',
    'snow_water_equivalent': 'SWE',
    'soil_suction': 'Psi_s',
    'pet' : 'PET',
    'Lai_500m' : 'LAI',
    # 'precip'
}

# TODO: combine with above.
# TODO: Add all flags. These are just the ones I'm considering "bad" for now. And 
# this is unfinished...
flag_dict = {
    'T_soil' : ['C01', 'C02'],
    'SWC' : ['C01', 'C03'],
    'P' : ['C01', 'C02'],
    'T_a' : ['C01', 'C02'],
    'Snow' : ['C01', 'C02'],
    'T_surf' : ['C01'],
    'SWE' : ['C01', 'C02'],
    'Psi_s' : ['C01', 'C02']
}

anc_dict = {
    # 'dPET' : ('dPET', 'pet'),
    # 'CHIRPS' : ('chirps05', 'precip'),
    # directory, suffix, column
    'PET' : ('dPET', 'dPET', 'pet'),
    'P' : ('CHIRPS', 'chirps05', 'precip'),
    'LAI' : ('LAI/MOD15A2H-061', 'MOD15A2H', 'Lai_500m')
}

# anc_dirs = {v[0] : k for k,v in anc_dict.items()}

# var_df = pd.DataFrame({
#     'variable' : [
#         'soil_temperature', 'soil_moisture', 'precipitation', 'air_temperature', 
#         'snow_depth', 'surface_temperature', 'snow_water_equivalent', 'soil_suction'
#     ],
#     'col' : ['T_soil', 'SWC', 'P', 'T_a', 'Snow', 'T_surf', 'SWE', 'Psi_s'],
#     # 'unit' : []
# })

#-------------------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------------------
def get_cols(col_list, variable=r'SWC', excludes=['_flag']):
    if variable:
        excludes_pattern = '|'.join(map(re.escape, excludes))
        pattern = re.compile(r'{}_(?!.*(?:{}))'.format(variable, excludes_pattern))
    else:
        excludes_pattern = '|'.join(map(re.escape, excludes))
        pattern = re.compile(f"^(?!.*(?:{excludes_pattern})).*$")
    cols = list(filter(pattern.match, col_list))
    return cols

def add_qc_flag(df, variable=None, excludes=['_flag'], suffix='_flag', cond='G'):
    cols = get_cols(df.columns, variable=variable, excludes=excludes+['_masked'])
    col_list = get_cols(df.columns, variable=variable, excludes=excludes)
    # for col in cols:
    mask = df[[col + suffix for col in cols]] == cond
    mask.rename(columns={col+suffix : col+'_QC' for col in cols}, inplace=True)

    return pd.concat([df[col_list], mask], axis=1)

def set_col_tz(df, tz_name, inplace=True):
    if not df.index.tz:
        utc = df.index.tz_localize('UTC')
        local = utc.tz_convert(tz_name)
        return df.set_index(local, inplace=inplace)
    else:
        return df

def check_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"Created directory: {out_dir}")
    else:
        print(f"Directory already exists: {out_dir}")

def calc_dz(depths : pd.Series, how='mid', lower_bound='last'):
    depths = depths.sort_values()
    mid = ((depths + depths.shift(-1))/2).reset_index(drop=True)
    # final = depths.iloc[-1]
    # Set lower bound to last value
    if lower_bound == 'last':
        final = depths.iloc[-1]
    # Set lower bound to dz_n/2 + last value
    elif not lower_bound:
        final = depths.iloc[-1] + (depths.iloc[-1]-mid.iloc[-2])
    # set lower bound to passed value
    else:
        final = lower_bound
    # Get lower bounds
    lower =  pd.Series(mid.to_list()[:-1]+[final]).round(6)
    # Calculate lower bounds
    dz = lower.diff().fillna(mid).round(6)
    # Convert to dicts
    lower_dict = dict(zip(depths, lower))
    dz_dict = dict(zip(depths, dz))
    return lower_dict, dz_dict

#-------------------------------------------------------------------------------
# CLASSES
#-------------------------------------------------------------------------------

class Station(components.Station):

    ismn_data = ISMN_Interface(ismn_dir, parallel=False) # for some reason, parallel isn't working

    station_vars = [
        'station', 'network', 'latitude', 'longitude', 'elevation', 
        'climate_KG', 'lc_2000', 'lc_2005', 'lc_2010',
    ]
    def __init__(self, config, station, network, keep_data=False):
        self.config = config
        self.id = station
        self.network = network
        self.var_info = META[(META.network == network) & (META.station == station)].sort_values(by='sensor')
        self.meta = self._get_meta()
        self.lon, self.lat, self.elev = self.meta['longitude'],self.meta['latitude'],self.meta['elevation']

        self._keep_data = keep_data

        super().__init__(name=self.id, lon=self.lon, lat=self.lat, elev=self.elev)

        self.sensors = self._get_sensors()
        self.variables = self.get_variables()
        self._var_instruments = self._get_var_instruments()
        self.instruments = self._get_instruments()

        self._data = None
        self._sensor_cols = None
        self._daily = None
        self.fid = f'{self.network}_{self.id}'
        # self.sensors = station.sensors  # This is messy, but allows access to everything else that's needed

        # self.sensors = OrderedDict([
        #     (k, Sensor(station.sensors[k])) for k in sorted(station.sensors)
        # ])

        # self.instruments = sorted(list(set(self.metadata['instrument'].values())))
        # self._set_tz_name()

    # def create_sensor(self, sensor):
    #     return Sensor(sensor)

    def _get_meta(self):
        meta = self.var_info[self.station_vars].iloc[0].to_dict()
        # meta['instruments'] = self.instruments
        # meta['variables'] = self.variables
        return meta
    
    # @property
    # def instruments(self):
    def _get_instruments(self):
        # names = sorted(list(set(self.metadata['instrument'].values())))
        # instruments = OrderedDict([])
        instruments = {}
        for var in sorted(self._var_instruments.keys()):
            instruments.update(self.add_instruments(var))
        # instruments = self.add_instruments()
        instruments = OrderedDict([
            (k, instruments[k]) for k in sorted(self.var_info.instrument.unique())
        ])

        return instruments

    @property
    def data(self):
        # if self._data is None:
        #     instrument_data = [inst.data for inst in self.iter_instruments()]
        #     self._data = pd.concat(instrument_data, axis=1)
        # return self._data
        return self.read_data()

    # @data.setter
    def read_data(self):
        if self._data is None:
            instrument_data = [inst.data for inst in self.iter_instruments()]
            data = pd.concat(instrument_data, axis=1)
            if self._keep_data:
                self._data = data
            else:
                return data
        # else:
        return self._data
    # def agg_data()


    def _get_sensors(self):
        # TODO: Don't use existing Station class; init sensors manually
        sensors = self.ismn_data[self.network][self.id].sensors
        # sensor_init_list=['instrument', 'variable', 'depth', 'filehandler','keep_loaded_data', 'name']
        # kwargs_dict = {
        #     s : {k:v for k,v in sensor.__dict__.items() if k in sensor_init_list} for s in sorted(sensors)
        # }
        # sensors = OrderedDict([
        #     (k, Sensor(**kwargs)) for k,kwargs in kwargs_dict.items()
        # ])
        sensors = OrderedDict([
            (k, self.add_sensor(sensors[k])) for k in sorted(sensors)
        ])
        return sensors
    
    # def add_sensor(self, instrument, variable, depth, filehandler=None, name=None, keep_loaded_data=False):
    #     return Sensor(
    #         instrument=instrument, variable=variable, depth=depth, filehandler=filehandler, 
    #         keep_loaded_data=keep_loaded_data, name=name
    #     )
    def add_sensor(self, sensor : components.Sensor):
        sensor_init_list=['instrument', 'variable', 'depth', 'filehandler','keep_loaded_data', 'name']
        kwargs = {'station' : self.id} | {k:v for k,v in sensor.__dict__.items() if k in sensor_init_list}
        return Sensor(**kwargs)

    # @property
    def _get_var_instruments(self):
        vi = self.var_info.groupby('variable')['instrument'].apply(
                lambda x: sorted(x.unique().tolist())
            ).to_dict()
        return vi

    def add_instruments(self, variable=None):
        inames = self._var_instruments.get(variable)
        # if variable:
        #     inames = self._var_instruments.get(variable)
        # else:
        #     inames = sorted(self.var_info.instrument.unique())
        instruments = OrderedDict([(
                name, self.add_instrument(instrument_name=name, instrument_label=i+1)
            ) for i,name in enumerate(inames)
        ])
        return instruments

    def add_instrument(self, instrument_name, instrument_label):
        sensors = self.get_isensors(instrument_name=instrument_name)
        return Instrument(
                name=instrument_name, station=self.id, 
                instrument_label=instrument_label, sensors=sensors
            )

    def get_isensors(self, instrument_name):
        sensors = [i for i in self.iter_sensors(instrument=instrument_name)]
        return sensors
    
    @property
    def sensor_cols(self):
        return self._set_sensor_cols()
    
    # @sensor_cols.setter
    def _set_sensor_cols(self):
        if self._sensor_cols is None:
            self._sensor_cols = self.get_sensor_cols()
        return self._sensor_cols
    
    def get_sensor_cols(self, **filter_kwargs):
        sensor_cols = {}
        for i in self.iter_instruments(**filter_kwargs):
            sensor_cols.update(i.get_sensor_cols(**filter_kwargs))
        return sensor_cols


    #TODO: Refactor with iter_instruments instead (probably slower, but matches 
    # structure of ismn code).

    def _set_tz_name(self):
        
        self._tz_name = tzf.TimezoneFinder().timezone_at(
            lng=self.lon, lat=self.lat
        )
        # return tz_name
    
    
    def iter_instruments(self, **filter_kwargs):
        for instrument in self.instruments.values():
            for sensor in instrument.sensors.values():
                if sensor.eval(**filter_kwargs):
                    yield instrument
                    break
    
    def calc_sm_dz(self, how='mid', lower_bound='last'):
        depths = self.var_info[
            self.var_info.variable == 'soil_moisture'
        ].sort_values(by='depth_to').drop_duplicates('depth_to').depth_to
        depths.reset_index(drop=True, inplace=True)
        lower_dict, dz_dict = calc_dz(depths, how=how, lower_bound=lower_bound)
        return lower_dict, dz_dict

    def add_sm_dz(self, how='mid', lower_bound='last'):
        if 'soil_moisture' in self.variables:
            if 'lower' not in META.columns:
                META['lower'] = np.nan
            if 'dz' not in META.columns:
                META['dz'] = np.nan
            lower, dz = self.calc_sm_dz(how=how, lower_bound=lower_bound)
            sm_view = META.loc[self.var_info.index][META.variable == 'soil_moisture']
            META.loc[sm_view.index,'lower'] = sm_view.depth_to.map(lower)
            META.loc[sm_view.index,'dz'] = sm_view.depth_to.map(dz)



    @property
    def daily(self):
        return self._set_daily()
    
    def _set_daily(self, set_tz=False, dropna=True, suff='DD_all'):
        if self._daily is None:
            # self.config.logger.info(f"Looking for existing file...")
            print(f"Looking for existing file...")
            if os.path.exists(os.path.join(self.config.get('daily_dir'), self._daily_file(suff))):
                # self.config.logger.info(f"Found daily file. Importing {self._daily_file(suff)}.")
                print(f"Found daily file. Importing {self._daily_file(suff)}.")
                daily = self.read_daily(suff=suff)
            else:
                daily = self.agg_daily(set_tz=set_tz, dropna=dropna)
            if self._keep_data:
                self._daily = daily
            return daily
        else:
            return self._daily

    def read_daily(self, suff='DD_all'):
        filepath = os.path.join(self.config.get('daily_dir'), self._daily_file(suff))
        daily = pd.read_csv(filepath, index_col=0)
        daily.index = pd.to_datetime(daily.index)
        return daily
    
    def _daily_file(self, suff='DD_all'):
        return f'{self.fid}_{suff}.csv'

    def agg_daily(self, set_tz=False, dropna=True):
        # 1. Add QC flag for aggregation.
        data = add_qc_flag(self.data, variable=None)
        # 2. Convert to local tz.
        if set_tz:
            if not hasattr(self, "_tz_name"):
                self._set_tz_name()
            set_col_tz(data, self._tz_name, inplace=True)
        # 3. Avg daily values
        daily = data.groupby(data.index.date).mean()
        p_cols = get_cols(data.columns, variable=r'P')
        if p_cols:
            daily[p_cols] = data[p_cols].groupby(data[p_cols].index.date).sum()
        if dropna:
            daily = daily.dropna(how='all')
        # Convert index to pandas dt (to match read_daily)
        daily.index = pd.to_datetime(daily.index)
        # self.daily = daily
        return daily
    
    def get_anc_data(self, variables=['PET','P', 'LAI']):
        # if sub_dirs is None:
        #     sub_dirs = [d.name for d in os.scandir(anc_dir) if d.is_dir()]
        # sub_dirs = [os.path.join(anc_dir, d) for d in sub_dirs]
        anc_dfs = [self._get_anc_data(v) for v in variables if not v is 'LAI'] 
        if 'LAI' in variables:
            lai = self.get_lai()
            # lai = lai.resample('D').asfreq()
            # lai['LAI'] = lai.LAI.ffill()
            anc_dfs.append(lai)
        
        anc_data = pd.concat(anc_dfs, axis=1)

        return anc_data

    def _get_anc_data(
        self, var, dt_col='time', pre_cols=['total60d','total_month'],
        drop_cols = ['network','station','latitude','longitude']
    ):
        # Get filename + column
        sub_dir, suff, col = anc_dict.get(var, ('','',''))
        file = os.path.join(anc_dir, sub_dir, f'{self.fid}_{suff}.csv')
        # Read file
        if os.path.exists(file):
            print(f"Found ancillary data: {file}.")
            df = pd.read_csv(file)
            df[dt_col] = pd.to_datetime(df[dt_col]) #, format='%Y-%m-%d')
            df.set_index(dt_col, inplace=True)
            df.rename(columns=var_dict, inplace=True)
            # return df
        else:
            print(f"File not found: {file}.")
            df = pd.DataFrame()
        if pre_cols:
            df.rename(columns={c : f'{col}_{c}' for c in pre_cols}, inplace=True)
        if drop_cols:
            df.drop(columns=drop_cols, inplace=True, errors='ignore')
        return df
    
    def get_lai(self, ):
        lai = self._get_anc_data('LAI', dt_col='Date')
        if not lai.empty:
            lai['LAI_raw'] = lai.LAI.copy()
            # Mask bad data
            lai['LAI'] = lai.LAI[lai.LAI < 248]
            # Fill missing values
            lai = lai.resample('D').asfreq()
            lai['LAI'] = lai.LAI.ffill()
        return lai

    

    # def mask_data(self):


    # def save(self, )

    def save(self, df, suff, out_dir):
        check_dir(out_dir)
        filename = f'{self.fid}_{suff}.csv'
        print(f"Saving {filename} to {out_dir}")
        df.to_csv(os.path.join(out_dir, filename))
        print(f"Saved {filename} to {out_dir}")
    



class Instrument():
    instrument_vars = Station.station_vars + ['instrument']
    def __init__(
            self, name : str, station : str, instrument_label=None, 
            sensors : list=None, keep_data : bool=False
        ):

        self.id = name
        self.station = station
        self.label = instrument_label if instrument_label is not None else 1

        # self.sensors = sensors
        # self.sensors = sensors if sensors else OrderedDict([])
        self.sensors = OrderedDict([
            (sensor.name, sensor) for sensor in sensors
        ])
        # self.variable = sensors[0].variable
        self.variables = list(set([s.variable for s in self.iter_sensors()]))
        self.depths = self.get_depths()
        
        self._meta = None
        self._var_info = None
        self._data = None
        self._sensor_cols = None
        self._keep_data = keep_data

    def __repr__(self):
        return f"Instrument '{self.id}' with Sensors: {[s.name for s in self.sensors.values()]}"

    # @property
    # def meta(self):
    #     if self._meta is None:
    #         self._meta = self.var_info[self.instrument_vars].iloc[0].to_dict()
    #         self._meta['n_sensors'] = 

    @property
    def var_info(self):
        if self._var_info is None:
            self._var_info = META[
                (META.station == self.station) & (META.instrument == self.id)
            ].sort_values(by='sensor')
        return self._var_info

    @property
    def data(self):
        if self._data is None:
            self._data = pd.concat(self.get_sensor_data(), axis=1)
        return self._data 
    
    @property
    def data(self):
        # if self._data is None:
        #     instrument_data = [inst.data for inst in self.iter_instruments()]
        #     self._data = pd.concat(instrument_data, axis=1)
        # return self._data
        return self.read_data()

    # @data.setter
    def read_data(self):
        if self._data is None:
            data = pd.concat(self.get_sensor_data(), axis=1)
            if self._keep_data:
                self._data = data
            else:
                return data
        else:
            return self._data
    
    @property
    def sensor_cols(self):
        if self._sensor_cols is None:
            self._sensor_cols = self.get_sensor_cols()
        return self._sensor_cols

    def get_sensor_cols(self, keys='sensor', **filter_kwargs):
        # sensor_cols = {
        #     sensor._col_name(i=self.label, d=d+1) : sensor.name for d,sensor in enumerate(
        #         self.iter_sensors(**filter_kwargs)
        #     )
        # }
        sensor_cols = {
            sensor.name : sensor._col_name(i=self.label, d=d+1) for d,sensor in enumerate(
                self.iter_sensors(**filter_kwargs)
            )
        }
        if not keys == 'sensor':
            sensor_cols = {v:k for k,v in sensor_cols.items()}
        return sensor_cols

    @property
    def n_sensors(self):
        return len(list(self.iter_sensors()))

    def get_sensor_data(self):
        sensor_data = [
            # sensor.rename_cols(i=self.label, d=d+1) for d,sensor in enumerate(self.iter_sensors())
            sensor.get_data(
                col_args={'i':self.label, 'd':d+1},
            ) for d,sensor in enumerate(self.iter_sensors())
        ]
        return sensor_data

    def iter_sensors(self, **filter_kwargs):
        """
        Iterates over all sensors in this station and yields those that
        comply with the passed filter settings (or all).

        Parameters
        ----------
        Keyword arguments are used to check all sensors at all stations,
        only stations that have at least one matching sensor are returned.

        For a description of possible filter kwargs, see
        :func:`ismn.components.Sensor.eval`

        Yields
        ------
        sensors : Sensor
            (Filtered) Sensors at the Station.
        """
        for sensor in self.sensors.values():
            if sensor.eval(**filter_kwargs):
                yield sensor
    
    def get_depth_ids(self):
        depths = sorted(list(set([(d[0],d[1]) for d in self.depths])))
        depth_map = {tup:i for i,tup in enumerate(depths)}
        return depth_map

    def get_depths(self, variable=None):
        return [s.depth for s in self.iter_sensors(variable=variable)]


class Sensor(components.Sensor):
    def __init__(
            self, station, instrument, variable, depth : components.Depth, name=None, filehandler=None, keep_loaded_data=None
        ):
        super().__init__(
            instrument=instrument, variable=variable, depth=depth, name=name, 
            filehandler=filehandler, keep_loaded_data=keep_loaded_data
        )
        # self.meta = META[(META.station == station) & (META.sensor == self.id)].squeeze().to_dict()
        self.station = station
        self.dz = 0.05 
        self.z = self.depth.extent
        # self.lat, self.lon = self.meta['latitude'], self.meta['longitude']
        self._meta = None

    @property
    def meta(self):
        if self._meta is None:
            # meta = self.metadata.to_pd().reset_index().pivot(
            #     index='variable', columns='key', values='data'
            # ).val.to_dict()
            # meta['depth_from'] = self.depth.start
            # meta['depth_to'] = self.depth.end
            # if self.metadata['clay_fraction'].depth:
            #     meta['tdepth_from'] = self.metadata['clay_fraction'].depth.start
            #     meta['tdepth_to'] = self.metadata['clay_fraction'].depth.end
            # else:
            #     meta['tdepth_from'],meta['tdepth_to'] = (np.nan,np.nan)
            # meta['name'] = self.id
            self._meta = META[
                (META.station == self.station) & (META.sensor == self.name)
            ].squeeze().to_dict()

        return self._meta

    # TODO: Add function to get soil texture class

    # def __gt__(self, other):
    #     if isinstance(other,Sensor) & (self.instrument == other.instrument):
    #         (self.depth[0] > other.depth[0]) & (self.depth[1] > other.depth[1])

    def eval(
        self,
        variable=None,
        depth=None,
        filter_meta_dict=None,
        check_only_sensor_depth_from=False,
        instrument=None,
    ):
        """
        Evaluate whether the sensor complies with the passed metadata
        requirements.

        Parameters
        ----------
        variable : str or list[str], optional (default: None)
            Check if the variable name matches, e.g. soil_moisture.
            One or multiple of :const:`ismn.const.VARIABLE_LUT`
        depth : Depth or list or tuple, optional (default: None)
            Check if the passed depth encloses the sensor depth.
            A list/tuple must contain 2 values where the first is the depth start
            and the second is the end. Start must be closer to 0 than end (or equal).
            A negative depth range is above the surface.
        filter_meta_dict : dict, optional (default: None)
            Additional metadata keys and values for which the file list is filtered
            e.g. {'lc_2010': [10, 130]} or
                 {'climate_KG': 'Dwa', 'lc_2010': [10, 130] }
            to filter for a multiple landcover classes and a climate class.
        check_only_sensor_depth_from : bool, optional (default: False)
            Ignores the sensors depth_to value and only checks if depth_from of
            the sensor is in the passed depth (e.g. for cosmic ray probes).

        Returns
        -------
        flag : bool
            Indicates whether metadata for this Sensor matches with the passed
            requirements.
        """
        flag = super().eval(
            variable=variable, depth=depth, filter_meta_dict=filter_meta_dict, 
            check_only_sensor_depth_from=check_only_sensor_depth_from
        )
        if instrument is not None:
            if self.instrument != instrument:
                flag = False
        return flag
    
    # @property
    # def col_name(self):
    #     return self._col_name
    
    # @col_name.setter
    def _col_name(self, i, d, r=None):
        new_col = var_dict.get(self.variable) + f'_{i}_{d}'
        # new_col += f'_{i}_{d}'
        if r:
            new_col += f'_{r}'
        self.col_name = new_col
        return self.col_name
    
    def get_data(
            self, col_args={'i' : '', 'd' : ''}, 
            exclude='_flag', suff='_masked', cond='G'
        ):
        data = self.rename_cols(**col_args)
        self.mask_data(data, exclude=exclude, suffix=suff, cond=cond)
        return data

    def mask_data(self, data, exclude='_flag', suffix='_masked', cond='G'):
        data[self.col_name + suffix] = data[self.col_name][
            data[self.col_name + exclude] == cond
        ]
        


    def rename_cols(self, i, d, r=None):
        # variable = list(df.columns)[0]
        # new = get_col_name(variable, i, d, r)
        # cols = [col.replace(variable, new) for col in df.columns]
        cols = [col.replace(self.variable, self._col_name(i, d, r)) for col in self.data.columns]
        return self.data.rename(columns=dict(zip(self.data.columns, cols)))

    def _set_tz(self, tz_name=None, lat=None, lon=None):
        if not tz_name:
            if not lat or lon:
                lat,lon = self.meta['latitude'], self.meta['longitude']
            tz_name = tzf.TimezoneFinder().timezone_at(
                lng=lon, lat=lat
            )
        utc = self.data.index.tz_localize('UTC')
        local = utc.tz_convert(tz_name)
        return local