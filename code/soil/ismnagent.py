#-------------------------------------------------------------------------------
# IMPORTS
#-------------------------------------------------------------------------------
# from SMAPgrid import SMAPgrid
import warnings
from datetime import datetime
import os
import getpass
import numpy as np
import pandas as pd
import logging

# import station
# from station import Station, META, ismn_data
from . import station
from .station import META, Station #, ismn_data

from drydowns import ISMNSensorData, DrydownModelHandler
# from drydowns import ISMNSoilData, DrydownModel
# from drydowns.ismndata import ISMNSoilData
# from .model import DrydownModel
# from .towerseparator import TowerEventSeparator

from drydowns.mylogger import getLogger

# Create a logger
log = getLogger(__name__)


#-------------------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------------------

def create_output_dir(parent_dir):
    username = getpass.getuser()
    formatted_now = datetime.now().strftime("%Y-%m-%d")
    output_dir = rf"{parent_dir}/{username}_{formatted_now}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        log.info(f"Directory '{output_dir}' created.")
    return output_dir



# def get_cols(
#     tower : fluxtower.FluxTower, 
#     var_cols : list = ['SWC', 'P', 'TA', 'VPD', 'LE', 'ET', 'e_a']
# ):
#     col_list = []
#     col_dict = {}
#     for var in var_cols:
#         col_list += tower.get_var_cols(variable=var, exclude='ERA')
#         cols = tower.get_var_cols(variable=var, exclude='ERA')
#         col_dict[var] = {
#             'var_cols' : sorted([col for col in cols if 'QC' not in col]), 
#             'qc_cols' : sorted([col for col in cols if 'QC' in col])
#         }

#     return col_list, col_dict

# def get_grps(tower, cols):
#     vi = tower._get_var_info(as_dict=False)
#     grps = vi[vi.VARNAME.isin(cols)].GROUP_ID.unique()
#     return grps

#-------------------------------------------------------------------------------
# VARIABLES
#-------------------------------------------------------------------------------

# meta = fluxtower.flx_tower.META

#-------------------------------------------------------------------------------
# TOWERAGENT CLASS
#-------------------------------------------------------------------------------

class ISMNAgent:
    def __init__(self, cfg=None, logger=None, station_ids : list = None, output_dir=None):
        # config
        self.cfg = cfg
        # logger
        self.logger = logger

        # data_ids
        # If a list of tower ids is provided, use that.
        if station_ids:
            self.data_ids = station_ids
        # Otherwise, use all the towers with soil moisture data.
        else:
            self.data_ids = pd.Series(
                list(zip(
                    META[META.variable == 'soil_moisture'].network,
                    META[META.variable == 'soil_moisture'].station
                ))
            ).unique()
        # filenames
        # self._filenames = sorted(os.listdir(cfg["PATHS"]["data_dir"]))
        # self._filenames = sorted(os.listdir(cfg["data_dir"]))

        # verbose
        # self.verbose = cfg["MODEL"]["verbose"].lower() in ["true", "yes", "1"]
        self.verbose = cfg.get("verbose").lower() in ["true", "yes", "1"]

        # output_dir
        if output_dir:
            self._output_dir = output_dir
        else:
            # self._output_dir = create_output_dir(parent_dir=cfg["PATHS"]["output_dir"])
            # self._output_dir = cfg.get('PATHS', 'output_dir')
            self._output_dir = cfg.get('output_dir')

    def initialize(self):
        None

    def run(self, did, return_data=False):
        """
        Run the analysis for one tower (which may contain multiple soil moisture
        columns)
        """
        try: 
            if self.verbose:
                log.info(f"Currently processing tower {did}")
            
            # 1. Initialize tower
            log.info(f"Initalizing Station {did}")
            station = self._init_station(did)

            # 2. Get list of soil moisture sensors
            sm_sensors = [i for i in station.iter_sensors(variable='soil_moisture')]
            # cols, col_dict = get_cols(tower)
            # sm_cols = col_dict['SWC']['var_cols']
            # grps = get_grps(tower, sm_cols)

            log.info(f"Found {len(sm_sensors)} soil moisture columns for {did} : {sm_sensors}")

            # 3. Run the analysis for each soil moisture column
            data = []
            results = []
            # for col in sm_cols:
            #     output = self.run_sensor(tower, col, return_data=return_data)
            # for grp in grps:
            for sensor in station.iter_sensors(variable='soil_moisture'):
                output = self.run_sensor(station, sensor, return_data=return_data)
                if return_data:
                    data.append(output[0])
                    results.append(output[1])
                else:
                    results.append(output)
            
            # If all results are NOne, return
            if all([r is None for r in results]):
                log.warning(f"No drydown events detected for {did}")
                if return_data:
                    return data, None
                else:
                    return None
             

            log.info(
                f"Finished with tower {did}: {len(results)}/{len(sm_sensors)} sensors analyzed"
            )

            results = pd.concat(results, ignore_index=True)
            # results['IGBP'] = tower.metadata.get('IGBP', None)
            # results['LAT'] = tower.coords[0]
            # results['LON'] = tower.coords[1]
            # results['MAT'] = tower.metadata.get('MAT', np.nan)
            # results['MAP'] = tower.metadata.get('MAP', np.nan)


            if return_data:
                return data, results
            else:
                return results

        except Exception as e:
            print(f"Error in thread: {did}")
            print(f"Error message: {str(e)}")


    # def run_sensor(self, tower, swc_col, return_data=False):
    def run_sensor(self, station, sensor, return_data=False):
        """ Run the analysis for one soil moisture sensor"""
        try:
            
            # 1. Initialize sensor data object
            log.info(f"Initializing sensor {sensor.name}")
            data = ISMNSensorData(self.cfg, station=station, sensor_name=sensor.name)

            # 2. Separate events
            log.info(f"Separating events for sensor {sensor.name}")
            data.separate_events()

            if not data.events:
                log.warning(f"No event drydown was detected for {sensor.name}")
                if return_data:
                    return data, None
                else:
                    return None
            else:
                log.info(f"Found {len(data.events)} events for sensor {sensor.name}")

            # 3. Fit drydown models
            log.info(f"Fitting drydown models for sensor {sensor.name}")
            # model = DrydownModel(self.cfg, data, data.events)
            # model.fit_models(output_dir=self._output_dir)
            # # 4. Return results
            # results = model.return_result_df()
            handler = DrydownModelHandler(self.cfg, data, data.events)
            handler.fit_events()
            results = handler.get_results()

            log.info(
                f"Drydown model analysis completed for {sensor.name}: {len(results)}/{len(data.events)} events fitted"
            )

            results['LC'] = sensor.meta.get('lc_2010')
            results['LAT'] = sensor.meta.get('latitude')
            results['LON'] = sensor.meta.get('longitude')
            results['depth'] = [
                (sensor.meta.get('depth_from'),sensor.meta.get('depth_to'))
            ]*len(results)
            results['texture'] = sensor.meta.get('soil_texture')

            if return_data:
                return data, results
            else:
                return results

        except Exception as e:
            print(f"Error in thread: {sensor.name}")
            print(f"Error message: {str(e)}")



    def _init_station(self, did):
        """Initialize a tower object from the Site ID (str)"""
        # # filename = next((t for t in self._filenames if tid in t),None)
        
        # # tower = fluxtower.FluxNetTower(os.path.join(self.cfg["PATHS"]["data_dir"], filename))
        # tower = fluxtower.FluxNetTower(os.path.join(self.cfg.get("data_dir"), filename))
        # tower.add_vp_cols()
        # tower.add_et_cols()
        station = Station(self.cfg, station=did[1], network=did[0], keep_data=True)
        # station.daily
        return station

    def finalize(self, results):
        """Finalize the analysis from all the pixels

        Args:
            results (list): concatenated results returned from serial/multi-threadding analysis
        """
        try:
            df = pd.concat(results)
        except:
            df = results

        # date = datetime.now().strftime("%Y-%m-%d")
        date = datetime.now().strftime('%d%b').lower()
        out_bn = self.cfg.get(
            'output_fid',
            'ismn_results'
        )
        fid = f"{out_bn}_{date}"

        self.save(df, fid=fid)


    def save(self, df, fid='ismn_results'):
        """Save the results to a csv file"""
        filename = f"{fid}.csv"
        log.info(f"Saving {filename} to {self._output_dir}")
        df.to_csv(os.path.join(self._output_dir, filename), index=False)
        log.info(f"Saved {filename} to {self._output_dir}")

