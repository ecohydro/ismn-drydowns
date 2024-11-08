
#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Thu 04 Apr 24 18:34:12'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           drydown.py
Compatibility:  Python 3.12.0

Description:    Class + functions for handling and plotting fitted drydowns.

URL:            https://

Requires:       list of libraries required

Dev ToDo:       None

AUTHOR:         Bryn Morgan
ORGANIZATION:   University of California, Santa Barbara
Contact:        bryn.morgan@geog.ucsb.edu
Copyright:      (c) Bryn Morgan 2024


"""


# IMPORTS
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import var_dict

from soil.station import Station
from drydowns import ISMNSensorData, figs
from drydowns.config import Config


config = Config('soil/config.ini').config

# FUNCTIONS



def calc_theta_exp(t, theta_0, theta_w, tau, theta_star=None):
    """
    Calculate the drydown curve for soil moisture over time using non-linear plant stress model.

    Parameters:
        t (int): Timestep, in day.
        delta_theta (float): Shift/increment in soil moisture after precipitation, in m3/m3.
        theta_w (float, optional): Wilting point soil moisture content, equal to s_star * porosity, in m3/m3. Default is 0.0.
        tau (float): decay rate, in 1/day.

    """
    # delta_theta = theta_0 - theta_w
    # t_star = (theta_0 - theta_w) / tau
    # theta = np.where()
    if theta_star is not None:
        t_star = - np.log((theta_star - theta_w) / (theta_0 - theta_w)) * tau
        k = (theta_star - theta_w) / tau
        theta = np.where(t < t_star, -k * t + theta_0, (theta_0 - theta_w) * np.exp(-t / tau) + theta_w)
    else:
        theta = (theta_0 - theta_w) * np.exp(-t / tau) + theta_w
    return theta

def calc_theta_lin(t, theta_0, k, theta_star=1.0, theta_w=0.0, z=0.05, buffer=0.0001):
    """
    Calculate the drydown curve for soil moisture over time using linear model.

    Parameters:
        t (int): Timestep, in day.
        k (float): Product of soil thickness (z) and maximum rate of change in normalized soil moisture (k), equivalent to maximum ET rate (ETmax), in m3/m3/day.
        delta_theta (float): Shift/increment in soil moisture after precipitation, in m3/m3. It is equal to theta_0 - theta_w.
        theta_star (float, optional): Critical soil moisture content, equal to s_star * porosity, in m3/m3. Default is 1.0.
        theta_w (float, optional): Wilting point soil moisture content, equal to s_star * porosity, in m3/m3. Default is 0.0.

    Returns:
        float: Rate of change in soil moisture (dtheta/dt) for the given timestep, in m3/m3/day.
    """
    t_star = np.round((theta_0 - theta_star)/k,6) # Time it takes from theta_0 to theta_star 
    theta = (theta_0 - theta_w) * np.exp((-k/(theta_star - theta_w)) * t) + theta_w
    theta = np.where(theta > theta_star, -k * t + theta_0, theta)
    # theta = np.where(theta < theta_w, (-0.05*z*1000)*t, theta)
    theta = np.where(theta < theta_w+buffer, np.nan, theta)

    return theta


def calc_theta_q(t, theta_0, k, q, theta_star=1.0, theta_w=0.0, z= 0.05, buffer=0.0001):
    """
    Calculate the drydown curve for soil moisture over time using non-linear plant stress model.

    Parameters:
        t (int): Timestep, in day.
        k (float): Product of soil thickness (z) and maximum rate of change in normalized soil moisture (k), equivalent to maximum ET rate (ETmax), in m3/m3/day.
        q (float): Degree of non-linearity in the soil moisture response.
        delta_theta (float): Shift/increment in soil moisture after precipitation, in m3/m3. It is equal to theta_0 - theta_w.
        theta_star (float, optional): Critical soil moisture content, equal to s_star * porosity, in m3/m3. Default is 1.0.
        theta_w (float, optional): Wilting point soil moisture content, equal to s_star * porosity, in m3/m3. Default is 0.0.

    Returns:
        float: Rate of change in soil moisture (dtheta/dt) for the given timestep, in m3/m3/day.
    """
    t_star = np.round((theta_0 - theta_star)/k,6) # Time it takes from theta_0 to theta_star 
    # t_w = (theta_0 - theta_w)/k # Time it takes from theta_0 to theta_w
    # delta_theta = theta_0 - theta_w

    b = (theta_0 - theta_w) ** (1 - q)

    a = (1 - q) / ((theta_star - theta_w) ** q)
    theta = (-k * a * t + b) ** (1 / (1 - q)) + theta_w
    theta = np.where(theta > theta_star, -k * t + theta_0, theta)
    # theta = np.where(t < t_star, -k * t + theta_0, (-k * a * t + b) ** (1 / (1 - q)) + theta_w)
    # theta = np.where(theta < theta_w+buffer, (-0.05*z*1000)*t, theta)
    # theta = np.where(theta < theta_w+buffer, np.nan, theta)
    # theta = np.where(t <= t_star, -k * t + theta_0, (-k * a * t + b) ** (1 / (1 - q)) + theta_w)

    return theta


def calc_dtheta_dt(theta, k, q, theta_w=0.0, theta_star=1.0, buffer=0.0001):
    """
    Calculate the rate of change in soil moisture over time using non-linear plant stress model.

    Parameters:
        t (int): Timestep, in day.
        k (float): Product of soil thickness (z) and maximum rate of change in normalized soil moisture (k), equivalent to maximum ET rate (ETmax), in m3/m3/day.
        q (float): Degree of non-linearity in the soil moisture response.
        theta (float): Soil moisture content at the given timestep, in m3/m3.
        theta_star (float, optional): Critical soil moisture content, equal to s_star * porosity, in m3/m3. Default is 1.0.
        theta_w (float, optional): Wilting point soil moisture content, equal to s_star * porosity, in m3/m3. Default is 0.0.

    Returns:
        float: Rate of change in soil moisture (dtheta/dt) for the given timestep, in m3/m3/day.
    """
    dtheta_dt = np.where(
        theta > theta_star, 
        k, 
        k * ( ( (theta - theta_w) / (theta_star - theta_w) )**q )
    )
    theta_h = theta_w * 0.5
    e_w = dtheta_dt[theta > theta_w][-1]
    
    dtheta_dt = np.where(
        theta < theta_w, e_w * ((theta - theta_h)/(theta_w-theta_h)), dtheta_dt
    )
    # dtheta_dt = k * ( ( (theta - theta_w) / (theta_star - theta_w) )**q )
    return dtheta_dt





# def plot(
#         ax, df, x, y, label=None, xlim=None, ylim=None,
#         kwargs={'color' : 'k', 'alpha' : 0.8}, legend=True
#     ):
#     ax.plot(
#         # self.df.t, self.df['theta'+suff], '.-', label=labels.get(suff), 
#         df[x], df[y], '.-', label=label,
#         **kwargs #color='k', alpha=0.8
#     )
#     ax.set_xlabel(
#         # var_dict.get(re.match(r'^([a-zA-Z]+(?:_mm)?)(?:_obs)?$', x).group(1))['label']
#         var_dict.get(x.rsplit('_', 1)[0])['label']
#     )
#     ax.set_ylabel(
#         # var_dict.get(re.match(r'^([a-zA-Z]+(?:_mm)?)(?:_obs)?$', y).group(1))['label']
#         var_dict.get(y.rsplit('_', 1)[0])['label']
#     )
#     if xlim:
#         ax.set_xlim(xlim)
#     if ylim:
#         ax.set_ylim(ylim)
#     if legend:
#         ax.legend(loc='upper right', bbox_to_anchor=(0.98, 1))

#     return ax



# DRYDOWN CLASS

class Drydown():
    labels = {
        '_obs' : 'Observed',
        '_exp' : 'Linear',
        '_q' : r"$q$",
    }


    def __init__(
        self, 
        t, theta_obs, pet,
        z, theta_w, theta_star, 
        et=None,
        exp_params=None, q_params=None,
        site_info=None,
        event_id=None,
        start_date=None, end_date=None, 
    ):
        self.id = event_id
        self.start_date = start_date
        self.end_date = end_date

        self.t = t
        self.z = z
        
        # self.k = k
        # self.ET_max = k * (z * 1000)

        self._theta_obs = theta_obs
        self.theta_w = theta_w
        self.theta_star = theta_star

        self.exp_params = self._unpack_params(exp_params, prefix='exp_')
        self.q_params = self._unpack_params(q_params, prefix='q_')
        self.q = self.q_params.get('q', 1.0)
        
        self.pet = pet
        self.et = et
        self.set_et_max()
        self.site_info = site_info

        self.df = self.create_df()

    def create_df(self):
        data = {
            't' : self.t,
            'theta_obs' : self._theta_obs,
            'dtheta_obs' : self.calc_dtheta_dt(self._theta_obs, method='obs'),
            'dtheta_exp' : self.calc_dtheta_dt(self._theta_obs, method='exp'),
            'dtheta_q' : self.calc_dtheta_dt(self._theta_obs, method='q'),
        }
        df = pd.DataFrame(data)
        # Get/calculate theta_exp
        theta_exp = self.exp_params.get(
            'theta_opt', 
            self.calc_theta_exp(
                self.t, #self.exp_params['theta_0'], self.theta_w, self.exp_params['tau']
            )
        )
        # Get/calculate theta_q
        theta_q = self.q_params.get(
            'theta_opt', 
            self.calc_theta_q(
                self.t, #self.q_params['k'], self.q_params['q'], self.q_params['theta_0'], 
                # self.theta_star, self.theta_w
            )
        )
        df.insert(2, 'theta_exp', theta_exp)
        df.insert(3, 'theta_q', theta_q)

        for suff in ['_obs', '_exp', '_q']:
            df['dtheta_mm'+suff] = df['dtheta' + suff] * self.z * 1000
        
        if self.et is not None:
            df['et'] = self.et

        return df

    def _unpack_params(self, params, prefix='exp_'):
        try :
            params_ed = {}
            for key in ['theta_0', 'theta_opt', 'k', 'ET_max', 'r_squared']:
                params_ed.update({key : params.get(prefix + key)})
            # params_ed =  {
            #     'theta_0' : params[prefix + 'theta_0'],
            #     'theta_opt' : params[prefix + 'theta_opt'],
            #     'k' : params[prefix + 'k'],
            #     'ET_max' : params[prefix + 'ET_max'],
            #     'r2' : params[prefix + 'r_squared'],
            # }
            if prefix == 'exp_':
                params_ed.update({
                    'tau' : params[prefix + 'tau'],
                    'theta_w' : params[prefix + 'theta_w'],
                })
            elif prefix == 'q_':
                params_ed.update({
                    'q' : params[prefix + 'q'],
                    # 'k' : params[prefix + 'k']
                })
        except:
            params_ed = {}

        return params_ed
    
    def set_et_max(self):
        self.ET_max = self.q_params['k'] * (self.z * 1000)


    def calc_theta_q(self, t):
        theta_0 = self.q_params['theta_0']
        k = self.q_params['k']

        return calc_theta_q(
            t=t, theta_0=theta_0, k=k, q=self.q, 
            theta_star=self.theta_star, theta_w=self.theta_w
        )



    def calc_theta_exp(self, t):
        theta_0 = self.exp_params['theta_0']
        theta_w = self.theta_w
        tau = self.exp_params['tau']

        return calc_theta_exp(t=t, theta_0=theta_0, theta_w=theta_w, tau=tau)
    # def calc_theta_q(self, t):
    #     theta_0 = self.q_params['theta_0']
    #     k = self.q_params['k']

    #     delta_theta = theta_0 - self.theta_w
    #     b = delta_theta ** (1 - self.q)
    #     a = (1 - self.q) / ((self.theta_star - self.theta_w) ** self.q)
    #     return (-k * a * t + b) ** (1 / (1 - self.q)) + self.theta_w
    

    # def calc_theta_exp(self, t, theta_0, theta_w, tau):
    #     delta_theta = theta_0 - theta_w
    #     return delta_theta * np.exp(-t / tau) + theta_w


    def calc_dtheta_dt(self, theta, method='q'):

        if method in ['q', 'nonlinear']:
            k = self.q_params['k']
            dtheta_dt = self._dtheta_dt_analytical(
                self.t, k, self.q, theta=theta, 
                theta_w=self.theta_w, theta_star=self.theta_star
            )
        elif method in ['exp','exponential', 'linear']:
            k = self.exp_params['k']
            dtheta_dt = self._dtheta_dt_analytical(
                self.t, k, q=1,
                theta=theta, theta_w=self.theta_w, theta_star=self.theta_star
            )
        else:
            dtheta_dt = self._dtheta_dt_obs(theta)
        
        return dtheta_dt

    def _dtheta_dt_analytical(self, t, k, q, theta, theta_w=0.0, theta_star=1.0):
        dtheta_dt = k * ( ( (theta - theta_w) / (theta_star - theta_w) )**q )
        return dtheta_dt

    def _dtheta_dt_obs(self, theta):
        return pd.Series(theta).diff().to_numpy()*-1
    
    def plot_theta(
            self, ax, suff=['_obs', '_exp', '_q'], xlim=None, ylim=None, 
            legend=True, kwargs={'s' : 12, 'linewidth' : 0, 'alpha' : 0.8}, buffer=5,
        ):
        for s in suff:
            figs.plot(
                ax=ax, df=self.df, x='t', y='theta'+s, 
                label=self.labels.get(s), xlim=xlim, ylim=ylim, legend=legend,
                kwargs=kwargs
            )

        return ax
    
    def plot_dtheta(
            self, ax, suff=['_obs', '_exp', '_q'], units='mm', xlim=None, ylim=None,
            legend=True, kwargs={'s' : 12, 'linewidth' : 0, 'alpha' : 0.8}, plot_et=True,
        ):
        if plot_et and 'et' in self.df.columns:
            figs.plot(
                ax=ax, df=self.df, x='theta_obs', y='et', label='ET', xlim=xlim, 
                ylim=ylim, legend=legend, kwargs={'color' : 'r', 'alpha' : 0.8}
            )
        for s in suff:
            if units == 'mm':
                ds = '_mm' + s
            else:
                ds = s
            figs.plot(
                ax=ax, df=self.df, x='theta'+s, y='dtheta'+ds, 
                label=self.labels.get(s), xlim=xlim, ylim=ylim, legend=legend,
                kwargs=kwargs
            )
        return ax

    def plot_drydown(
            self, axs, suff=['_obs', '_exp', '_q'], units='mm', xlim=None, ylim=None,
            kwargs={'s' : 12, 'linewidth' : 0, 'alpha' : 0.8}, plot_et=True,
        ):
        self.plot_theta(
            axs[0], suff=suff, xlim=xlim, ylim=ylim, legend=False, kwargs=kwargs, 
        )
        self.plot_dtheta(
            axs[1], suff=suff, units=units, xlim=xlim, ylim=ylim, legend=True, kwargs=kwargs, plot_et=plot_et
        )
        axs[1].legend(loc='upper left', bbox_to_anchor=(0.01, 1))
        # lat long with 2 digits
        if self.site_info:
            title = f"{self.id[0]}, {self.site_info['IGBP']}, z = {self.z} m "
            if self.site_info.get('LAT') and self.site_info.get('LON'):
                title += f"({self.site_info['LAT']:.2f}, {self.site_info['LON']:.2f})"
            plt.gcf().suptitle(title)
        return axs


    def get_sensor_data(self):
        station = Station(
            config=config['ISMN'],
            station_id=self.id[0],
            network=self.site_info['network'],
            keep_data=True,
        )
        self._sens_data = ISMNSensorData(config=config, station=station, sensor_name=self.id[1])
        return self._sens_data


    def get_event_data(self, buffer=5):
        # self._event_data = self._sens_data.get_event_data(
        #     start_date=self.start_date, end_date=self.end_date, buffer=buffer
        # )



        return self._event_data

