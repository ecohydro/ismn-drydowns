#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Thu 04 Apr 24 18:38:08'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           utils_figs.py
Compatibility:  Python 3.12.0

Description:    Utility functions for creating figures.

URL:            https://

Requires:       list of libraries required

Dev ToDo:       None

AUTHOR:         Bryn Morgan
ORGANIZATION:   University of California, Santa Barbara
Contact:        bryn.morgan@geog.ucsb.edu
Copyright:      (c) Bryn Morgan 2024


"""


# IMPORTS
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import numpy as np
import pandas as pd
import getpass
import re
import seaborn as sns
import cartopy.crs as ccrs

# SETTINGS

mpl.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
font_files = mpl.font_manager.findSystemFonts(
    fontpaths=[os.path.join('home', getpass.getuser(), 'Fonts')]
)

for font_file in font_files:
    mpl.font_manager.fontManager.addfont(font_file)



mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['font.sans-serif'] = 'Myriad Pro'
# mpl.rcParams['font.sans-serif'] = 'Source Sans Pro'
mpl.rcParams['font.size'] = 12.0
mpl.rcParams['axes.titlesize'] = 12.0

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

# FUNCTIONS

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]


def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp


def normalize(x, xmin=None, xmax=None):
    if xmin is None:
        xmin = np.min(x)
    if xmax is None:
        xmax = np.max(x)
    return (x - xmin)/(xmax - xmin)

def rescale(x, xmin=None, xmax=None, new_min=0.25, new_max=2.2):
    return normalize(x, xmin, xmax)*(new_max - new_min) + new_min


def get_labs(df, col):
    key = f'{col}_bin' if 'bin' not in col else col
    labs = bin_dict.get(re.findall(r'.*?bin', key)[0])[f'{re.findall(r'bin.*$', key)[0]}'][:-1]
    # if not isinstance(labs, np.ndarray):
    #     labs = np.array(labs)
    # Format labels to be 2 digits
    # labs = [f'{i:.2f}'.rstrip('0').rstrip('.') for i in labs]
    return labs

def get_ticks(df, x, xmax):
    xlabs = get_labs(df, x)
    if xmax > len(xlabs):
        if not isinstance(xlabs, np.ndarray):
            xlabs = np.array(xlabs)
        xticks = rescale(xlabs, df[x].min(), df[x].max(), new_min=0, new_max=xmax)
    else:
        xticks = np.arange(len(xlabs))
    xlabs = [f'{i:.2f}'.rstrip('0').rstrip('.') for i in xlabs]
    return xticks, xlabs



def get_ax_ticks(df, ax, x=None, y=None):
    # if labs is None:
    xticks = get_ticks(df, x, np.max(ax.get_xlim())) if x else (None, None)
    yticks = get_ticks(df, y, np.max(ax.get_ylim())) if y else (None, None)
    
    return xticks, yticks


# VARIABLES

map_dict = {
    'Mercator' : {
        'crs' : ccrs.Mercator(),
        # 'extent' : [-20037508.34, 20037508.34, -20037508.34, 20037508.34],
        'extent' : [-179.4, 179.4, -56, 72],
        # 'projection' : 'mercator',
        # 'transform' : ccrs.PlateCarree(),
    },
    'PlateCarree' : {
        'crs' : ccrs.PlateCarree(),
        'extent' : [-180, 180, -56, 84],
        # 'projection' : 'cyl',
        # 'transform' : None,
    },
    'Robinson' : {
        'crs' : ccrs.Robinson(),
        'extent' : [-180, 180, -56, 84],
        'projection' : 'robin',
        # 'transform' : ccrs.PlateCarree(),
    },
}

pt_size_dict = {
    'q_med' : {
        'size_norm' : (0., 4.),
        's_labs' : np.array([0, 0.5, 1, 1.5,  2, 4, ]),
        'label' : r"$q$"
    },
    'n_events' : {
        'size_norm' : (1, 150),
        's_labs' : np.array([1, 10, 50, 150, 250]),
    },
    'LAI_max' : {
        'size_norm' : (1., 7.),
        's_labs' : np.array([1., 3., 5., 7.]),
        'label' : r'$\mathrm{LAI}_{\mathrm{max}}$''\n'r'(m$^2$ m$^{-2}$)'
    }
}


def get_sizes(var, df, size_norm, s_labs, leg_sizes=(4,50), **kwargs):
    s_min, s_max = leg_sizes
    n_min, n_max = size_norm
    # Sizes for plot
    s_norm = (df[var] - n_min) / (n_max - n_min)
    sizes = s_min + (s_max - s_min) * s_norm
    # Resize points for legend
    s_leg = s_min + (s_max - s_min) * ((s_labs - n_min) / (n_max - n_min))
    s_leg = np.sqrt(s_leg) * 2
    
    return sizes, s_leg


def add_size_legend(var, s_labs, s_leg, ax, **kwargs):
    # if not leg_ax:
    #     try:
    #         # leg_ax = fig.add_axes(ax.get_position(), frameon=False)
    #     except Exception as e:
    #         print(f'Must pass one of leg_ax or ax: {e}')
    #         return
    leg = ax.legend(
        handles = [
            mpl.lines.Line2D(
                [], [], linestyle='None', markeredgewidth=0, markerfacecolor='k', 
                marker='.', markersize=si
            ) for si in s_leg
        ],
        labels=[s for s in s_labs],
        **kwargs,
        title = pt_size_dict[var].get('label', '')
        # multialignment='center',
    )
    leg.get_title().set_multialignment('center')
    leg = ax.get_legend()
    return leg

# Get projection


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
    'AI_60d' : {
        'symbol' : r"AI$_{60}$",
        # 'label' : r"Aridity Index (60 days)",
        'label' : r"AI (60 days)",
        'units' : r"(P$_{60}$/PET$_{60}$)",
        'bins' : [0, 0.05, 0.2, 0.5, 0.65, 0.8, 1000.], 
        'bin_labels' : ['0-0.05','0.05-0.2','0.2-0.5','0.5-0.65','0.65-0.8', '>0.8']
    },
    'AI' : {
        'symbol' : r"AI",
        # 'label' : r"Aridity Index (P/PET)",
        'label' : r"AI (P/PET)",
        'units' : r"(P/PET)",
        'bins' : [0 ,0.05, 0.2, 0.5, 0.65, 0.8, 1., 1000.], 
        'bin_labels' : ['0-0.05','0.05-0.2','0.2-0.5','0.5-0.65','0.65-0.8', '0.8-1.0', '>1.0']
    },
    'LAI' : {
        'symbol' : r"LAI",
        'label' : r"LAI (m$^2$ m$^{-2}$)",
        'units' : r"(m$^2$ m$^{-2}$)",
        'bins' : [0, 0.5, 1, 2, 3, 5, 10],
        'bin_labels' : ['0-0.5','0.5-1','1-2','2-3','3-5','>5']
    },
    'MAP' : {
        'symbol' : r"MAP",
        'label' : r"Mean annual precipitation (mm)",
        'units' : r"(mm)",
        'bins' : [0,300,600,900,1200,1500,2000, 6000],
        'bin_labels' : ['0-300','300-600','600-900','900-1200','1200-1500','1500-2000','>2000']
    },
    'perc_woodyveg' : {
        'symbol' : r"Woody cover",
        'label' : r"Woody cover (%)",
        'units' : r"%",
        'bins' : [0, 10, 20, 30, 40, 50, 100],
        'bin_labels' : ['0-10%','10-20%','20-30%','30-40%','40-50%','>50%'],
    },
    # 'AI_bin' : {
    #     'symbol' : r"AI",
    #     'label' : r"Aridity Index",
    #     'units' : r"(MAP/MAE)",
    #     'order' : 
    # },
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

# #7FB800

igbp_dict = {
    'BAR' : {
        'label' : 'Barren',
        'color' : '#8B4513'
    },
    'NONVEG' : {
        'label' : 'Bare soil',
        'color' : '#8B4513'
    },
    'CRO' : {
        'label' : 'Croplands',
        'color' : '#F7C906'
    },
    'GRA' : {
        'label' : 'Grasslands',
        # 'color' : '#13BFB2'
        'color' : '#92BA31'
    },
    'SAV' : {
        'label' : 'Savannas',
        # 'color' : '#92BA31'
        # 'color' : '#13BFB2'
        'color' : '#15BFB2'
    },
    'WSA' : {
        'label' : 'Woody savannas',
        'color' : '#4C6903'
    },
    'MF' : {
        'label' : 'Mixed forest',
        'color' : '#00796b'
    },
    'FOR' : {
        'label' : 'Forest',
        'color' : '#388e3c'
    },
    'BLF' : {
        'label' : 'Broadleaf forests',
        'color' : '#388e3c'
        # 'color' : '#578C29',
        # 'color' : '#4F993A',
    },
    'DBF' : {
        'label' : 'Deciduous broadleaf forests',
        'color' : '#388e3c'
    },
    'EBF' : {
        'label' : 'Evergreen broadleaf forests',
        'color' : '#00796b'#'#388e3c'
    },
    'DNF' : {
        'label' : 'Deciduous Needleleaf forests',
        'color' : '#1b5e20'
    },
    'ENF' : {
        'label' : 'Evergreen needleleaf forests',
        'color' : '#1b5e20'
        # 'color' : '#0e6c46',
    },
    'EVF' : {
        'label' : 'Evergreen forests',
        'color' : '#1b5e20'
    },
    'OSH' : {
        'label' : 'Open shrublands',
        'color' : '#C99728'
    },
    'CSH' : {
        'label' : 'Closed shrublands',
        'color' : '#C99728'
    },
    'SHR' : {
        'label' : 'Shrublands',
        'color' : '#C99728'
    },
    'WET' : {
        'label' : 'Wetlands',
        'color' : '#808080'
    },
    'MOS' : {
        'label' : 'Mosaic',
        'color' : '#808080'
    },
    'HERB' : {
        'label' : 'Grasslands',
        # 'color' : '#13BFB2'
        'color' : '#92BA31'
    },
    'WOOD' : {
        'label' : 'Woodlands + forests',
        # 'color' : '#00796b'
        'color' : '#01796B'
    }
}

igbp_cats = {
    'og' : ['GRA', 'CSH', 'OSH', 'SAV', 'WSA', 'DBF', 'EBF', 'ENF'],
    '7' : ['GRA', 'SHR', 'SAV', 'WSA', 'DBF', 'EBF', 'ENF'],
    '6s' : ['GRA', 'SHR', 'SAV', 'WSA', 'BLF', 'ENF'],
    '6b' : ['GRA', 'SHR', 'SAV', 'DBF', 'EBF', 'ENF'],
    '5' : ['GRA', 'SHR', 'SAV', 'BLF', 'ENF'],
    '2' : ['HERB', 'WOOD'],
    '3' : ['HERB', 'SAV', 'WOOD'],
    'LAI' : ['HERB', 'SAV', 'WOOD', 'NONVEG'],
}

# c_list = ['#1B3A39', '#5BB8B0', '#82B960', '#B9C354', '#EDC01D']
# c_listr = ['#EDC01D', '#82B960', '#5BB8B0', '#1B3A39']
c_list1 = ['#004446','#026f6f', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
cmap = get_continuous_cmap(c_list1)



bin_dict = {
    'AI_60d_bin' : {
        'bin' : [0, 0.2, 0.5, 0.65, 0.8, 1., 1.5, 2., 1000.],
        'bin_labels' : ['0-0.2','0.2-0.5','0.5-0.65','0.65-0.8','0.8-1.0','1.0-1.5','1.5-2.0','>2.0'],
        'edges' : [0, 0.2, 0.5, 0.65, 0.8, 1., 1.5, 2., 3.2],
        'bin2' : [0, 0.25, 0.5, 1., 2., 1000.],
        'bin2_labels' : ['0-0.25','0.25-0.5','0.5-1.0','1.0-2.0','>2.0'],
        'edges2' : [0, 0.25, 0.5, 1., 2., 3.5],
        # 'bin2' : [0, 0.05, 0.2, 0.5, 0.65, 0.8, 1, 1.5, 2, 1000.],
        # 'bin2_labels' : ['0-0.05','0.05-0.2','0.2-0.5','0.5-0.65','0.65-0.8','0.8-1.0','1.0-1.5','1.5-2.0','>2.0'],
        # 'edges2' : [0, 0.05, 0.2, 0.5, 0.65, 0.8, 1., 1.5, 2., 4.],
        # 'bin2' : [0,0.05,0.2,0.5,0.65,0.8,1., 1000.],
        # 'bin2_labels' : ['0-0.05','0.05-0.2','0.2-0.5','0.5-0.65','0.65-0.8', '0.8-1.0', '>1.0'],
        'label' : 'AI (preceding 60 days)',
        'cmap' : cmap
    },
    'AI_bin' : {
        'bin' : [0, 0.2, 0.5, 0.65, 0.8, 1., 1000.],
        'bin_labels' : ['0-0.2','0.2-0.5','0.5-0.65','0.65-0.8','0.8-1.0','>1.0'],
        'bin2' : [0, 0.05, 0.2, 0.5, 0.65, 0.8, 1.0, 1000.],
        'bin2_labels' : ['0-0.05','0.05-0.2','0.2-0.5','0.5-0.65','0.65-0.8', '0.8-1.0', '>1.0'],
        'label' : 'Aridity Index',
        'cmap' : cmap
    },
    'AI_rel_bin' : {
        'bin' : [0, 0.2, 0.5, 0.65, 0.8, 1., 1.5, 2., 1000.],
        'bin_labels' : ['0-0.2','0.2-0.5','0.5-0.65','0.65-0.8','0.8-1.0','1.0-1.5','1.5-2.0','>2.0'],
        'label' : 'Relative Aridity Index',
        'cmap' : cmap
    },
    'AIi_60d_bin' : {
        'bin' : [0, 0.5, 1/1.5, 1., 1.25, 1.5, 2., 5., np.inf],
        'bin_labels' : ['0-0.5','0.5-0.67','0.67-1.0','1.0-1.25','1.25-1.5','1.5-2.0','2.0-5.0','>5.0'],
        'edges' : [0.1, 0.5, 1/1.5, 1., 1.25, 1.5, 2., 5., 15],
        'bin2' : [0, 0.5, 1., 2., 4., np.inf],
        'bin2_labels' : ['0-0.5','0.5-1.0','1.0-2.0','2.0-4.0','>4.0'],
        'edges2' : [0.1, 0.5, 1., 2., 4., 12.],
        # 'bin2' : [0, 0.05, 0.2, 0.5, 0.65, 0.8, 1, 1.5, 2, 1000.],
        # 'bin2_labels' : ['0-0.05','0.05-0.2','0.2-0.5','0.5-0.65','0.65-0.8','0.8-1.0','1.0-1.5','1.5-2.0','>2.0'],
        # 'edges2' : [0, 0.05, 0.2, 0.5, 0.65, 0.8, 1., 1.5, 2., 4.],
        # 'bin2' : [0,0.05,0.2,0.5,0.65,0.8,1., 1000.],
        # 'bin2_labels' : ['0-0.05','0.05-0.2','0.2-0.5','0.5-0.65','0.65-0.8', '0.8-1.0', '>1.0'],
        'label' : 'AI (preceding 60 days)',
        'cmap' : cmap
    },
    'AIi_bin' : {
        'bin' : [0, 1., 1.25, 1.5, 2., 5., np.inf],
        'bin_labels' : ['0-1.0','1.0-1.25','1.25-1.5','1.5-2.0','2.0-5.0','>5.0'],
        'bin2' : [0, 0.05, 0.2, 0.5, 0.65, 0.8, 1.0, 1000.],
        'bin2' : [0, 1., 1.25, 1.5, 2., 5., 20., np.inf],
        'bin2_labels' : ['0-1.0','1.0-1.25','1.25-1.5','1.5-2.0','2.0-5.0','5.0-20.0','>20.0'],
        'label' : 'Aridity Index',
        'cmap' : cmap
    },
    'AIi_rel_bin' : {
        'bin' : [0, 0.5, 1/1.5, 1., 1.25, 1.5, 2., 5., np.inf],
        'bin_labels' : ['0-0.5','0.5-0.67','0.67-1.0','1.0-1.25','1.25-1.5','1.5-2.0','2.0-5.0','>5.0'],
        'label' : 'Relative Aridity Index',
        'cmap' : cmap
    },
    'z_AI_60d_bin' : {
        'bin' : [-1, -0.5, -0.25, 0, 0.25, 0.5, 1., 100.],
        'bin_labels' : ['<-0.5','-0.5 - -0.25','-0.25-0','0-0.25','0.25-0.5','0.5-1', '>1'],
        'edges' : [-1, -0.5, -0.25, 0, 0.25, 0.5, 1., 2.],
        'label' : 'AI_{60} (z-score)',
        'cmap' : cmap,
    },
    'z_LAI_bin' : {
        'bin' : [-1, -0.5, -0.25, 0, 0.25, 0.5, 1., 100.],
        'bin_labels' : ['<-0.5','-0.5-0.25','-0.25-0','0-0.25','0.25-0.5','0.5-1', '>1'],
        'edges' : [-1, -0.5, -0.25, 0, 0.25, 0.5, 1., 3.],
        'label' : 'LAI (z-score)',
        'cmap' : 'summer_r',
    },
    'LAI_bin' : {
        'bin' : [0,0.5,1,2,3,5,10],
        'bin_labels' : ['0-0.5','0.5-1','1-2','2-3','3-5','>5'],
        'edges' : [0, 0.5, 1, 2, 3, 5, 7],
        'bin2' : [0,0.5,1,2,3,4,10],
        'bin2_labels' : ['0-0.5','0.5-1','1-2','2-3','3-4','>4'],
        'edges2' : [0, 0.5, 1, 2, 3, 4, 7],
        'label' : r'LAI (m$^2$ m$^{-2}$)',
        'cmap' : 'summer_r'
    },
    'MAP_bin' : {
        'bin' : [0,300,600,900,1200,1500, 6000],
        'bin_labels' : ['0-300','300-600','600-900','900-1200','1200-1500','>1500 mm'],
        'bin_log' : [0,125, 250, 500, 1000, 2000, 6000],
        'bin_log_labels' : ['0-125','125-250','250-500','500-1000','1000-2000','>2000 mm'],
        'label' : 'MAP (mm)',
        # 'cmap' : cmap
        'cmap' : c_list1[::-1]
    },
    'frequency_bin' : {
        'bin' : [0, 0.1, 0.2, 0.3, 0.4, 1.],
        'bin_labels' : ['0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','>0.4'],
        'bin_mth' : [0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.],
        'bin_mth_labels' : ['0-0.1','0.1-0.2','0.2-0.3','0.3-0.4', '0.4-0.5', '>0.5'],
        'label' : r'Rainfall frequency ($\lambda$, d$^{-1}$)',
        'cmap' : 'YlGnBu'
    },
    'intensity_bin' : {
        'bin' : [0, 5, 10, 15, 20, 100],
        'bin_labels' : ['0-5','5-10','10-15','15-20','>20'],
        'bin_mth' : [0, 5, 10, 15, 20, 25, 100],
        'bin_mth_labels' : ['0-5','5-10','10-15','15-20', '20-25', '>25'],
        'label' : r'Rainfall intensity ($\alpha$, mm d$^{-1}$)',
        'cmap' : 'YlGnBu'
    },
    'coeff_var_bin' : {
        'bin' : [0, 0.2, 0.25, 0.3, 0.35, 0.4, 1.],
        'bin_labels' : ['0-0.2','0.2-0.25','0.25-0.3','0.3-0.35','0.35-0.4','>0.4'],
        'bin_mth' : [0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 100.],
        'bin_mth_labels' : ['0-0.4','0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.9-1.0','>1.0'],
        'label' : 'Rainfall coefficient of variation',
        'cmap' : 'YlGnBu'
    },
    'perc_woodyveg_bin' : {
        # 'bin' : [0, 10, 20, 30, 40, 50, 60, 70 , 100],
        # 'bin' : [0, 10, 20, 30, 40, 50, 60, 100],
        'bin' : [0, 10, 20, 30, 40, 50, 100],
        'bin_labels' : [
            # '0-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%','80-90%','>90%'
            # '0-10%','10-20%','20-30%','30-40%','40-50%','50-60%','>60%'
            '0-10%','10-20%','20-30%','30-40%','40-50%','>50%'
        ],
        'edges' : [0, 10, 20, 30, 40, 50, 70],
        # 'bin2' : [0, 10, 20, 30, 40, 60, 100],
        # 'bin2_labels' : [
        #     '0-10%','10-20%','20-30%','30-40%','40-60%','>60%'
        # ],
        'bin2' : [0,10,20,30,40,60,100],
        'bin2_labels' : [
            '0-10%','10-20%','20-30%','30-40%','40-60%','>60%'
        ],
        'edges2' : [0, 10, 20, 30, 40, 60, 80],
        'label' : 'Woody cover (%)',
        'cmap' : 'summer_r'
    },
    'perc_woody_bin' : {
        'bin' : [0, 10, 20, 30, 40, 50, 100],
        'bin_labels' : ['0-10%','10-20%','20-30%','30-40%','40-50%','>50%'],
        'bin2' : [0, 10, 20, 30, 40, 50, 60, 100],
        'bin2_labels' : ['0-10%','10-20%','20-30%','30-40%','40-50%','50-60%','>60%'],
        'label' : 'Woody cover (% of vegetated)',
        'cmap' : 'summer_r'
    }
}
bin_dict['iperc_woodyveg_bin'] = bin_dict['perc_woodyveg_bin'].copy()
bin_dict['iperc_woodyveg_avg_bin'] = bin_dict['perc_woodyveg_bin'].copy()


cax_dict = {
    0 : [0.95, 0.2, 0.01, 0.68],
    2 : [0.915, 0.08, 0.023, 0.88],
    3 : [0.915, 0.15, 0.02, 0.78],
    5 : [0.95, 0.2, 0.01, 0.68],
    6 : [0.915, 0.08, 0.023, 0.88],
}

adj_dict = {
    0 : {'left':0.031, 'right':0.94, 'bottom':0.145, 'top':0.925, 'wspace':0.105},
    2 : {'left':0.054, 'right':0.895, 'bottom':0.085, 'top':0.957, 'wspace':0.1},
    3 : {'left':0.058, 'right':0.895, 'bottom':0.155, 'top':0.922, 'wspace':0.1},
    5 : {'left':0.031, 'right':0.94, 'bottom':0.145, 'top':0.925, 'wspace':0.105},
    6 : {'left':0.054, 'right':0.895, 'bottom':0.085, 'top':0.957, 'wspace':0.1},
}

size_dict = {
    2 : (6.25, 2.5),
    3 : (6.25, 2.2),
    5 : (12, 2.5),
    6 : (6.25, 4),
}

hue_dict = bin_dict.copy()
hue_dict.update({
    'texture' : {
        'hue_order' : [
            'Silty Clay', 'Clay', 'Clay Loam', 'Silt Loam', 'Silty Clay Loam', 
            'Sandy Clay', 'Loam', 'Sandy Clay Loam', 'Sandy Loam', 'Loamy Sand', 'Sand'
        ],
        'cmap' : 'husl'
    }
})


def get_palette(by):
    if 'IGBP' in by:
        hue_order = igbp_cats.get(
            re.findall(r"(?<=_)[^_]*(?!.*_)", by)[0] if '_' in by else by, 'og'
        )
        pal = [igbp_dict[k]['color'] for k in hue_order]
    elif 'bin' in by:
        # by = by[1:] if by[0] == 'i' else by
        hue_order = bin_dict.get(re.findall(r'.*?bin', by)[0])[f'{re.findall(r'bin.*$', by)[0]}_labels']
        n = len(hue_order)
        cmap1 = bin_dict.get(re.findall(r'.*?bin', by)[0])['cmap']
        if isinstance(cmap1, list):
            cmap1 = cmap1
        elif isinstance(cmap1, str):
            cmap1 = list(sns.color_palette(cmap1, n_colors=n).as_hex()) 
        else:
            cmap1 = [cmap1(plt.Normalize(0, n-1)(i)) for i in range(n)]    # 'Blues'
        pal = cmap1
    else:
        hue_order = hue_dict[by]['hue_order']
        hue_order = hue_dict.get(by, {}).get(
            'hue_order', 
            pd.Series(df[by].unique()).dropna().sort_values(ignore_index=True).tolist()
        )
        cmap = hue_dict.get(by, {}).get('cmap', 'husl')
        pal = list(sns.color_palette(cmap, n_colors=len(hue_order)).as_hex())
    return hue_order, pal

def get_fig_dets(n, width=None):
    fig_size = size_dict.get(n, (n*2.4, 2.5))
    if width:
        fig_size = (width, width/(fig_size[0]/fig_size[1]))
    if n > 5 and n % 2 == 0:
        r, c = 2, n//2
    else:
        r, c = 1, n
    cax_dims = cax_dict.get(n, [0.95, 0.2, 0.01, 0.68])
    adj_kwargs = adj_dict.get(n, {'left':0.031, 'right':0.94, 'bottom':0.145, 'top':0.925, 'wspace':0.105})
    return fig_size, r, c, cax_dims, adj_kwargs


# PLOT FUNCTIONS
def plot_box(ax, data, c, positions, flierprops, width=0.75, **kwargs):

    bp = ax.boxplot(
        data, patch_artist=True,
        positions=positions,
        flierprops=flierprops,
        widths=width, **kwargs
    )
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=c)

    for patch in bp['boxes']:
        patch.set(facecolor=c, alpha=0.75)       


    return bp

def boxplots(ax, data, y, x, by=None, var_list=None, color_list=None, spacing=2, other_vals=None, **kwargs):
    if not by:
        x_list = [data[data[x] == val][y] for val in var_list]
        positions = np.arange(len(x_list))
        bps = [plot_box(
            ax=ax, data=x_list[i], c=color_list[i], positions=[positions[i]], 
            flierprops=dict(color=color_list[i], marker='.', markeredgecolor=color_list[i]),
            **kwargs
        ) for i in np.arange(len(x_list))]
        xticks = np.arange(len(x_list))

    else:
        bps = []
        for i,var in enumerate(var_list):
            x_list = [data[(data[by] == var) & (data[x] == val)][y] for val in data[x].unique()]

            c = color_list[i]
            positions = [i+j*(len(var_list) + spacing) for j in np.arange(len(x_list))]

            bp = plot_box(
                ax=ax, data=x_list, c=color_list[i], positions=positions, 
                flierprops=dict(color=c, marker='.', markeredgecolor=c), **kwargs
            )
            if other_vals is not None:
                ax.plot(positions, other_vals[i], '.', markersize=7, color='k')

            bps.append(bp)
        
    # xtick1 = np.mean(np.arange(len(var_list)))
        xticks = [
            np.mean(np.arange(len(var_list))) + j*(
                len(var_list) + spacing) for j in np.arange(len(x_list))
        ]
    ax.set_xticks(xticks)

    buffer = xticks[0] + spacing

    ax.set_xlim(xticks[0]-buffer, xticks[-1]+buffer)
    
    return bps

def boxplot(
        ax, data, y, x, by, var_list, color_list, spacing=2, multi=True,
        other_vals=None, **kwargs
    ):

    bps = []
    for i,var in enumerate(var_list):
        # x_list = [data[(data[by] == var) & (data[x] == val)][y] for val in data[x].unique()]

        c = color_list[i]
        if multi:
            x_list = [data[(data[by] == var) & (data[x] == val)][y] for val in data[x].unique()]
            positions = [i+j*(len(var_list) + spacing) for j in np.arange(len(x_list))]
        else:
            x_list = [data[(data[by] == var)][y]]
            positions = [i]
        

        bp = plot_box(
            ax=ax, data=x_list, c=color_list[i], positions=positions, 
            flierprops=dict(color=c, marker='.', markeredgecolor=c),
            **kwargs
        )
        if other_vals is not None:
            ax.plot(positions, other_vals[i], '.', markersize=7, color='k')

        bps.append(bp)
    
    # xtick1 = np.mean(np.arange(len(var_list)))
    xticks = [
        np.mean(np.arange(len(var_list))) + j*(
            len(var_list) + spacing) for j in np.arange(len(x_list))
    ]

    buffer = xticks[0] + spacing if multi else 0

    ax.set_xticks(xticks)

    ax.set_xlim(xticks[0]-buffer, xticks[-1]+buffer)
    
    return bps

def plot_heatmaps(
        df, axs, x, y, z, by, cmap, cax, vmin=None, vmax=None,
        fs=9, adj_kwargs={'left' : 0.1, 'right' : 0.9, 'bottom' : 0.1, 'top' : 0.9},
        annot_kwargs={'annot' : True, 'annot_kws' : {'size' : 6.7}},
        contours=False, contour_kwargs={'levels' : 10}, contour_label=True,
        xlab='AI (preceding 60 days)', ylab='LAI', clab=None,
        xticks=None, yticks=None, cticks=([0.5, 1.0, 1.5, 2.0], [0.5, 1.0, 1.5, 2.0]),
        aggfunc='median', mask_count=0, invert_x=True, invert_y=True, ticks=False,
    ):

    xvals = pd.Series(df[x].unique()).sort_values(ignore_index=True).dropna()
    yvals = pd.Series(df[y].unique()).sort_values(ignore_index=True).dropna()[1:]
    X, Y = np.meshgrid(xvals, yvals)

    if 'IGBP' in by:
        by_order = igbp_cats.get(
                re.findall(r"(?<=_)[^_]*(?!.*_)", by)[0] if '_' in by else by, igbp_cats['og']
        )
    else:
        by_order = bin_dict.get(re.findall(r'.*?bin', by)[0])[f'{re.findall(r'bin.*$', by)[0]}_labels']


    for i,(val, ax) in enumerate(zip(by_order, axs.flatten())):
        dfi = df[df[by] == val]
        
        data = dfi.pivot_table(index=y, columns=x, values=z, aggfunc=aggfunc)
        data = data.reindex(
            index=yvals,
            columns=xvals,
        )

        data1 = dfi.pivot_table(index=y, columns=x, values=z, aggfunc='count')
        data1 = data1.reindex(
            index=yvals, columns=xvals,
        )
        # ax1.plot_surface(
        #     X, Y, Z,
        #     alpha=0.5, color=pal[i]
        # )
        if mask_count:
            data = data.mask(data1 < mask_count)

        im =  sns.heatmap(
            data,
            ax=ax, 
            vmin=vmin, vmax=vmax,
            cbar=i == 0, cbar_ax=None if i else cax,
            # cbar = i == len(by_order)-1, 
            # cbar_ax=cax, #None if i == len(by_order)-1 else cax,
            cmap=cmap,
            **annot_kwargs,
        )
        if contours:
            cont = ax.contour(
                np.arange(0.5, len(xvals)), np.arange(0.5, len(yvals)), data, 
                colors='k', linewidths=0.5, **contour_kwargs
            )
            if contour_label:
                ax.clabel(cont, cont.levels, inline=True, fmt='%1.1f', fontsize=fs-1)
        ax.set_ylabel("")
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_xticklabels([])

        if invert_x:
            ax.invert_xaxis()
        if invert_y:
            ax.invert_yaxis()
        try:
            ax.set_title(igbp_dict[val]['label'], fontsize=fs, pad=5)
        except:
            ax.set_title(val, fontsize=fs, pad=5)

        if not ticks:
            ax.tick_params(left=False, bottom=False, pad=-0.2)
        else:
            ax.tick_params(width=0.5, length=2,)

    r,c = axs.shape if len(axs.shape) > 1 else (1, axs.shape[0])
    # r = r[0] if r else 1

    if r == 1:
        yaxs = [axs[0]]
        xaxs = axs
    elif r == 2:
        yaxs = axs[:,0]
        xaxs = axs[1,:]

    if xticks is None:
        xticks = get_ticks(df, x, np.max(ax.get_xlim()))
    if yticks is None:
        yticks = get_ticks(df, y, np.max(ax.get_ylim()))
        if yticks[0][0] < 0:
            yticks = (yticks[0][1:], yticks[1][1:])
        else:
            yticks = (yticks[0][:-1], yticks[1][1:])

    for ax in yaxs:
        ax.set_ylabel(ylab, labelpad=2, fontsize=fs)
        ax.set_yticks(
            *yticks, rotation=0, ha='right', rotation_mode='anchor', fontsize=fs-0.5
        )
    for ax in xaxs:
        ax.set_xlabel(xlab, fontsize=fs)
        ax.set_xticks(
            *xticks, rotation=0, ha='center', fontsize=fs-0.5
        )
    # cb = plt.colorbar(im, cax=cax, orientation='vertical')
    # cb.outline.set_visible(False)

    plt.subplots_adjust(
        **adj_kwargs
    )
    cax.set_ylabel(clab, fontsize=fs, labelpad=2)
    cax.yaxis.set_tick_params(labelsize=fs-0.5, pad=2)
    cax.set_yticks(*cticks)

    return axs, cax




def label_line(
    line, label_text, near_i=None, near_x=None, near_y=None, rotation_offset=0, offset=(0,0),
    fontsize=9., loc=('center', 'center')
):
    """call 
        l, = plt.loglog(x, y)
        label_line(l, "text", near_x=0.32)
    """
    def put_label(i):
        """put label at given index"""
        i = min(i, len(x)-2)
        dx = sx[i+1] - sx[i]
        dy = sy[i+1] - sy[i]
        rotation = np.rad2deg(np.arctan2(dy, dx)) + rotation_offset
        pos = [(x[i] + x[i+1])/2. + offset[0], (y[i] + y[i+1])/2 + offset[1]]
        ax.text(
            pos[0], pos[1], label_text, size=fontsize, rotation=rotation, transform_rotates_text=True,
            color = line.get_color(),
            ha=loc[0], va=loc[1], bbox = dict(ec='1',fc='1', pad=0.3)
        )

    x = line.get_xdata()
    y = line.get_ydata()
    ax = line.axes
    if ax.get_xscale() == 'log':
        sx = np.log10(x)    # screen space
    else:
        sx = x
    if ax.get_yscale() == 'log':
        sy = np.log10(y)
    else:
        sy = y

    # find index
    if near_i is not None:
        i = near_i
        if i < 0: # sanitize negative i
            i = len(x) + i
        put_label(i)
    elif near_x is not None:
        for i in range(len(x)-2):
            if (x[i] < near_x and x[i+1] >= near_x) or (x[i+1] < near_x and x[i] >= near_x):
                put_label(i)
    elif near_y is not None:
        for i in range(len(y)-2):
            if (y[i] < near_y and y[i+1] >= near_y) or (y[i+1] < near_y and y[i] >= near_y):
                put_label(i)
    else:
        raise ValueError("Need one of near_i, near_x, near_y")
    
