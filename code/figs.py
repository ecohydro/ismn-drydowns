#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Wed 01 May 24 14:52:16'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           figs.py
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


#%% IMPORTS

import os
import sys
import re
import numpy as np
import pandas as pd
import datetime
import pytz
import cartopy.crs as ccrs
import cartopy
import itertools
from scipy import stats
from scipy.optimize import curve_fit

from xarray import DataArray


import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib import rcParams
from mpl_toolkits import mplot3d

# %matplotlib qt

import seaborn as sns


from utils_figs import igbp_dict, var_dict, get_continuous_cmap, get_ticks, get_fig_dets
from utils_figs import igbp_cats, bin_dict, cmap, plot_heatmaps, get_palette, label_line, rescale, boxplots

# from station import config, META, Station, ismn_data
from drydown import Drydown, calc_dtheta_dt, calc_theta_q, calc_theta_exp
from soil.station import config, META, Station#, ismn_data
import drydowns
from drydowns.config import Config

import utils_figs as utils
from utils_read import read_output, clean_data

#%% PLOTTING PARAMS

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Myriad Pro'
# mpl.rcParams['font.sans-serif'] = 'Source Sans Pro'
mpl.rcParams['font.size'] = 10.0
mpl.rcParams['axes.titlesize'] = 10.0
awl = 0.6
for p in ['axes.linewidth', 'xtick.major.width', 'ytick.major.width']:
    mpl.rcParams[p] = awl

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
mpl.rcParams['mathtext.fontset'] = 'custom'


cmap = get_continuous_cmap(['#F0CC47', '#82B960', '#5BB8B0',  '#0D5654'])
cmap = get_continuous_cmap(['#F0CC47', '#0D5654'])

c_list = ['#1B3A39', '#5BB8B0', '#82B960', '#B9C354', '#EDC01D']
c_listr = ['#EDC01D', '#82B960', '#5BB8B0', '#1B3A39']
cmap = get_continuous_cmap(c_listr)


#%% PATHS

config = Config('soil/config.ini')
proj_dir = config.config['ISMN']['project_dir']   
data_dir = os.path.join(proj_dir, 'data')
res_dir = os.path.join(proj_dir, 'outputs')

chirps_dir = os.path.join(res_dir, 'chirps')

fig_dir = os.path.join(proj_dir, 'results') # 'agg_curves_igbp')
# fig_dir = '/Users/brynmorgan/dev/dissertation/figs/ch4'
# fig_dir = '/Users/brynmorgan/dev/ismn-drydowns-manuscript/figs'
fig_dir = os.path.join(proj_dir, 'results', 'gpp')

site_cols = [
    'network', 'station', 'latitude', 'longitude', 'elevation', 
    'MAP', 'MAT', 'AI', 'climate_KG', 
    'LAI_max', 'LAI_mean',
    'iperc_woodyveg', 'iperc_woodyveg_avg',
    *[f'{pre}{suff}' for pre in [
        'CCI_IGBP','iIGBP','iIGBP_ed','mIGBP'
    ] for suff in ['','_5','_6s','_6b','_7']], 'iIGBP_3', 'mIGBP_3'
]


#%% DATA IMPORT

file = '../data/ismn_results_star_th0_04may.csv'

ismn_t04 = read_output(file)

raw_dfs = {
    'ismn_t04' : ismn_t04,
}


filt_dfs = {
    k : clean_data(
            v, r2_thresh=0.7, duration=(6, 100), rfrac=0.2
        ) for k,v in raw_dfs.items()
}   

key = 'ismn_t04'
df = filt_dfs[key]


#%% FUNCTIONS


exp_cols = [col for col in df.columns if col.startswith('exp_') and '_opt' not in col]
q_cols = [col for col in df.columns if col.startswith('q_') and '_opt' not in col]
# cols = ['duration', 'min_sm', 'max_sm', 'theta_fc', 'pet'] + exp_cols + q_cols
cols = ['duration', 'z_m', 'theta_fc', 'min_sm', 'max_sm', 'pet'] + exp_cols + q_cols 


def get_data(df, i, by, pre='q_'):
    
    data = df.iloc[i].to_dict()
    
    # theta_w = data['min_sm']
    theta_w = data.get(f'theta_w', data['min_sm'])
    theta_star = data.get(f'{pre}theta_star', data['theta_fc'])
    theta = np.linspace(theta_star, theta_w, 25)

    dtheta = calc_dtheta_dt(
        theta = theta, k = data[f'{pre}k'], q = data[f'{pre}q'],
        theta_w = theta_w, theta_star = theta_star
    )
    # dtheta_mm = dtheta*0.05*1000
    dtheta_mm = calc_dtheta_dt(
        theta = theta, k = data[f'{pre}ET_max'], q = data[f'{pre}q'],
        theta_w = theta_w, theta_star = theta_star,
    )
    dtheta_lin = calc_dtheta_dt(
        theta = theta, k = data[f'{pre}k'], q = 1.0,
        theta_w = theta_w, theta_star = theta_star,
    )
    dtheta_mm_lin = calc_dtheta_dt(
        theta = theta, k = data[f'{pre}ET_max'], q = 1.0,
        theta_w = theta_w, theta_star = theta_star,
    )
    out = {
        'theta' : theta,
        's' : theta / data['theta_fc'],
        's_norm' : (theta - theta_w)/(data['theta_fc'] - theta_w),
        'dtheta' : dtheta,
        'dtheta_mm' : dtheta_mm,
        'theta_norm' : (theta - theta_w)/(theta_star - theta_w),
        'dtheta_mm_norm' : dtheta_mm/data[f'{pre}ET_max'],
        'dtheta_norm' : dtheta/data[f'{pre}k'],
        'dtheta_lin' : dtheta_lin,
        'dtheta_mm_lin' : dtheta_mm_lin,
        'q' : data[f'{pre}q'],
    }
    out = df.loc[i, by].to_dict() | out
    df_out = pd.DataFrame(out)

    return df_out

def get_agg_data(df, by, how='median'):
    dfg = df.groupby(by=by)[cols].agg(how).reset_index()
    df_out = pd.concat(
        [get_data(dfg, i, by) for i, row in dfg.iterrows()]
    )
    return df_out

def get_site_data(
    df, by=site_cols, var='q_q', how=['median','count'], 
    col_dict={'median':'q_q', 'count':'n_events'},
    crs=ccrs.PlateCarree(), 
):
    dfg = df.groupby(by=by)[var].agg(how).reset_index()
    dfg = dfg.rename(columns=col_dict)

    dfg = gpd.GeoDataFrame(
        dfg, geometry=gpd.points_from_xy(dfg.longitude, dfg.latitude)
    ).set_crs(epsg=4326)

    dfg = dfg.to_crs(crs.proj4_init)

    return dfg

def mask_data(df, col_dict={'LAI' : 0.5, 'MAP' : 1200.}):
    mask = df
    for k,v in col_dict.items():
        if 'LAI' in k:
            mask = mask[mask[k] >= v]
        elif 'season' in k:
            mask = mask[mask[k] == v]
        else:
            mask = mask[mask[k] <= v] 
    return mask

def annotate_vline(
    ax, text, x, ylocs, x_buff=0, y_buff=0, color='k', va='bottom',
    line_kwargs={'linestyle':'--', 'linewidth':0.5},
    text_kwargs={'fontsize': 8.5, 'ha':'center',}
):
    if va == 'bottom':
        text_kwargs['va'] = 'bottom'
        ytext = ylocs[1] + y_buff
    else:
        text_kwargs['va'] = 'top'
        ytext = ylocs[0] - (y_buff + 0.02)
    
    ax.vlines(x=x, ymin=ylocs[0], ymax=ylocs[1], color=color, **line_kwargs)
    ax.text(
        x=x+x_buff, y=ytext, s=text, color=color, **text_kwargs
    )

# def annotate_thetas(ax, dfi, x='theta', y='dtheta_mm', vals=[0], ys=[0.35], x_buff=0):
#     x_mins = dfi.groupby(by)[x].min().values[vals]
#     x_maxes = dfi.groupby(by)[x].max().values[vals]
#     y_maxes = dfi.groupby(by)[y].max().values[vals]

def annotate_curves(
    axs, dfi, pal, theta_vals=[0], q_vals=[0, -1], 
    x_buff=0, y_buff=0.03, x_norms=[0.3, 0.7]
):

    maxes = dfi.groupby(by)[['theta','dtheta_mm']].max()
    x_mins = dfi.groupby(by).theta.min()#.values[q_vals]
    qs = dfi.groupby(by).q.median()
    # x_maxes = dfi.groupby(by).theta.max().values[q_vals]
    # y_maxes = dfi.groupby(by).dtheta_mm.max().values[q_vals]

    ys = [0.35]


    for i,j in enumerate(theta_vals):
        theta_star = maxes['theta'].iloc[j]
        et_max = maxes['dtheta_mm'].iloc[j]
        # theta_w
        annotate_vline(
            ax=axs[0], text=r"$\theta_{\mathrm{w}}$", 
            x=x_mins.iloc[j], ylocs=[0, ys[i]], 
            x_buff=x_buff, y_buff=y_buff, color=pal[j]
        )
        # theta_star
        annotate_vline(
            ax=axs[0], text=r"$\theta_{\ast}$",
            # x=x_maxes[i], ylocs=[y_maxes[i]-ys[0], y_maxes[i]],
            x = theta_star, ylocs = [et_max - ys[0], et_max],
            x_buff=0, y_buff=y_buff, color=pal[j], va='top'
        )
        ys += [ys[i] + 0.2]


    for i,j in enumerate(q_vals):
        theta_star = maxes['theta'].iloc[j]
        et_max = maxes['dtheta_mm'].iloc[j]
        theta_w = x_mins.iloc[j]
        q = qs.iloc[j]
        y_mm = calc_dtheta_dt(
            theta = (x_norms[i]*(theta_star-theta_w) + theta_w), k = et_max, 
            q=q, theta_w=theta_w, theta_star=theta_star
        ) / et_max
        if q < 1:
            ylocs = (y_mm, y_mm + 0.3)
            va = 'bottom'
        else:
            ylocs = (y_mm - 0.3, y_mm)
            va = 'top'
        annotate_vline(
            ax=axs[1], text=f"$q = {q.round(2)}$", 
            x=x_norms[i], ylocs=ylocs, 
            x_buff=0, y_buff=0.02, color=pal[j], va=va
        )


def plot_curves(
        df, axs, by, how, hue_order, pal, leg_kwargs, adj_kwargs,
        letts=['(a)','(b)'], fs=10, annotate=False
    ):
    dfi = get_agg_data(df, by=[by], how=how)

    try:
        labs = [igbp_dict[k]['label'] for k in hue_order]
    except:
        labs = hue_order

    for ax, norm in zip(axs, ['','_norm']):
        x = 'theta' + norm
        y = 'dtheta_mm' + norm

        if norm == '_norm':
            ax.plot(
                [0,1], [0,1], 'k--', linewidth=0.75
            )

        # Plot absolute curves
        sns.lineplot(
            ax=ax,
            x = dfi[x], 
            y = dfi[y], 
            hue = dfi[by],
            hue_order= hue_order,
            palette = pal,
            legend=False if norm == '' else True,
            # label=labs,
        )
        lp = 3

        if norm == '_norm':
            # ax.set_xlabel(r"$\frac{\theta - \theta_w}{\theta_* - \theta_w}$", fontsize=fs+1)
            # ax.set_ylabel(r"$\frac{\mathrm{ET}}{\mathrm{ET_{max}}}$", fontsize=fs+1)
            ax.set_xlabel(r"$\frac{\theta - \theta_w}{\theta_{\ast} - \theta_w}$", fontsize=fs+1, labelpad=lp)
            ax.set_ylabel(r"$\frac{\mathrm{ET}}{\mathrm{ET_{max}}}$", fontsize=fs+1, labelpad=lp)
            ax.set_yticks([0, 0.5, 1])
            ax.set_xticks([0, 0.5, 1])
        else:
            # ax.set_xlabel(r"$\theta$ (m$^3$ m$^{-3}$)", fontsize=fs+1)
            ax.set_xlabel(r"$\theta$ (m$^3$ m$^{-3}$)", fontsize=fs, labelpad=lp)
            # ax.set_ylabel(r"$\frac{d\theta}{dt}$ (mm d$^{-1}$)", fontsize=fs+1)
            ax.set_ylabel(r"$\mathrm{ET}$ (mm d$^{-1}$)", fontsize=fs, labelpad=lp)

        ax.tick_params(axis='both', which='major', labelsize=fs)

    if letts:
        for lett,ax in zip(letts, axs):
            ax.text(
                0.03, 0.97, lett, 
                transform=ax.transAxes, 
                fontsize=fs, 
                verticalalignment='top'
            )

    axs[1].legend(
        **leg_kwargs
    )

    axs[1].get_legend_handles_labels()
    axs[1].legend(
        handles=axs[1].get_legend_handles_labels()[0],
        labels=labs,
        **leg_kwargs
    )

    plt.tight_layout()

    plt.subplots_adjust(
        **adj_kwargs
    )

    return axs

#%%

# FIG. 1. CONCEPTUAL --> see code/conceptual.py

# FIG. 2. MAPS --> see code/ternary.py


#%% FIG. 3. AGGREGATED CURVES


save = False
# bys = ['perc_woodyveg_bin', 'MAP_bin']
bys = ['MOD_IGBP_3', 'MAP_bin']

how = 'median'

df = df[
    df.MOD_IGBP_3.isin(['HERB','SAV','WOOD'])
    & (df.LAI >= 0.5)
]

ratio = 6.25/4.2
cm = 1/2.54

width = 11.4*cm
fig_size = (width, width/ratio)

leg_kwargs = {
    'frameon':False,
    'handletextpad':0.7,
    'handlelength':1.4,
    'bbox_to_anchor':(0.99, 1.02),
    'loc':'upper left',
    'fontsize':9,
}

# fig_size = (6.25, 4.2)  # dissertation

fs = 9.
# 4-panel figure
fig, axs = plt.subplots(2, 2, figsize=fig_size)

# adj_kwargs = {'left':0.088, 'right':0.83, 'bottom':0.128, 'top':0.999, 'wspace':0.3, 'hspace':0.1}
adj_kwargs = {'left':0.078, 'right':0.785, 'bottom':0.16, 'top':0.999, 'wspace':0.425, 'hspace':0.1}


for i,by in enumerate(bys):
    axsi = axs[i]

    hue_order, pal = get_palette(by)

    dfi = get_agg_data(df, by=[by], how=how)
    vals = [0]

    annotate_curves(
        axs=axsi, dfi=dfi, pal=pal, theta_vals=vals, x_buff= 0 if i == 0 else 0.0075
    )

    letts = ['(a)','(b)'] if i == 0 else ['(c)','(d)']
    # letts = ['A','B'] if i == 0 else ['C','D']

    axsi = plot_curves(
        df, axsi, by=by, how=how, hue_order=hue_order, pal=pal,
        leg_kwargs=leg_kwargs, adj_kwargs=adj_kwargs, letts=letts, fs=fs
    )
    axsi[0].set_xlim(0.00025, 0.33)
    axsi[0].set_xticks([0.1, 0.2, 0.3])

    if i == 0:
        for ax in axsi:
            ax.set_xlabel('')
            ax.set_xticklabels([])

if save:
    by_names = []        
    for by in bys:
        by_names.append(
            by.lower() if 'IGBP' in by else re.findall(r'^(.*?)(?:_bin|$)', by)[0].lower()
        )
    # suff = '_'.join([k.lower() for k in mask])
    # suff = '_' + suff if suff else ''
    suff = ''
    name = f'{how}_curves_{key[5:]}_{'_'.join(by_names)}{suff}' + '_sci'
    
    plt.savefig(os.path.join(fig_dir, f'{name}.pdf'), dpi=1200, transparent=True,) # bbox_inches='tight'
    plt.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=1200, transparent=True,) # bbox_inches='tight'


#%% FIG. 4. HEATMAPS (SEE )

# fig_dir = os.path.join(proj_dir, 'results') # 'agg_curves_igbp')
save = False

key = 'ismn_t04'
df = filt_dfs[key]

# df = df[df.season == 'wet']
# df = df[df.LAI >= 0.5]

col = 'MOD_IGBP_3'

aggfunc = 'median'

y = 'LAI_bin'
x = 'AI_60d_bin'

xlab = 'AI (preceding 60 days)'

df = df[df.MOD_IGBP.isin(['GRA','SHR','SAV','WSA','EBF','DBF','ENF','MF', 'BLF', 'OSH','CSH','MOS'])]


val = 'q_q'

if 'IGBP' in col:
    igbp_order = igbp_cats.get(
        re.findall(r"(?<=_)[^_]*(?!.*_)", col)[0] if '_' in col else col, igbp_cats['og']
    )
else:
    colk = col[1:] if col[0] == 'i' else col
    colk = colk[:-8]+colk[-4:] if colk[-8:-4] == '_avg' else colk
    igbp_order = bin_dict.get(re.findall(r'.*?bin', colk)[0])[f'{re.findall(r'bin.*$', colk)[0]}_labels']

n = len(igbp_order)
cm = 1/2.54
fig_size, r, c, cax_dims, adj_kwargs = get_fig_dets(n, width=17.8*cm)
# fig_size, r, c, cax_dims, adj_kwargs = get_fig_dets(n)
adj_kwargs = {'left':0.053, 'right':0.9025, 'bottom':0.135, 'top':0.922, 'wspace':0.1}
cax_dims = [0.925, 0.135, 0.02, 0.7875]

# cax_dims = [0.915, 0.155, 0.02, 0.765]

c_list1 = ['#004446','#026f6f', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
cmap = get_continuous_cmap(c_list1[::-1])

# cmap = 'BrBG_r'

fs = 9

def get_vrange(df, val=val, col='IGBP', aggfunc='median'):
    df_agg = df.groupby(by=col)[val].agg(aggfunc)
    return df_agg.min(), df_agg.max()

vmin, vmax = get_vrange(df, col=col, val=val, aggfunc=aggfunc)
# vmin, vmax = get_vrange(df, col=[col,x,y], val=val, aggfunc=aggfunc)
# vmin, vmax = (0, 200)
if val == 'q_q':
    vmin, vmax = (0.25, 2.2)
    cticks = ([0.5, 1.0, 1.5, 2.0], [0.5, 1.0, 1.5, 2.0])
    clab = f'{aggfunc.capitalize()} $q$'
    # vmin, vmax = (0.25, 1.75)
elif val == 'q_theta_0':
    vmin, vmax = (0.14, 0.28)
    cticks = ([0.15, 0.175, 0.2, 0.225, 0.25, 0.275], [0.15, 0.175, 0.2, 0.225, 0.25, 0.275])
    clab = f'{aggfunc.capitalize()} $\\theta_0$'
elif val == 'q_ET_max':
    vmin, vmax = (1.5, 3.)
    cticks = ([1.5, 2.0, 2.5, 3.0], [1.5, 2.0, 2.5, 3.0])
    clab = f'{aggfunc.capitalize()} $ET_{max}$'

inv_x = True if not 'AIi' in x else False


fig, axs = plt.subplots(r, c, figsize=fig_size, sharex=True, sharey=True)
cax = fig.add_axes(cax_dims)  # Adjust the position and size as needed

axs, cax = plot_heatmaps(
    df, axs, x, y, val, col, cmap, cax, vmin=vmin, vmax=vmax, cticks=cticks, aggfunc=aggfunc, fs=fs,
    clab=clab, xlab=xlab, ylab='LAI (m$^2$ m$^{-2}$)', mask_count=5, invert_x=inv_x
)

for ax,lett in zip(axs, ['a','b','c','d','e','f']):
    ax.text(0.0, 1.09, f'({lett})', transform=ax.transAxes, fontsize=fs, va='top')

plt.subplots_adjust(
    **adj_kwargs
)


if save:
    xystr = f'{re.findall(r'(.*?)_bin', x)[0]}_{re.findall(r'(.*?)_bin', y)[0]}'
    what = col.lower() #'igbp' + str(len(igbp_order))
    name = f'heatmap_{val[2:]}_{key[5:]}_{what}_{xystr}_{aggfunc}'  
    plt.savefig(os.path.join(fig_dir, f'{name}.pdf'), dpi=1200, transparent=True,) # bbox_inches='tight'
    plt.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=1200, transparent=True,) # bbox_inches='tight'



#%% 

# FIG. 4. HEATMAPS + SURFACES --> see code/surfaces.py

# FIG. 5. CONTOURS --> see code/surfaces.py




# %% APPENDIX FIGS

#%% FIG S1. COMPARISON BETWEEN LINEAR + NONLINEAR r2

save = False


mako = sns.color_palette("mako", 5).as_hex()


fig_size = (3.125, 2.75)

key = 'ismn_t04'
df = filt_dfs[key]


x = 'exp_r_squared'
y = 'q_r_squared'
by = 'MOD_IGBP_6s'

df = df[
    # (df.z_m <= 0.3) &
    # (df.MAP <= 1200.)
    (df.LAI >= 0.5)
    & (df[x].notna())
    & (df[y].notna())
]


by = 'MOD_IGBP_3'
by_suff = re.findall(r"(?<=_)[^_]*(?!.*_)", by)[0] if '_' in by else by

if 'IGBP' in by:
    hue_order = igbp_cats.get(
        by_suff, 'og'
    )
    pal = [igbp_dict[k]['color'] for k in hue_order]
# else:
#     hue_order = bin_dict.get(re.findall(r'.*?bin', by)[0])[f'{re.findall(r'bin.*$', by)[0]}_labels']
#     n = len(hue_order)
#     cmap1 = [cmap(plt.Normalize(0, n-1)(i)) for i in range(n)]    # 'Blues'
#     pal = cmap1

# c_list = ['#1B3A39', '#5BB8B0', '#82B960', '#B9C354', '#EDC01D']
# c_listr = ['#EDC01D', '#82B960', '#5BB8B0', '#1B3A39']
# cmap = get_continuous_cmap(c_listr)

df = df[df[by].isin(hue_order)]


leg_kwargs = {
    'fontsize' : 9,
    'frameon' : False,
    'handletextpad' : 0.5,
    # 'loc' : 'upper left',
    # 'bbox_to_anchor' : (1.0, 1)
}
adj_kwargs = {'left':0.19, 'right':0.98, 'top':0.995, 'bottom':0.170}


fig = plt.figure(figsize=fig_size)
ax = fig.add_subplot(1,1,1)

ax.plot(
    [0.5,1], [0.5,1], 'k--'
)

sns.scatterplot(
    x = df[x],
    y = df[y],
    # hue = df[by],
    # hue_order= hue_order,
    # palette = pal,
    # alpha=0.8, s=14
    # s=5, color='.15'
    c=np.log10(stats.gaussian_kde(df[[x,y]].T)(df[[x,y]].T)),
    s=5, cmap='mako',
    # hue_norm=(0, 0.005)
)
# Plot trendline
sns.regplot(
    x = df[x],
    y = df[y],
    scatter=False,
    # color='k',
    color=mako[1],
    line_kws={'linewidth':1},
    ci=None
)
# sns.histplot(x=df[x], y=df[y], bins=500, pthresh=0.9, cmap='mako')
# sns.kdeplot(x=df[x], y=df[y], levels=5, color='k', linewidths=0.5)

# ax.set_xlabel(r"$R^2_{\mathrm{exp}}$")
# ax.set_ylabel(r"$R^2_{\mathrm{q}}$")
ax.set_xlabel(r"$R^2_{\mathrm{linear}}$")
ax.set_ylabel(r"$R^2_{q}$")

# ax.set_xlim((0.56, 1.02))
# ax.set_ylim((0.81,1.0075))
ax.set_xlim((0.81, 1.005))
ax.set_ylim((0.89,1.0025))

# ax.set_yticks(np.arange(0.85, 1.01, 0.05))
ax.set_yticks(np.arange(0.9, 1.01, 0.05))

ax.legend(
    bbox_to_anchor=(0.98, 1), loc='upper left',
    **leg_kwargs
)

plt.tight_layout()
plt.subplots_adjust(**adj_kwargs)


if save:
    # by = f'igbp_{by_suff}' #+ str(len(hue_order))
    name = f'r2_{key[5:]}'
    plt.savefig(os.path.join(fig_dir, f'{name}.pdf'), dpi=1200, transparent=True,) # bbox_inches='tight'
    plt.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=1200, transparent=True,) # bbox_inches='tight'


# %% FIG S2. BOX PLOTS

save = False

# from utils_figs import boxplot

# fig_size = (6.25, 3)
cm = 1/2.54
width = 17.8*cm

fig_size = (width, width*0.48)

fs = 10

key = 'ismn_t04'
df = filt_dfs[key]

df1 = df


y = 'q_q'

adj_kwargs = {'left':0.055, 'right':0.995, 'bottom':0.237, 'top':0.94, 'wspace':0.076}

# hue_order = igbp_order

# fig = plt.figure(figsize=fig_size)
# ax1 = fig.add_subplot(1, 2, 1)
# ax2 = fig.add_subplot(1, 2, 2)

# fig, axs = plt.subplots(1, 2, figsize=fig_size, sharey=True)
# Set up gridspec
fig = plt.figure(figsize=fig_size)
gs = fig.add_gridspec(1, 2, width_ratios=[1, 2])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
# ax1 = axs[0]
# ax2 = axs[1]

ax1.hlines(1, -1, 3, color='k', linestyle='--', linewidth=0.5)
ax2.hlines(1, -1, 6, color='k', linestyle='--', linewidth=0.5)


x = 'MOD_IGBP_3'
by = x
hue = x
hue_order, pal = get_palette(by)
x_labs = ['Grasslands', 'Savannas', 'Woodlands']
df = df[df[by].isin(hue_order) & (df.MAP_bin.notna())]


bps = boxplots(
    ax1, df, y=y, x=x, var_list=hue_order, color_list=pal, spacing=1, width=0.65,# multi=False, #leg_kwargs=None, #adj_kwargs=adj_kwargs
)
by2 = 'MAP_bin'
hue_order2, pal2 = get_palette(by2)

bps2 = boxplots(
    ax2, df, y=y, x=by2, var_list=hue_order2, color_list=pal2, spacing=1, width=0.65 # multi=False, #leg_kwargs=None, #adj_kwargs=adj_kwargs
)


# ax1.set_xticklabels([igbp_dict[k]['label'] for k in hue_order], fontsize=9) #rotation=45, ha='right')
ax1.set_xticklabels(['Grasslands', 'Savannas', 'Woodlands'], fontsize=9, rotation=45, ha='right', rotation_mode='anchor')
ax1.tick_params(axis='x', pad=1.5)

for ax in [ax1,ax2]:
    ax.set_ylim((-1.5, 6.3))

ax1.set_yticks([0, 2, 4, 6], labels=[0, 2, 4, 6], fontsize=9)
ax1.set_xlim((-0.65, 2.65))

ax2.set_xlim((-0.65, 5.65))
ax2.set_yticklabels([])
ax2.set_xticklabels(hue_order2[:-1]+['>1500'], fontsize=9, rotation=45, ha='right', rotation_mode='anchor')
ax2.tick_params(axis='x', pad=1.5)
# ax2.set_xticks(np.arange(0, 6)-0.5, labels=['0','300','600','900','1200','>1500'], fontsize=9)
# ax2.set_xticklabels([f'{int(k.left)}-{int(k.right)}' for k in hue_order],) #rotation=45, ha='right')
ax1.set_ylabel(r"$q$")

ax1.set_xlabel("")
ax2.set_xlabel("MAP (mm)")
# plt.title('Median $q$ by site')
# plt.tight_layout()
plt.subplots_adjust(**adj_kwargs)

for ax,lett in zip([ax1,ax2], ['(a)','(b)']):
    ax.text(0.0, 1.08, lett, transform=ax.transAxes, fontsize=fs, va='top')

sig_list = ['$**$', '$**$', '++', '\n+', '++\n+', '$***$',]

for i,xval in enumerate(ax2.get_xticks()):
    # ax2.annotate(xval, -0.25, '*', ha='center', fontsize=9)
    ax2.annotate(sig_list[i], (xval, -0.25), ha='center', va='top', fontsize=8, xycoords='data')
for i, xval in enumerate(ax1.get_xticks()):
    ax1.annotate('$***$', (xval, -0.25), ha='center', va='top', fontsize=8, xycoords='data')

if save:
    name = f'boxplot_{key[5:]}'
    plt.savefig(os.path.join(fig_dir, f'{name}.pdf'), dpi=600, transparent=True,) # bbox_inches='tight'
    plt.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=600, transparent=True,) # bbox_inches='tight'


# %%
