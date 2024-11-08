
#%%
# import os
# import numpy as np
# import pandas as pd

# from scipy import stats

# import matplotlib.pyplot as plt
# import matplotlib as mpl
# from matplotlib import rcParams


import mpltern
import cartopy.io.shapereader as shpreader
import geopandas as gpd

import utils_figs as utils

# FUNCTIONS
#%%






fig_size = (6,5.25)
adj_kwargs = {
    'left' : 0.1, 'right' : 0.8, 'top' : 0.8,'bottom' : 0.1, 'wspace' : 0.1,
}

fs = 10.


#%%
# '#76B200', '#00B9A3','#006A58'

df = filt_dfs['ismn_t04']

save = False

# fig_size = (3.25, 3.25)
# fig_size = (3.85, 3.05)
fig_size = (3.5, 3)

adj_kwargs = {
    # 'left' : 0.0, 'right' : 1.0, 'top' : 0.997, 'bottom' : 0.003, 'wspace' : 0.1, 'hspace' : 0.1,
    # 'left' : 0.0, 'right' : 1.0, 'top' : 0.997, 'bottom' : 0.003, 'wspace' : 0.1, 'hspace' : 0.1,
    'left' : 0.0, 'right' : 1.0, 'top' : 0.997, 'bottom' : 0.003, 'wspace' : 0.1, 'hspace' : 0.1,
}

fs = 8.5

# width = 2.4
# fig_size = (6.25, 3.25)
# fig_size = (width,width/2.3)
# crs = ccrs.Mercator()
# proj_name = 'Mercator'
# crs = ccrs.Robinson()
# proj_name = 'Robinson'
# crs = ccrs.PlateCarree()
proj_name = 'PlateCarree'
crs = utils.map_dict.get(proj_name).get('crs') #if not us else ccrs.LambertConformal()


ratio = 6.25/2.95
cm = 1/2.54

width = 17.8*cm
fig_size = (width, width/ratio)

# fig_size = (6.25, 2.95)     # dissertation
# fig_size = (6.25, 2.83)
# Set up grid spec
fig = plt.figure(figsize=fig_size)
# gs = fig.add_gridspec(4, 3, width_ratios=[0.38, 0.3, 0.03], height_ratios=[0.09,0.4, 0.58,0.09], **adj_kwargs)
gs = fig.add_gridspec(4, 4, width_ratios=[0.555, 0.005, 0.426, 0.03], height_ratios=[0.15,0.4, 0.58,0.07], **adj_kwargs)
#   width_ratios=[0.57, 0.43], height_ratios=[0.58, 0.42], **adj_kwargs)



ax = fig.add_subplot(gs[1:3, 0], projection='ternary')

axm = fig.add_subplot(gs[:2, 2:], projection=crs)
ax_us = fig.add_subplot(gs[2:, 2], projection=ccrs.LambertConformal())

# ax_us = fig.add_subplot(gs[1, 1], projection=ccrs.LambertConformal())




hlz = df.groupby(by=site_cols).q_q.agg(['median','count']).reset_index()

hlz.rename(columns={
    'q_q' : 'n_events', 'median' : 'q_med', 'count' : 'n_events'
}, inplace=True)

hlz['MAT_ed'] = hlz.MAT.where((hlz.MAT > 0) & (hlz.MAT <=30.), 0.)


mask = (hlz['iIGBP_6s'].isin(['GRA','SHR','SAV','WSA','EBF','DBF','ENF','MF', 'BLF', 'OSH','CSH']))


# hue = 'iIGBP_6s'
hue = 'iIGBP_3'
# hue = 'q_med'
# hue_lab = '$q$'
# hue = 'iperc_woodyveg'
# hue_lab = 'Woody cover (%)'
size = 'LAI_max'

if 'IGBP' in hue:
    hue_order, pal = get_palette(hue)
    mask = (hlz[hue].isin(hue_order))
    hue_args = {
    'c' : hlz[mask][hue].apply(lambda x: igbp_dict[x]['color']), 
    }   
else:
    hue_order = None
    pal = cmap if hue == 'q_med' else 'summer_r'
    vmin, vmax = (0.25, 2.2) if hue == 'q_med' else (0., 80.)
    hue_args = {
        'c' : hlz[mask][hue],
        'cmap' : pal,
        'norm' : mpl.colors.Normalize(vmin=vmin, vmax=vmax),
    }


sizes, s_leg = utils.get_sizes(size, hlz, leg_sizes=(5, 40), **utils.pt_size_dict[size])
s_labs = utils.pt_size_dict[size]['s_labs']


def norm(x, labs=None):
    if labs is None:
        labs = x
    return (x - labs.min()) / (labs.max() - labs.min())

def rescale(x, xmin, xmax):
    return x * (xmax - xmin) + xmin

t_labs = np.array([0.1875, 1.5, 3, 6, 12, 24, 36])
# t_labs = np.array([0.09375, 1.5, 3, 6, 12, 24, 36])
t_divs = norm(-np.log10(t_labs))

map_labs = np.array([62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000])
# map_labs = np.array([31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000])
map_divs = norm(np.log10(map_labs))

pet_labs = np.array([0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32])
pet_divs = norm(np.log10(pet_labs))
ai_labs = 1/np.flip(pet_labs)
ai_divs = norm(-np.log10(ai_labs))

logmap = np.log10(hlz.MAP)
logmat = -np.log10(hlz.MAT_ed)
logpet = np.log10(1/hlz.AI)
logai = -np.log10(hlz.AI)

norm_map = norm(logmap, np.log10(map_labs))
norm_mat = norm(logmat, -np.log10(t_labs))
norm_pet = norm(logpet, np.log10(pet_labs))
norm_ai = norm(logai, -np.log10(ai_labs))

# norm_ai = 1 - norm_map - norm_mat
# norm_mat = 1 - norm_map - norm_pet



# ax = fig.add_subplot(projection='ternary')
# ax = fig.add_axes(ax1.get_position(), projection='ternary')
# fig.add_subplot(1, 1, 1, projection="ternary")
for t in t_divs[1:-1]:
    ax.axhline(
        t, xmin=-.1, xmax=1.1, clip_on=False,
        color='gray', alpha=0.3, linestyle='--', linewidth=0.75
    )
# ax.set_tlabel('Mean annual temperature (°C)')
# ax.set_llabel('PET ratio')
ax.set_llabel(r"Aridity Index $\left(\frac{P}{PET}\right)$", fontsize=fs,)
ax.set_rlabel('MAP (mm)', fontsize=fs)
ax.taxis.set_label_position('tick1')
ax.laxis.set_label_position('tick1')
ax.raxis.set_label_position('tick2')

ax.grid(
    linestyle='--', alpha=0.5, color='gray'
)

sc = ax.scatter(
    norm_mat[mask], 
    norm_ai[mask], 
    norm_map[mask], 
    alpha = 0.7, 
    edgecolor='none',
    s=sizes[mask],
    **hue_args
)


# Color ticks, grids, tick-labels
# ax.taxis.set_ticks(t_divs, t_labs)
ax.taxis.set_ticks([])
# ax.axis.set_ticks_position('tick')
# ax.laxis.set_ticks(pet_divs, pet_labs)
ax.laxis.set_ticks(ai_divs, [f'{i:.3f}'.rstrip('0').rstrip('.') for i in ai_labs], fontsize=fs)
ax.raxis.set_ticks(map_divs, [f'{i:.2f}'.rstrip('0').rstrip('.') for i in map_labs], fontsize=fs)
ax.raxis.set_ticks_position('tick2')

if 'IGBP' in hue:
    ax.legend(
        handles=[mpl.lines.Line2D(
            [], [], linestyle='None', marker='o', markeredgewidth=0, markerfacecolor=igbp_dict[k]['color'], 
        ) for k in hue_order],
        labels=[igbp_dict[k]['label'] for k in hue_order],
        bbox_to_anchor = (1.25,1.1), loc = 'upper right',
        prop = {'size': fs},
        frameon = False,
        framealpha = 0,
        handletextpad = 0.25,
        handlelength = 1, 
    )
else:
    cbar = fig.colorbar(sc, ax=ax, fraction=0.02, pad=-0.04)
    cbar.set_label(hue_lab, fontsize=fs)
    cbar.ax.tick_params(labelsize=fs, width=0.5)
    cbar.outline.set_visible(False)
    cbar.ax.set_anchor((0.0, 1.0))
    cbar.ax.set_box_aspect(15)

    if hue == 'q_med': #or 'iperc_woodyveg' in hue:    
        cbar.set_ticks(np.arange(0.5, 2.1, 0.5))
        cbar.set_ticklabels([f'{i:.1f}' for i in np.arange(0.5, 2.1, 0.5)])
        cbar.set_label(hue_lab, fontsize=fs, rotation=0, labelpad=10)
    


# plt.subplots_adjust(**adj_kwargs)

cat_labs = [
    'Hyper-\narid', 'Per-\narid', 'Arid', 'Semi-\narid', 
    'Sub-\nhumid', 'Humid', 'Per-\nhumid', 'Hyper-\nhumid'
]
cat_ticks = np.linspace(0, 1, len(ai_labs))
cat_locs = cat_ticks[1:] - np.diff(cat_ticks)[0]/2

# y-axis ticks (MAT)
ax1 = fig.add_axes(ax.get_position(), frameon=False)
ax1.set_xlim((-0, 1.0))
ax1.set_xticks(
    # np.linspace(0, 1, len(ai_labs)),
    cat_locs,
    labels=cat_labs,
    fontsize=fs-1,
)
# ax1.set_xticks(map_divs, labels=map_labs)
# ax1.set_yticks([])
ax1.set_xlabel('')
ax1.set_ylabel('Mean annual temperature (°C)', fontsize=fs)

ax1.set_ylim((0, 1.))
ax1.set_yticks(
    t_divs[1:-1], labels=[f'{i:.3f}'.rstrip('0').rstrip('.') for i in t_labs][1:-1], fontsize=fs,
)
ax1.yaxis.set_label_position('right')
ax1.yaxis.set_label_coords(1.185, 0.34)
ax1.spines['right'].set_position(('outward', 18))
# ax1.tick_params(axis='y', direction='out')
ax1.yaxis.tick_right()
ax1.tick_params(axis='both', length=0,)
ax1.tick_params(axis='y', labelsize=fs)


leg_kwargs = {
    'prop' : {'size': fs},
    'frameon' : False,
    'framealpha' : 0,
    'handletextpad' : 0.25,
    'handlelength' : 1, 
    'fontsize' : fs,
    'loc' : 'upper left',
    'bbox_to_anchor' : (-0.135,1.15),
    'title_fontsize' : fs,
}

leg = utils.add_size_legend(size, s_labs, s_leg, ax=ax1, **leg_kwargs)

# plt.subplots_adjust(**adj_kwargs)
gs.update(**adj_kwargs)





def plot_map(
    df, ax, by, hue_order, cmap, size, size_lims,
    extent, adj_kwargs, leg_kwargs, sleg_kwargs, 
    coast_args={'resolution':'110m', 'linewidth':0.5, 'zorder':2}
):
    # size_lims = (4, 50)
    sizes, s_leg = utils.get_sizes(size, df, size_lims=size_lims, **utils.pt_size_dict[size])
    size_norm = utils.pt_size_dict[size]['size_norm']
    s_labs = utils.pt_size_dict[size]['s_labs']

    labs = [igbp_dict[k]['label'] for k in hue_order] if 'IGBP' in hue else hue_order

    if 'IGBP' in by:
        mask = (df[by].isin(hue_order))
        df = df[mask]

    sns.scatterplot(
        ax=ax,
        x = df.geometry.x,
        y = df.geometry.y,
        hue=df[by], 
        hue_order=hue_order,
        palette=cmap, 
        alpha=0.7,
        size = df[size],
        sizes=size_lims,
        size_norm=size_norm,
        **{'edgecolor' : 'none'},
    )
    if coast_args:
        ax.coastlines(**coast_args)
    # ax.set_global()
    ax.axis('off')
    # Set extent
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    ax.set_yticks([])
    ax.set_ylabel('')


    # Legend
    h = ax.get_legend_handles_labels()
    hp = h[0][1:len(hue_order)+1]

    ax.legend(
        handles=hp, labels=labs,
        **leg_kwargs
    )
    plt.subplots_adjust(**adj_kwargs)

    if 'event' in size:
        s_labs = list(utils.pt_size_dict[size]['s_labs'])
        s_labs[-1] = f'{str(s_labs[-1])} events'

    ax1 = fig.add_axes(ax.get_position(), frameon=False)

    leg = utils.add_size_legend(size, s_labs, s_leg, ax=ax1, **sleg_kwargs)
    ax1.tick_params(axis='both', length=0)
    ax1.set_yticks([])
    ax1.set_xticks([])
    ax1.set_ylabel('')

    return ax

# width = 2.4
# # fig_size = (6.25, 3.25)
# # fig_size = (width,width/2.3)
# # crs = ccrs.Mercator()
# # proj_name = 'Mercator'
# # crs = ccrs.Robinson()
# # proj_name = 'Robinson'
# # crs = ccrs.PlateCarree()
# proj_name = 'PlateCarree'
# crs = utils.map_dict.get(proj_name).get('crs') if not us else ccrs.LambertConformal()

key = 'ismn_t04'
df = filt_dfs[key]


dfg = get_site_data(df, by=site_cols, var='q_q', how=['median','count'], crs=crs)
dfu = get_site_data(df, by=site_cols, var='q_q', how=['median','count'], crs=ccrs.LambertConformal())



extent = utils.map_dict.get(proj_name).get('extent', [-180, 180, -56, 84])

us_extent = utils.map_dict.get(proj_name).get('us_extent', [-121.5, -71.5, 24, 51])

leg_kwargs = {
    'frameon':False,
    'handletextpad': -0.075,
    'labelspacing': 0.4,
    'markerscale': 0.8,
    'fontsize': 9,
    'title_fontsize': 9,
    'loc':'lower left',
    'bbox_to_anchor':(-0.02, -0.02),
}
sleg_kwargs = leg_kwargs.copy()
# sleg_kwargs.update({'loc':'lower right', 'bbox_to_anchor':(-0.02, .32)})
sleg_kwargs.update({'loc':'lower right', 'bbox_to_anchor':(1.17, 0.0),'title_fontsize':fs})

adj_kwargs = {'left':0.005, 'bottom':0.005, 'right':0.995, 'top':0.995, 'wspace':0.3}


# size_lims = (4, 40)
size_lims = (1,15)


if 'IGBP' in hue:
    hue_order, pal = get_palette(hue)
    mask = (dfg[hue].isin(hue_order))

    dfg = dfg[mask]
else:
    hue_order = None
    pal = cmap if hue == 'q_med' else 'summer_r'
    vmin, vmax = (0.25, 2.2) if hue == 'q_med' else (0., 80.)




# mask = (dfg.iIGBP.isin(hue_order))

# dfg = dfg[mask]

# labs = [igbp_dict[k]['label'] for k in hue_order] if 'IGBP' in hue else hue_order


# fig, ax = plt.subplots(1, 1, figsize=fig_size, subplot_kw={'projection': crs})
# ax.add_feature(cartopy.feature.OCEAN, zorder=0, color='lightblue')
# if us:
#     ax.add_feature(cartopy.feature.BORDERS, linewidth=0.5, zorder=1)

axm = plot_map(
    dfg, ax=axm, by=hue, hue_order=hue_order, cmap=pal, size=size, size_lims=size_lims,
    extent=extent, adj_kwargs=adj_kwargs, leg_kwargs=leg_kwargs, sleg_kwargs=sleg_kwargs
)
# for axi in fig.get_axes():
#     axi.legend().remove()

ax_us = plot_map(
    dfu, ax=ax_us, by=hue, hue_order=hue_order, cmap=pal, size=size, size_lims=size_lims,
    extent=us_extent, adj_kwargs=adj_kwargs, leg_kwargs=leg_kwargs, sleg_kwargs=sleg_kwargs,
    coast_args=None
)
shp = shpreader.natural_earth(
    resolution='110m', category='cultural', name='admin_0_countries'
)
# Get us
# countries = [c for c in shpreader.Reader(shp).records()]
# i_us = [i for c in shpreader.Reader(shp).records() if 'United States' in c.attributes['ADMIN']]
i_us = [4]
us = [g for i,g in enumerate(shpreader.Reader(shp).geometries()) if i in i_us]
us = [ccrs.LambertConformal().project_geometry(u, src_crs=ccrs.Geodetic()) for u in us]
# coast = shpreader.natural_earth(
#     resolution='110m', category='physical', name='coastline'
# )
# [g.difference(us) for g in shpreader.Reader(coast).geometries()]

ax_us.add_geometries(
    # shpreader.Reader(shp).geometries()[4],
    us,
    crs = ccrs.LambertConformal(),
    facecolor='none', edgecolor='black', linewidth=0.5
)
# ax_us.add_feature(cartopy.feature.BORDERS, linewidth=0.5, zorder=1,)# resolution='50m')



axes = [ax for ax in fig.get_axes()]

for i in [1, 2, 4]:
    axes[i].get_legend().remove()

for i,lett,loc in zip([3, 4, 5], ['a', 'b', 'c'], [(-0.1, 1.18), (0.08, 1.1), (0.08, 1.06)]):
    axes[i].annotate(
        # f'({lett})',
        f'{lett.capitalize()}', 
        xy=loc, xycoords='axes fraction', fontsize=fs, ha='left', va='top'
    )

plt.subplots_adjust(**adj_kwargs)


# %%

if save:
    name = f'hlz_map_igbp_{hue}_{size}' #+ '_sci'
    plt.savefig(os.path.join(fig_dir, f'{name}.pdf'), dpi=1200, transparent=True)
    plt.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=1200, transparent=True)

# %%
