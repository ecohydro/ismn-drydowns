#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Wed 25 Sep 24 15:06:22'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           surfaces.py
Compatibility:  Python 3.12.0

Description:    Plot surfaces

URL:            https://

Requires:       list of libraries required

Dev ToDo:       None

AUTHOR:         Bryn Morgan
ORGANIZATION:   University of California, Santa Barbara
Contact:        bryn.morgan@geog.ucsb.edu
Copyright:      (c) Bryn Morgan 2024


"""


#%% IMPORTS

import itertools
from scipy import stats
from scipy.optimize import curve_fit
from utils_surfaces import fit_surfaces, calc_midpoints, func, create_dfz, get_traces, calc_r2, contour



#%%


key = 'ismn_t04'
df = filt_dfs[key]
mask = (df.LAI > 0.5) #& (df.perc_nonveg < 100.) #& (df.MOD_IGBP_6s == 'GRA') & (df.AI_60d <= 3)

df = df[mask]
# df = df[df.MOD_IGBP.isin(['GRA','SHR','SAV','WSA','EBF','DBF','ENF','MF', 'BLF', 'OSH','CSH','MOS'])]
# df = df[~df.MOD_IGBP.isin(['URB','WET'])]

# # z-score LAI and AI
# df['z_LAI'] = stats.zscore(df.LAI)
# df['z_AI_60'] = stats.zscore(df.AI_60d, nan_policy='omit')
# # Bin z-scores

by = 'MOD_IGBP_3'
# by = 'perc_woody_bin'
# by = 'iperc_woodyveg_avg_bin'
# by = 'perc_woody_bin2'
x = 'AI_60d_bin'
y = 'LAI_bin'
# x = 'AI_60d'
# y = 'LAI'
by_order, _ = get_palette(by)

if 'bin' in x:
    x_name = re.findall(r'(.*?)_bin', x)[0]
    y_name = re.findall(r'(.*?)_bin', y)[0]
    x_mid, x_edges = calc_midpoints(x, return_edges=True)
    y_mid, y_edges = calc_midpoints(y, return_edges=True)
else:
    x_name = x
    y_name = y

min_events = 5


# x_mid = df.groupby(by=x)[x_name].median().to_list()
# y_mid = df.groupby(by=y)[y_name].median().to_list()

dfr = fit_surfaces(df, x, y, by, func, min_events=min_events, use_meds=False)

# x1 = np.linspace(x_edges[0], x_edges[-1], 100)
# y1 = np.linspace(y_edges[1], y_edges[-1], 100)
# if not 'AIi' in x:
#     x1 = np.linspace(x_edges[0], x_mid[-1], 100)
# else:
x_vals = x_edges
x_vals[-1] = x_mid[-1]

# x1 = np.unique(np.concatenate([np.linspace(x_vals[i], x_vals[i+1], 13) for i in range(len(x_vals)-1)]))

x1 = np.linspace(x_edges[0], x_mid[-1], 100)
y1 = np.linspace(y_edges[1], y_mid[-1], 100)
# x1 = x_edges
# y1 = y_edges[1:]

dfz = create_dfz(dfr, by, x1, y1, x_name=x_name, y_name=y_name, z_name='q')

# dfz.rename(columns={by: by[1:]}, inplace=True)
# dfr.rename(columns={by: by[1:]}, inplace=True)

#%% SURFACES

save = False

# x = 'AI_60d_bin2'
# y = 'LAI_bin'
# x = 'AIi_60d'
# y = 'LAI'
# x = x_name
z = 'q'


# c_list = ['#1B3A39', '#5BB8B0', '#82B960', '#B9C354', '#EDC01D']
# c_listr = ['#EDC01D', '#82B960', '#5BB8B0', '#1B3A39']
# cmap = get_continuous_cmap(c_listr)
# c_list1 = ['#026f6f', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
# c_list1 = ['#025951', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
c_list1 = ['#004446','#026f6f', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
cmap = get_continuous_cmap(c_list1[::-1])

# cmap = 'husl'
n = len(dfr)
pal1 = list(sns.color_palette('summer_r', n_colors=n+2).as_hex())[:n]

aggfunc = 'median'


inv_x = False if 'AIi' in x else True

fs = 9
# fig_size = (n*2.4, 2.5)
# rs, cs = 1, n
vmin, vmax = 0.25, 2.2
# cax_dims = [0.95, 0.2, 0.01, 0.68]
# adj_kwargs = {'left':0.031, 'right':0.94, 'bottom':0.145, 'top':0.925, 'wspace':0.105}

# fig_size, r, c, cax_dims, adj_kwargs = get_fig_dets(n)
cm = 1/2.54
fig_size, r, c, cax_dims, adj_kwargs = get_fig_dets(n, width=17.8*cm)
# fig_size = (fig_size[0], fig_size[1]*1.1)
adj_kwargs['top'] = 0.952
adj_kwargs['bottom'] = 0.098
adj_kwargs['left'] = 0.062

# xticks = get_ticks(dfz, x_name, len(x1))
# xticks[0][0] = 0.5
xlabs = [0, 0.5, 1.0, 1.5, 2.0, 2.5]
xticks = (utils.rescale(np.array(xlabs), xmin=0, xmax=np.max(x1), new_min=0, new_max=len(x1)), xlabs)
xticks[0][0] = 0.75
# xticks[1] = [f'{i:.2f}'.rstrip('0').rstrip('.') for i in xlabs]
ylabs = np.arange(1, 6, 1.)
yticks = (utils.rescale(ylabs, xmin=np.min(y1), xmax=np.max(y1), new_min=0, new_max=len(y1)), ylabs)

if r == 1:
    adj_kwargs['top'] = 0.912
    adj_kwargs['bottom'] = 0.18
    adj_kwargs['left'] = 0.064

adj_kwargs = {'left': 0.059, 'bottom': 0.15, 'top':0.922, 'right':0.9, 'wspace':0.1}
cax_dims = [0.925, 0.15, 0.02, 0.7715]

fig, axs = plt.subplots(r, c, figsize=fig_size, sharex=True, sharey=True)
# cax_dims = [0.915, 0.18, 0.02, 0.73]    # dissertation
cax = fig.add_axes(cax_dims)


axs, cax = plot_heatmaps(
    dfz, axs, x_name, y_name, 'q', by, cmap, cax, vmin, vmax, fs=9, 
    adj_kwargs=adj_kwargs, annot_kwargs={'annot' : False},
    contours=True, clab = f'{aggfunc.capitalize()} $q$', invert_x=inv_x,
    contour_kwargs={'levels' : [0.6, 0.8, 1., 1.2, 1.4, 1.6, 2., 2.5]}, contour_label=True,
    ticks=True, xticks=xticks, yticks=yticks,
    ylab='LAI (m$^2$ m$^{-2}$)', xlab='AI (preceding 60 days)'
)

# x10 = np.linspace(x_edges[0], x_mid[-1], 1000)
# y10 = np.linspace(y_edges[1], 7, 1000)

# traces = get_traces(
#     create_dfz(
#         dfr, by, x=x10, y=y10, x_name=x_name, y_name=y_name, z_name=z, derivatives=False
#     ), 
#     dfr, by, x_name, y_name, z
# )

# for i,(wc,ax) in enumerate(zip(by_order, axs.flatten())):
#     lab = igbp_dict.get(wc, {}).get('label', wc)
#     lab = re.match('^\S*',lab)[0]
#     trc = traces[traces[by] == wc]
#     trc_sc = rescale(trc['trace'], new_min=0, new_max=np.max(ax.get_ylim()))
#     x_sc = rescale(trc[x_name], new_min=0, new_max=np.max(ax.get_xlim()))
#     line = ax.plot(x_sc, trc_sc, label=wc, color='k')
#     # label_line(line[0], lab, near_x=xc_locs[i], rotation_offset=-180)

for i,(ax,lett) in enumerate(zip(axs.flatten(),['d','e','f'])):
    title = f'{ax.get_title()} ($R^2 $ = {np.round(dfr.r2.iloc[i],2)})'
    ax.text(0.0, 1.09, f'({lett})', transform=ax.transAxes, fontsize=fs, va='top')

    if lett == 'f':
        title = '  ' + title
    ax.set_title(
        title, fontsize=fs, pad=5
    )
plt.subplots_adjust(**adj_kwargs)

if save:
    name = f'predicted_q_{key[5:]}_{by}_{x_name}_{y_name}' + '_lett'
    plt.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=1200, transparent=True)
    plt.savefig(os.path.join(fig_dir, f'{name}.pdf'), transparent=True)

#%% HEATMAPS + SURFS

save = False


key = 'ismn_t04'
df = filt_dfs[key]

col = 'MOD_IGBP_3'

aggfunc = 'median'
# aggfunc = lambda x : x.var()/x.mean()

y = 'LAI_bin'
x = 'AI_60d_bin'
# x = 'AI_rel_bin'
# xlab = "AI$_{60}$"
xlab = 'AI (preceding 60 days)'


# x = 'AI_bin'
# # xlab = r"AI ($\frac{P}{PET}$)"
# xlab = r"Aridity index"
df = df[df.MOD_IGBP.isin(['GRA','SHR','SAV','WSA','EBF','DBF','ENF','MF', 'BLF', 'OSH','CSH','MOS'])]
# df = df[df.perc_nonveg < 50.]

val = 'q_q'

if 'IGBP' in col:
    igbp_order = igbp_cats.get(
        re.findall(r"(?<=_)[^_]*(?!.*_)", col)[0] if '_' in col else col, igbp_cats['og']
    )
else:
    colk = col[1:] if col[0] == 'i' else col
    colk = colk[:-8]+colk[-4:] if colk[-8:-4] == '_avg' else colk
    igbp_order = bin_dict.get(re.findall(r'.*?bin', colk)[0])[f'{re.findall(r'bin.*$', colk)[0]}_labels']

n = len(igbp_order) * 2
cm = 1/2.54
fig_size, r, c, cax_dims, adj_kwargs = get_fig_dets(n, width=17.8*cm)
fig_size = (7, 4.9)
# fig_size, r, c, cax_dims, adj_kwargs = get_fig_dets(n)
# adj_kwargs = {'left':0.053, 'right':0.9025, 'bottom':0.135, 'top':0.922, 'wspace':0.1}
cax_dims = [0.9265, 0.075, 0.02, 0.885]


c_list1 = ['#004446','#026f6f', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
# c_list1 = ['#025951', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
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
# fig, axs = plt.subplots(1, len(igbp_order), figsize=(2.4*len(igbp_order), 2.5), sharex=True, sharey=True)
# fig, axs = plt.subplots(2, 3, figsize=(6.25, 4), sharex=True, sharey=True)
fig, axs = plt.subplots(r, c, figsize=fig_size, )#sharex=True, sharey=True)

axs_heat = axs[0]
axs_surf = axs[1]
# cbar_ax = fig.add_axes([0.95, 0.2, 0.01, 0.68])  # Adjust the position and size as needed
# cbar_ax = fig.add_axes([0.915, 0.08, 0.023, 0.88])  # Adjust the position and size as needed
# cax_dims = [0.915, 0.155, 0.02, 0.765]
cax = fig.add_axes(cax_dims)  # Adjust the position and size as needed

# cax = divider.append_axes("right", size="5%", pad=0.15)
# plt.colorbar(im, cmap=cmap, cax=cax, orientation='vertical', ticks=cticks[0], label=clab)
# # cb = plt.colorbar(im, cax=cax, )
# # cb.outline.set_visible(False)

axs_heat, cax = plot_heatmaps(
    df, axs[0], x, y, val, col, cmap, cax, vmin=vmin, vmax=vmax, cticks=cticks, aggfunc=aggfunc, fs=fs,
    clab=clab, xlab=xlab, ylab='LAI (m$^2$ m$^{-2}$)', mask_count=5, invert_x=inv_x
)

for ax,lett in zip(axs_heat, ['a','b','c','d','e','f']):
    # ax.text(0.0, 1.09, f'({lett})', transform=ax.transAxes, fontsize=fs, va='top')
    ax.text(0.0, 1.0875, f'{lett.capitalize()}', transform=ax.transAxes, fontsize=fs, va='top')

plt.subplots_adjust(
    **adj_kwargs
)

z = 'q'

c_list1 = ['#004446','#026f6f', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
cmap = get_continuous_cmap(c_list1[::-1])

# cmap = 'husl'
n = len(dfr)
pal1 = list(sns.color_palette('summer_r', n_colors=n+2).as_hex())[:n]

aggfunc = 'median'


inv_x = False if 'AIi' in x else True

fs = 9
# fig_size = (n*2.4, 2.5)
# rs, cs = 1, n
vmin, vmax = 0.25, 2.2
# cax_dims = [0.95, 0.2, 0.01, 0.68]
# adj_kwargs = {'left':0.031, 'right':0.94, 'bottom':0.145, 'top':0.925, 'wspace':0.105}

# fig_size, r, c, cax_dims, adj_kwargs = get_fig_dets(n)
cm = 1/2.54
fig_size, r, c, cax_dims, adj_kwargs = get_fig_dets(n, width=17.8*cm)
# fig_size = (fig_size[0], fig_size[1]*1.1)
adj_kwargs['top'] = 0.952
adj_kwargs['bottom'] = 0.098
adj_kwargs['left'] = 0.062

# xticks = get_ticks(dfz, x_name, len(x1))
# xticks[0][0] = 0.5
xlabs = [0, 0.5, 1.0, 1.5, 2.0, 2.5]
xticks = (utils.rescale(np.array(xlabs), xmin=0, xmax=np.max(x1), new_min=0, new_max=len(x1)), xlabs)
xticks[0][0] = 0.75
# xticks[1] = [f'{i:.2f}'.rstrip('0').rstrip('.') for i in xlabs]
ylabs = np.arange(1, 6, 1.)
yticks = (utils.rescale(ylabs, xmin=np.min(y1), xmax=np.max(y1), new_min=0, new_max=len(y1)), ylabs)

if r == 1:
    adj_kwargs['top'] = 0.912
    adj_kwargs['bottom'] = 0.18
    adj_kwargs['left'] = 0.064

# adj_kwargs = {'left': 0.059, 'bottom': 0.15, 'top':0.922, 'right':0.9, 'wspace':0.1}
adj_kwargs = {'left': 0.059, 'bottom': 0.075, 'top':0.962, 'right':0.9, 'wspace':0.1, 'hspace':0.375}
# cax_dims = [0.925, 0.5, 0.02, 0.7715]

# fig, axs = plt.subplots(r, c, figsize=fig_size, sharex=True, sharey=True)
# # cax_dims = [0.915, 0.18, 0.02, 0.73]    # dissertation
# cax = fig.add_axes(cax_dims)


axs_surfs, cax = plot_heatmaps(
    dfz, axs[1], x_name, y_name, 'q', by, cmap, cax, vmin, vmax, fs=9, 
    adj_kwargs=adj_kwargs, annot_kwargs={'annot' : False},
    contours=True, clab = f'{aggfunc.capitalize()} $q$', invert_x=inv_x,
    contour_kwargs={'levels' : [0.6, 0.8, 1., 1.2, 1.4, 1.6, 2., 2.5]}, contour_label=True,
    ticks=True, xticks=xticks, yticks=yticks,
    ylab='LAI (m$^2$ m$^{-2}$)', xlab='AI (preceding 60 days)'
)
for ax in axs_surfs[1:]:
    ax.set_yticks(
        yticks[0], labels=[], rotation=0, ha='right', rotation_mode='anchor', fontsize=fs-0.5
    )
    


for i,(ax,lett) in enumerate(zip(axs_surfs.flatten(),['d','e','f'])):
    title = f'{ax.get_title()} ($R^2 $ = {np.round(dfr.r2.iloc[i],2)})'
    if i == 2:
        title = ' ' + title
    # ax.text(0.0, 1.09, f'({lett})', transform=ax.transAxes, fontsize=fs, va='top')
    ax.text(0.0, 1.0875, f'{lett.capitalize()}', transform=ax.transAxes, fontsize=fs, va='top')

    if lett == 'f':
        title = '  ' + title
    ax.set_title(
        title, fontsize=fs, pad=5
    )
plt.subplots_adjust(**adj_kwargs)

if save:
    name = f'heatmaps_combined' + '_sci' #f'predicted_q_{key[5:]}_{by}_{x_name}_{y_name}' + '_lett'
    plt.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=1200, transparent=True)
    plt.savefig(os.path.join(fig_dir, f'{name}.pdf'), transparent=True)


#%% GRADIENTS
 
save = False

fig, axs = plt.subplots(r, c, figsize=fig_size, sharex=True, sharey=True)
cax = fig.add_axes(cax_dims)


axs, cax = plot_heatmaps(
    dfz, axs, x_name, y_name, 'slp', by, cmap, cax, vmin=0, vmax=1.2, fs=9, 
    adj_kwargs=adj_kwargs, annot_kwargs={'annot' : False},
    contours=True, cticks=(np.arange(0., 1.5, 0.5), None),
    clab = '$\\nabla q ($ AI$_{60}$, LAI$)$', invert_x=inv_x
)
plt.subplots_adjust(**adj_kwargs)


if save:
    name = f'gradient_{val[2:]}_{key[5:]}_{by}_{x_name}_{y_name}'
    plt.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=600, transparent=True)
    plt.savefig(os.path.join(fig_dir, f'{name}.pdf'), transparent=True)


# Plot dz/dx and dz/dy
fig, axs = plt.subplots(r, c, figsize=fig_size, sharex=True, sharey=True)
cax = fig.add_axes(cax_dims)

axs, cax = plot_heatmaps(
    dfz, axs, x_name, y_name, f'dq_d{x_name}', by, cmap, cax, vmin=-1, vmax=1, fs=9, 
    adj_kwargs=adj_kwargs, annot_kwargs={'annot' : False},
    contours=True, cticks=(np.arange(-1., 1.5, 0.5), None),
    clab = '$\\partial q / \\partial $ AI$_{60}$', invert_x=inv_x
)
fig, axs = plt.subplots(r, c, figsize=fig_size, sharex=True, sharey=True)
cax = fig.add_axes(cax_dims)

axs, cax = plot_heatmaps(
    dfz, axs, x_name, y_name, f'dq_d{y_name}', by, cmap, cax, vmin=-1, vmax=1, fs=9, 
    adj_kwargs=adj_kwargs, annot_kwargs={'annot' : False},
    contours=True, cticks=(np.arange(-1., 1.5, 0.5), None),
    clab = '$\\partial q / \\partial $ AI$_{60}$', invert_x=inv_x
)

