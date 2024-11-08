#-------------------------------------------------------------------------------
# FIG. 1. CONCEPTUAL FIGURE
#-------------------------------------------------------------------------------

# IMPORTS 
import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# %matplotlib qt

from drydown import calc_dtheta_dt, calc_theta_q
from utils_figs import get_continuous_cmap


# CONCEPTUAL FIGURE
save = False
# adj_kwargs = {'left': 0.163, 'bottom': 0.165, 'top':0.985, 'wspace':0.3}
# adj_kwargs = {'left':0.06, 'right':0.965, 'bottom':0.19, 'top':0.947, 'wspace':0.38, 'hspace':0.} # dissertation
adj_kwargs = {'left':0.055, 'right':0.97, 'bottom':0.19, 'top':0.947, 'wspace':0.38, 'hspace':0.}


ratio = 6.25/2
cm = 1/2.54

width = 17.8*cm
fig_size = (width, width/ratio)

# fig_size = (6.25, 2.)
# fig_size = (6.25, 2.83)
# Set up grid spec
fig = plt.figure(figsize=fig_size)
# gs = fig.add_gridspec(4, 3, width_ratios=[0.38, 0.3, 0.03], height_ratios=[0.09,0.4, 0.58,0.09], **adj_kwargs)
gs = fig.add_gridspec(3, 4, width_ratios=[1,1,0,1.08], height_ratios=[0.04,1.04,0.08], **adj_kwargs)
#   width_ratios=[0.57, 0.43], height_ratios=[0.58, 0.42], **adj_kwargs)



ax1 = fig.add_subplot(gs[1, 0], )
ax2 = fig.add_subplot(gs[1, 1])
# ax1.set_aspect('box')
# ax2.set_aspect('box')
# gs2 = gs[0, 3].subgridspec(1, 2, width_ratios=[1,0.05], hspace=0.1)

# ax = fig.add_subplot(gs2[0, 0])
# cax = fig.add_subplot(gs2[0, 1])
ax = fig.add_subplot(gs[1, 3])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)

# aax = fig.add_subplot(gs[1, 3])

theta_fc = 0.4
theta_0 = 0.4
theta_star = 0.36
theta_w = 0.125
# tau = 1.8
et_max = 2.5
k = et_max/(0.05*1000)

t_star = (theta_0 - theta_star)/k

t1 = np.linspace(0, t_star, 25)
t2 = np.linspace(0, 16, 261)
t = np.concatenate([t1, np.linspace(t_star, 16, len(t2))])

qpos = 1.5
qneg = 0.75



theta = np.linspace(theta_fc, theta_w-0.03, len(t))
# theta_exp = calc_theta_exp(t, theta_0, theta_w, tau, theta_star)
theta_exp = calc_theta_q(t2, theta_star, k, 1.0000000001, theta_star, theta_w)
theta_pos = calc_theta_q(t2, theta_star, k, qpos, theta_star, theta_w)
theta_neg = calc_theta_q(t2, theta_star, k, qneg, theta_star, theta_w)

np.concatenate([-k * t[t<=1] + theta_0, theta_exp])

df = pd.DataFrame({
    't' : t,
    'theta_exp' : np.concatenate([-k * t1 + theta_0, theta_exp]),
    'theta_pos' : np.concatenate([-k * t1 + theta_0, theta_pos]),
    'theta_neg' : np.concatenate([-k * t1 + theta_0, theta_neg]),
    'dtheta_exp' : calc_dtheta_dt(
        theta = theta, k = k, q = 1.0, theta_w = theta_w, theta_star = theta_star
    ),
    'dtheta_pos' : calc_dtheta_dt(
        theta = theta, k = k, q = qpos, theta_w = theta_w, theta_star = theta_star
    ),
    'dtheta_neg' : calc_dtheta_dt(
        theta = theta, k = k, q = qneg, theta_w = theta_w, theta_star = theta_star
    ),
    'dtheta_mm_exp' : calc_dtheta_dt(
        theta = theta, k = et_max, q = 1.0, theta_w = theta_w, theta_star = theta_star
    ),
    'dtheta_mm_pos' : calc_dtheta_dt(
        theta = theta, k = et_max, q = qpos, theta_w = theta_w, theta_star = theta_star
    ),
    'dtheta_mm_neg' : calc_dtheta_dt(
        theta = theta, k = et_max, q = qneg, theta_w = theta_w, theta_star = theta_star
    ),
})

fs = 9.



# fig = plt.figure(figsize=(5.25, 2.4))
# ax1 = fig.add_subplot(1,2,1)
# ax2 = fig.add_subplot(1,2,2)

lw = 0.5

ax1.hlines(
    y=[theta_star,theta_w], xmin=0, xmax=[(theta_0-theta_star)/k, np.max(t)], color='k', linestyle='--', linewidth=lw
)
ax2.hlines(
    y=et_max, xmin=0, xmax=1, color='k', linestyle='--', linewidth=lw
)
ax2.vlines(
    x=[theta_star,theta_star], ymin=[-1,et_max*0.35-0.15], ymax=[et_max*0.35-0.55,et_max], 
    color='k', linestyle='--', linewidth=lw
)


# colors = ['#82B960', '#5BB8B0', '#EDC01D']
colors = ['#89bc42', '#12bbaf', '#f4c807']
anns = ['$q = 1$','$q > 1$\nconservative','aggressive\n$q < 1$']
locs = [
    (theta_w + (theta_star-theta_w)*0.5, et_max*0.9), 
    (theta_star-(theta_star-theta_w)*0.2, et_max*0.35),
    (theta_w+(theta_star-theta_w)*0.2, et_max*0.65), 
]
qs = [1.0000000001, qpos, qneg]

for i,(suff,lab) in enumerate(zip(['exp','pos','neg'], anns)):
    ax1.plot(
        df.t, df[f'theta_{suff}'], label=lab, color=colors[i]
    )
    ax2.plot(
        theta, df[f'dtheta_mm_{suff}'], label=lab, color=colors[i]
    )
    if i == 1:
        ymin = locs[i][1] + 0.08
        ymax = calc_dtheta_dt(locs[i][0], et_max, qs[i], theta_w, theta_star)
    else:
        ymax = locs[i][1] - 0.08
        ymin = calc_dtheta_dt(locs[i][0], et_max, qs[i], theta_w, theta_star)

    ax2.vlines(
        x=locs[i][0], ymin=ymin, ymax=ymax, color=colors[i], linestyle='--', linewidth=1
    )
    tex = ax2.text(
        x=locs[i][0], y=locs[i][1], s=lab, color=colors[i], fontsize=fs-1, ha='center',
        va='bottom' if not i == 1 else 'top',
    )
    print(tex)

ax1.set_xlabel('$t$ (d)', fontsize=fs)
ax1.set_ylabel(r'$\theta$ (m$^3$ m$^{-3}$)', fontsize=fs, labelpad=-1)
ax2.set_xlabel(r'$\theta$ (m$^3$ m$^{-3}$)', fontsize=fs, labelpad=-1)
ax2.set_ylabel(r'$-\Delta z \, \, \frac{d\theta}{dt}$ (mm d$^{-1}$)', fontsize=fs, labelpad=-8)
# ax2.legend()
ax1.set_xlim(0, np.max(t))
ax1.set_ylim(theta_w-0.025, theta_fc)
ax2.set_xlim(np.min(theta), np.max(theta))
ax2.set_ylim(-0.1, et_max+0.25)

ax1.set_xticks([])
ax1.set_yticks(
    [theta_w, theta_star, theta_fc], 
    labels=[r'$\theta_{\mathrm{w}}$', r'$\theta_{\ast}$', r'$\theta_{\mathrm{fc}}$'],
    fontsize=fs
)
ax2.set_xticks(
    [theta_w, theta_star, theta_fc], 
    labels=[r'$\theta_{\mathrm{w}}$', r'$\theta_{\ast}$', r'$\theta_{\mathrm{fc}}$'],
    fontsize=fs
)
ax2.set_yticks(
    [0, et_max], 
    labels=['0 ', r'$\mathrm{ET}_{\mathrm{max}}$'],
    # pad=-2,
    # fontsize=fs
)
ax2.tick_params(axis='y', pad=-0.5)
for tick,size in zip(ax2.get_yticklabels(), [fs,fs-1]):
    tick.set_fontsize(size)




# adj_kwargs = {'left':0.065, 'right':0.99, 'bottom':0.195, 'top':0.972, 'wspace':0.3, 'hspace':0.3}

# HEATMAP

# fig_size = (3.1, 2.4)
# adj_kwargs = {'left': 0.163, 'bottom': 0.165, 'top':0.985, 'wspace':0.3}
fs = 9


size = 100

df = pd.DataFrame([[i + j for j in range(size)] for i in range(size)])

x1 = np.linspace(0, 1., 100)
y1 = np.linspace(0, 1., 100)
X, Y = np.meshgrid(x1, y1)
# Z = func((X,Y), 1., 1., -.5, 0.25)
Z = np.sqrt((X-1)**2 + Y**2)

# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection='3d')
# ax.plot_surface(X, Y, Z, cmap='viridis')

c_list = ['#1B3A39', '#5BB8B0', '#82B960', '#B9C354', '#EDC01D']
c_listr = ['#EDC01D', '#82B960', '#5BB8B0', '#1B3A39']
cmap = get_continuous_cmap(c_listr[::-1])
# c_list1 = ['#026f6f', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
c_list1 = ['#025951', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
# c_list1 = ['#004446','#026f6f', '#12bbaf', '#65bb5e', '#adbd25', '#f4c807']
cmap = get_continuous_cmap(c_list1)

# vmin, vmax = (0, 200)
# vmin,vmax = 0.2,2.2
vmin, vmax = np.min(Z), np.max(Z)

# fig, ax = plt.subplots(1, 1, figsize=fig_size)
# ax = ax3
ax.set_aspect('equal')

im = ax.imshow(
    # df, 
    Z,
    cmap=cmap,
    vmin=vmin, vmax=vmax
)



ax.spines[['top','bottom','left','right']].set_visible(False)

# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.15)
cb = plt.colorbar(im, cax=cax, )
cb.outline.set_visible(False)

# cax.set_yticks([0, 100, 200], labels=[r"$q < 1$", r"$q = 1$", r"$q > 1$"])
# cax.set_yticks([2.5, 100, 197.5], labels=[r"<1", r"1", r">1"], fontsize=fs)
cax.set_yticks(
    [0.025*np.max(Z), 0.5*np.max(Z), np.max(Z)*0.975], 
    # labels=[r"<1", r"1", r">1"], 
    labels=[r">1", r"1", r"<1"],
    fontsize=fs
)
cax.tick_params(axis='y', length=0)
cax.set_ylabel(r"$q$", fontsize=fs, rotation=0, labelpad=0)# labelpad=10)
ax.invert_yaxis()
ax.invert_xaxis()
cax.invert_yaxis()

ax.xaxis.set_ticks([5, 95], labels=[r'Arid', r'Humid'], fontsize=fs-1, rotation=0)
ax.yaxis.set_ticks([5, 95], labels=[r'Low', r'High'], fontsize=fs-1, rotation=0)

ax.tick_params(axis='both', length=0, labelsize=fs-1)

ax.set_xlabel('Aridity index', fontsize=fs, labelpad=-2)
# ax.set_ylabel('Leaf area index (m$^2$ m$^{-2}$)', fontsize=fs, labelpad=0)
ax.set_ylabel('LAI (m$^2$ m$^{-2}$)', fontsize=fs, labelpad=-7)

ax.annotate(
    '', xy=(0, -0.225), xycoords='axes fraction', xytext=(1, -0.225), 
    arrowprops=dict(arrowstyle='<|-', color='k', lw=0.5)
)
ax.annotate(
    '', xy=(-0.32, 0), xycoords='axes fraction', xytext=(-0.32, 1), 
    arrowprops=dict(arrowstyle='<|-', color='k', lw=0.5)
)
ax.annotate(
    'Hydroclimatic demand', xy=(0.5, -0.335), xycoords='axes fraction',
    ha='center', fontsize=fs
)
ax.annotate(
    'Ecological demand', xy=(-0.38, 0.5), xycoords='axes fraction',
    ha='center', va='center', rotation=90, fontsize=fs
)
ax.annotate(
    'conservative', xy=(0.075, 0.15), xycoords='axes fraction',
    ha='left', va='center', fontsize=fs-1, color='w'
)
ax.annotate(
    'aggressive', xy=(1-0.075, 1-0.15), xycoords='axes fraction',
    ha='right', va='center', fontsize=fs-1, color='w'
)

for axi,lett in zip([ax1,ax2,ax],['a','b','c']):
    axi.annotate(
        # f'({lett})',
        f'{lett.capitalize()}', 
        xy=(0.0, 1.05), xycoords='axes fraction', fontsize=fs
    )

# plt.subplots_adjust(**adj_kwargs)

# %%

if save:
    name = 'conceptual' #+ '_sci'
    plt.savefig(os.path.join(fig_dir, f'{name}.pdf'), dpi=1200, transparent=True,) # bbox_inches='tight'
    plt.savefig(os.path.join(fig_dir, f'{name}.png'), dpi=1200, transparent=True,) # bbox_inches='tight'

# %%
