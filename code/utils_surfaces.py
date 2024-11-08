#!usr/bin/env python
# -*- coding: utf-8 -*-
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

__author__ = 'Bryn Morgan'
__contact__ = 'bryn.morgan@geog.ucsb.edu'
__copyright__ = '(c) Bryn Morgan 2024'

__license__ = 'MIT'
__date__ = 'Tue 14 May 24 20:00:02'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''

"""

Name:           gradient.py
Compatibility:  Python 3.12.0

Description:    Utility functions for plotting surfaces/heatmaps.

URL:            https://

Requires:       list of libraries required

Dev ToDo:       None

AUTHOR:         Bryn Morgan
ORGANIZATION:   University of California, Santa Barbara
Contact:        bryn.morgan@geog.ucsb.edu
Copyright:      (c) Bryn Morgan 2024


"""

# IMPORT 
import re
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import itertools

from utils_figs import bin_dict, get_palette

abc = [chr(i) for i in range(ord('a'),ord('z')+1)]

# FUNCTIONS FOR SURFACE

def func(xy, a, b, c, d):
    x, y = xy
    return a*x + b*y + c*x*y + d

def saddle(x, y, a, b, c, d):
    return -b/c, -a/c, func((-b/c, -a/c), a, b, c, d)


def dz_dx(x, y, a, b, c, d):
    return a + c*y

def dz_dy(x, y, a, b, c, d):
    return b + c*x

def calc_grad(
        dzdx=None, dzdy=None, x=None, y=None, a=None, b=None, c=None, d=None
    ):
    if dzdx is None:
        dzdx = dz_dx(x, y, a, b, c, d)
    if dzdy is None:
        dzdy = dz_dy(x, y, a, b, c, d)
    return np.sqrt((dzdx)**2 + (dzdy)**2)

# def grad2(x, y, a, b, c, d):
#     grad = calc_grad(x, y, a, b, c, d)
#     dgdx = (c * (c * x + b)) / grad
#     dgdy = (c * (c * y + a)) / grad
#     return calc_grad(dzdx=dgdx, dzdy=dgdy)

def contour(x, z, a, b, c, d):
    return -(a*x + d - z) / (b + c * x)

def hessian(x, y, a, b, c, d):
    # dzdx = dz_dx(x, y, a, b, c, d)
    # dzdy = dz_dy(x, y, a, b, c, d)
    return np.array([[0, c], [c, 0]])


def get_steepest(x, y, coeffs, func):
    # Calculate the Hessian
    h = hessian(x, y, **coeffs)
    # Calculate the eigenvalues and eigenvectors
    eigvals, eigvects = np.linalg.eig(h)
    # Create the path of steepest descent
    t = np.linspace(-10, 10, 100)
    x1 = eigvects[0,0] * t
    y1 = eigvects[1,0] * t
    s1 = func((x1, y1), **coeffs)

    x2 = eigvects[0,1] * t
    y2 = eigvects[1,1] * t
    s2 = func((x2, y2), **coeffs)

    return (x1, y1, s1), (x2, y2, s2)


# GENERAL FUNCTIONS

def calc_midpoints(col, return_edges=False):
    key = f'edges{col[-1]}' if col[-1] != 'n' else 'edges'
    edges = bin_dict.get(re.findall(r'.*?bin', col)[0])[key].copy()
    midpoints = np.round([(edges[i] + edges[i+1]) / 2 for i in range(len(edges) - 1)], 5)
    if return_edges:
        return midpoints, edges
    return midpoints


def calc_r2(x, y, popt, func):
    y_opt = func(x, *popt)
    residuals = y - y_opt
    ss_res = np.sum(residuals**2)
    
    return 1 - ss_res / np.sum((y - np.mean(y))**2)

def fit_surface(df, x, y, by, func, z='q_q'):

    mask = df[[x, y, z]].isna().any(axis=1)

    xdata = (df[~mask][x], df[~mask][y])
    ydata = df[~mask][z]

    popt, pcov = curve_fit(func, xdata, ydata)

    r2 = calc_r2(xdata, ydata, popt, func)
    # y_opt = func(xdata, *popt)

    # # Calculate residuals
    # residuals = ydata - y_opt
    # ss_res = np.sum(residuals**2)
    
    # r2 = 1 - ss_res / np.sum((ydata - np.mean(ydata))**2)
    results = {
        by: df[by].iloc[0],
        # 'a': popt[0],
        # 'b': popt[1],
        # 'c': popt[2],
        # 'd': popt[3],
        **{abc[i]: popt[i] for i in range(len(popt))},
        'r2': r2
    }
    return results


def fit_surfaces(df, x, y, by, func, z='q_q', min_events=5, use_meds=False):

    if 'bin' in x:
        x_name = re.findall(r'(.*?)_bin', x)[0]
        y_name = re.findall(r'(.*?)_bin', y)[0]
        x_mid, x_edges = calc_midpoints(x, return_edges=True)
        y_mid, y_edges = calc_midpoints(y, return_edges=True)
    else:
        x_name = x
        y_name = y

    if use_meds:
        x_mid = df.groupby(by=x)[x_name].median().to_list()
        y_mid = df.groupby(by=y)[y_name].median().to_list()

    by_order, pal = get_palette(by)

    if 'bin' in x:
        dfg = df.groupby(
            by=[x, y, by]
        )[z].median().reset_index() 
        dfg = dfg.merge(
            df.groupby(
                by=[x, y, by])['q_q'].count().reset_index().rename(columns={z:'n_events'})
        )
        dfg.loc[dfg.n_events < min_events, z] = np.nan
        dfg[f"{x}_mid"] = dfg[x].replace(dict(zip(dfg[x].unique(), x_mid))).astype(float)
        dfg[f"{y}_mid"] = dfg[y].replace(dict(zip(dfg[y].unique(), y_mid))).astype(float)
    else:
        dfg = df
        
    x_col = f"{x}_mid" if 'bin' in x else x
    y_col = f"{y}_mid" if 'bin' in y else y
    

    results = [
        fit_surface(dfg[dfg[by] == cat], x_col, y_col, by, func, z=z) for cat in by_order
    ] 
    return pd.DataFrame(results)


def calc_derivatives(df, coeffs : dict, x_name='x', y_name='y', z_name='z'):

    # df = pd.DataFrame(list(itertools.product(x, y)), columns=[x_name, y_name])
    # df[z_name] = func((df[x_name], df[y_name]), **coeffs)
    df[f'd{z_name}_d{x_name}'] = dz_dx(df[x_name], df[y_name], **coeffs)
    df[f'd{z_name}_d{y_name}'] = dz_dy(df[x_name], df[y_name], **coeffs)
    
    # Calculate gradient
    # df['slp'] = np.sqrt(df[f'd{z_name}_d{x_name}']**2 + df[f'd{z_name}_d{y_name}']**2)
    df['slp'] = calc_grad(x=df[x_name], y=df[y_name], **coeffs)
    dz_dx_max = df[f'd{z_name}_d{x_name}'].max()
    dz_dy_max = df[f'd{z_name}_d{y_name}'].max()
    # Normalize gradient
    df['slp_norm'] = np.sqrt(
        (df[f'd{z_name}_d{x_name}'] / dz_dx_max)**2 + (df[f'd{z_name}_d{y_name}'] / dz_dy_max)**2
    )

    # Calc derivative of gradient
    # df['grad2'] = grad2(df[x_name], df[y_name], **coeffs)

    return df

def create_dfz(dfr, by, x, y, x_name='x', y_name='y', z_name='z', derivatives=True):

    dfs = []
    for i,row in dfr.iterrows():
        coeffs = row[['a', 'b', 'c', 'd']].to_dict()

        df = pd.DataFrame(list(itertools.product(x, y)), columns=[x_name, y_name])
        df[z_name] = func((df[x_name], df[y_name]), **coeffs)
        
        if derivatives:
            df = calc_derivatives(df, coeffs, x_name, y_name, z_name)

        df.insert(0, by, row[by])
        dfs.append(df)
    
    return pd.concat(dfs, ignore_index=True)

# GRADIENT DESCENT INTERSECTING THE SADDLE POINT

# Steps
# - Get saddle point
# - Check if saddle point is in domain.
# - If not, find global min.
# - If saddle point is in domain, find local minima (quick and dirty)


def get_minima(data):
    # Get global minimum
    z_min = data.min().min()
    i_min = np.where(data == z_min)
    x_min, y_min = data.columns[i_min[1][0]], data.index[i_min[0][0]]
    # Check index
    ix = -1 if i_min[0][0] == 0 else 0
    iy = -1 if i_min[1][0] == 0 else 0

    # Get local minimum in top row or left column (?)
    z_top = np.min([data.iloc[iy].min(),data.iloc[:,ix].min()])
    i_top = np.where(data == z_top)
    x_top, y_top = data.columns[i_top[1][0]], data.index[i_top[0][0]]

    return np.array((x_min, y_min, z_min)), np.array((x_top, y_top, z_top))


def plane(p1, p2, p3): 
    v12 = p2 - p1
    v13 = p3 - p1
    cp = np.cross(v12, v13)
    a, b, c = cp
    d = np.dot(cp, p3)
    # Print eq of plane with coeffs to 2 decimals
    print(f"{a:.2f}x + {b:.2f}y + {c:.2f}z = {d:.2f}")
    return a, b, c, d

def trace(x, a, b, c, d, a_p, b_p, c_p, d_p):
    # plne = plane(p1, p2, p3)
    # Find the intersection of the plane with the surface
    y = (d_p - c_p*d - (a_p + c_p*a)*x ) / (b_p + c_p*b + c_p*c*x)
    return y

def get_trace(data, coeffs):
    x1 = data.columns.values
    y1 = data.index.values
    # 1. Get the saddle point
    sp = saddle(x1, y1, **coeffs)
    # 2. Get local minima
    p1, p2 = get_minima(data)
    # 3. Get the plane
    plne = plane(p1, p2, sp)
    # 4. Trace the plane
    trce = trace(x1, *coeffs.values(), *plne)
    return trce

# THIS WORKS IF SADDLE IS IN DOMAIN!

def get_traces(dfz, dfr, by, x_name, y_name, z):
    traces = []
    for i,row in dfr.iterrows():
        coeffs = row[['a', 'b', 'c', 'd']].to_dict()

        data = dfz[dfz[by] == row[by]].pivot(index=y_name, columns=x_name, values=z)

        trce = get_trace(data, coeffs)
        traces.append(
            pd.DataFrame({by : row[by], x_name : data.columns, 'trace' : trce})
        )
    return pd.concat(traces)