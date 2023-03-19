import numpy as np
import os
import sys
import nibabel as nb
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import SUITPy as suit
import warnings

_base_dir = os.path.dirname(os.path.abspath(__file__))
_surf_dir = os.path.join(_base_dir, 'standard_mesh')

def plotmap(
        data,
        surf=None,
        underlay=None,
        undermap='gray',
        underscale=[-1.5, 1],
        overlay_type='func',
        threshold=None,
        cmap=None,
        cscale=None,
        label_names=None,
        borders=None,
        bordercolor = 'k',
        bordersize = 2,
        alpha=1.0,
        render='matplotlib',
        hover = 'auto',
        new_figure=False,
        colorbar=False,
        cbar_tick_format="%.2g",
        backgroundcolor = 'w',
        frame = None
        ):
    """Plot activity on a flatmap

    Args:
        data (np.array, giftiImage, or name of gifti file):
            Data to be plotted, should be a 28935x1 vector
        surf (str or giftiImage):
            surface file for flatmap, or ('fs32k_L','fs32k_R')
        underlay (str, giftiImage, or np-array):
            Full filepath of the file determining underlay coloring (default: sulc for standard surface)
        undermap (str)
            Matplotlib colormap used for underlay (default: gray)
        underscale (array-like)
            Colorscale [min, max] for the underlay (default: [-1, 0.5])
        overlay_type (str)
            'func': functional activation (default)
            'label': categories
            'rgb': RGB(A) values (0-1) directly specified. Alpha is optional
        threshold (scalar or array-like)
            Threshold for functional overlay. If one value is given, it is used as a positive threshold.
            If two values are given, an positive and negative threshold is used.
        cmap (str)
            A Matplotlib colormap or an equivalent Nx3 or Nx4 floating point array (N rgb or rgba values). (defaults to 'jet' if none given)
        label_names (list)
            labelnames (default is None - extracts from .label.gii )
        borders (str)
            Full filepath of the borders txt file 
        bordercolor (char or matplotlib.color)
            Color of border - defaults to 'k'
        bordersize (int)
            Size of the border points - defaults to 2
        cscale (int array)
            Colorscale [min, max] for the overlay, valid input values from -1 to 1 (default: [overlay.max, overlay.min])
        alpha (float)
            Opacity of the overlay (default: 1)
        render (str)
            Renderer for graphic display 'matplot' / 'plotly'. Dafault is matplotlib
        hover (str)
            When renderer is plotly, it determines what is displayed in the hover label: 'auto', 'value', or None
        new_figure (bool)
            If False, plot renders into matplotlib's current axis. If True, it creates a new figure (default=True)
        colorbar (bool)
            By default, colorbar is not plotted into matplotlib's current axis (or new figure if new_figure is set to True)
        cbar_tick_format : str, optional
            Controls how to format the tick labels of the colorbar, and for the hover label.
            Ex: use "%i" to display as integers.
            Default='%.2g' for scientific notation.

    Returns:
        ax (matplotlib.axis)
            If render is matplotlib, the function returns the axis
        fig (plotly.go.Figure)
            If render is plotly, it returns Figure object

    """
    if surf=='fs32k_L':
        surf = os.path.join(_surf_dir,'fs_L','fs_LR.32k.L.flat.surf.gii')
        if underlay is None: 
            underlay = os.path.join(_surf_dir,'fs_L','fs_LR.32k.L.shape.gii')
    elif surf=='fs32k_R':
        surf = os.path.join(_surf_dir,'fs_R','fs_LR.32k.R.flat.surf.gii')
        if underlay is None: 
            underlay = os.path.join(_surf_dir,'fs_R','fs_LR.32k.R.shape.gii')

    fig = suit.flatmap.plot(data,surf,
        underlay,undermap,underscale,
        overlay_type,threshold,cmap,cscale,label_names,
        borders,bordercolor,bordersize,
        alpha,render,hover,new_figure,colorbar,
        cbar_tick_format,backgroundcolor)
    return fig