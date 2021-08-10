import numpy as np
import os
import sys
import nibabel as nb
import matplotlib.pyplot as plt
import scipy.stats as ss
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
import warnings

_base_dir = os.path.dirname(os.path.abspath(__file__))
_surf_dir = os.path.join(_base_dir, 'standard_mesh')

def plotmap(data, surf, underlay = None,
        undermap = 'Greys', underscale = None, overlay_type = 'func', threshold = None,
        cmap = None, cscale = None, borders = None, alpha = 1.0,
        outputfile = None, render='matplotlib'):
    """
    Visualised cerebellar cortical acitivty on a flatmap in a matlab window
    INPUT:
        data (np.array, giftiImage, or name of gifti file)
            Data to be plotted 
        surf (str or giftiImage)
            Flat surface file for flatmap
        underlay (str, giftiImage, or np-array)
            Full filepath of the file determining underlay coloring (default: SUIT.shape.gii in SUIT pkg)
        undermap (str)
            Matplotlib colormap used for underlay (default: gray)
        underscale (array-like)
            Colorscale [min, max] for the underlay (default: [-1, 0.5])
        overlay_type (str)
            'func': functional activation 'label': categories 'rgb': RGB values (default: func)
        threshold (scalar or array-like)
            Threshold for functional overlay. If one value is given, it is used as a positive threshold.
            If two values are given, an positive and negative threshold is used.
        cmap (str)
            Matplotlib colormap used for overlay (defaults to 'jet' if none given)
        borders (str)
            Full filepath of the borders txt file (default: borders.txt in SUIT pkg)
        cscale (int array)
            Colorscale [min, max] for the overlay, valid input values from -1 to 1 (default: [overlay.max, overlay.min])
        alpha (float)
            Opacity of the overlay (default: 1)
        outputfile (str)
            Name / path to file to save figure (default None)
        render (str)
            Renderer for graphic display 'matplot' / 'opengl'. Dafault is matplotlib
    OUTPUT:
        ax (matplotlib.axis)
            If render is matplotlib, the function returns the axis
    """

    # load topology
    flatsurf = nb.load(surf)
    vertices = flatsurf.darrays[0].data
    faces    = flatsurf.darrays[1].data

    # Underlay 
    if underlay is not None:
        # Load underlay and assign color
        if type(underlay) is not np.ndarray:
            underlay = nb.load(underlay).darrays[0].data
        underlay_color = _map_color(faces, underlay, underscale, undermap)
    else: 
        # set the underlay to white
        underlay_color = np.ones((faces.shape[0],4)) 

    # Load the overlay if it's a string
    if type(data) is str:
        data = nb.load(data)

    # If it is a giftiImage, figure out colormap
    if type(data) is nb.gifti.gifti.GiftiImage:
        if overlay_type == 'label':
            cmap = get_gifti_colortable(data)
            data = data.darrays[0].data
        else:
            data = data.darrays[0].data

    # If 2d-array, take the first column only
    if data.ndim>1:
        data = data[:,0]
    
    # depending on data type - type cast into int
    if overlay_type=='label':
        i = np.isnan(data)
        data = data.astype(int)
        data[i]=0

    # map the overlay to the faces
    overlay_color  = _map_color(faces, data, cscale, cmap, threshold)

    # Combine underlay and overlay: For Nan overlay, let underlay shine through
    face_color = underlay_color * (1-alpha) + overlay_color * alpha
    i = np.isnan(face_color.sum(axis=1))
    face_color[i,:]=underlay_color[i,:]
    face_color[i,3]=1.0

    # If present, get the borders
    if borders is not None:
        borders = np.genfromtxt(borders, delimiter=',')

    # Render with Matplotlib
    #ax = _render_matplotlib(vertices, faces, face_color, borders)
    ax = _render_trisurf(vertices, faces, face_color, borders)
    return ax

def _map_color(faces, data, scale, cmap=None, threshold = None):
    """
    Maps data from vertices to faces, scales the values, and
    then looks up the RGB values in the color map

    Input:
        data (1d-np-array)
            Numpy Array of values to scale. If integer, if it is not scaled
        scale (array like)
            (min,max) of the scaling of the data
        cmap (str, or matplotlib.colors.Colormap)
            The Matplotlib colormap
        threshold (array like)
            (lower, upper) threshold for data display -
             only data x<lower and x>upper will be plotted
            if one value is given (-inf) is assumed for the lower
    """

    # When continuous data, scale and threshold
    if data.dtype.kind == 'f':
        # if threshold is given, threshold the data
        if threshold is not None:
            if np.isscalar(threshold):
                threshold=np.array([-np.inf,threshold])
            data[np.logical_and(data>threshold[0], data<threshold[1])]=np.nan

        # if scale not given, find it
        if scale is None:
            scale = np.array([np.nanmin(data), np.nanmax(data)])

        # Scale the data
        data = ((data - scale[0]) / (scale[1] - scale[0]))

    # Map the values from vertices to faces and integrate
    numFaces = faces.shape[0]
    face_value = np.zeros((3,numFaces),dtype = data.dtype)
    for i in range(3):
        face_value[i,:] = data[faces[:,i]]

    if data.dtype.kind == 'i':
        face_value,_ = ss.mode(face_value,axis=0)
        face_value = face_value.reshape((numFaces,))
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            face_value = np.nanmean(face_value, axis=0)

    # Get the color map
    if type(cmap) is str:
        cmap = plt.get_cmap(cmap)
    elif type(cmap) is np.ndarray:
        cmap = ListedColormap(cmap)
    elif cmap is None:
        cmap = plt.get_cmap('jet')

    # Map the color
    color_data = cmap(face_value)

    # Set missing data 0 for int or NaN for float to NaN
    if data.dtype.kind == 'f':
        color_data[np.isnan(face_value),:]=np.nan
    elif data.dtype.kind == 'i':
        color_data[face_value==0,:]=np.nan
    return color_data

def _render_matplotlib(vertices,faces,face_color, borders):
    """
    Render the data in matplotlib: This is segmented to allow for openGL renderer

    Input:
        vertices (np.ndarray)
            Array of vertices
        faces (nd.array)
            Array of Faces
        face_color (nd.array)
            RGBA array of color and alpha of all vertices
    """
    patches = []
    for i in range(faces.shape[0]):
        polygon = Polygon(vertices[faces[i],0:2], True)
        patches.append(polygon)
    p = PatchCollection(patches)
    p.set_facecolor(face_color)
    p.set_linewidth(0.0)

    # Get the current axis and plot it
    ax = plt.gca()
    ax.add_collection(p)
    xrang = [np.nanmin(vertices[:,0]),np.nanmax(vertices[:,0])]
    yrang = [np.nanmin(vertices[:,1]),np.nanmax(vertices[:,1])]

    ax.set_xlim(xrang[0],xrang[1])
    ax.set_ylim(yrang[0],yrang[1])
    ax.axis('equal')
    ax.axis('off')

    if borders is not None:
        ax.plot(borders[:,0],borders[:,1],color='k',
                marker='.', linestyle=None,
                markersize=2,linewidth=0)
    return ax

def _render_trisurf(vertices,faces,face_color,borders):
    """
    Render the data in matplotlib using plot_trisurf: This is segmented to allow for openGL renderer

    Input:
        vertices (np.ndarray)
            Array of vertices
        faces (nd.array)
            Array of Faces
        face_color (nd.array)
            RGBA array of color and alpha of all vertices
    """

    verticesTemp = np.zeros_like(vertices)
    verticesTemp[:,1] = vertices[:,0]
    verticesTemp[:,2] = vertices[:,1]
    verticesTemp[:,0] = 0
    vertices = verticesTemp

    # Get the current axis and plot it
    #ax = plt.gca(projection='3d')
    ax = plt.gca()
    p3dcollec = ax.plot_trisurf(vertices[:,0], vertices[:,1], vertices[:,2], triangles=faces, linewidth=0.0, antialiased=False, color='white')
    ax.view_init(elev=0,azim=0)
    p3dcollec.set_facecolor(face_color)
    xrang = [np.nanmin(vertices[:,1]),np.nanmax(vertices[:,1])]
    yrang = [np.nanmin(vertices[:,2]),np.nanmax(vertices[:,2])]

    ax.set_xlim(xrang[0],xrang[1])
    ax.set_ylim(yrang[0],yrang[1])
    ax.set_axis_off()

    if borders is not None:
        ax.plot(borders[:,0],borders[:,1],color='k',
                marker='.', linestyle=None,
                markersize=2,linewidth=0)
    return ax