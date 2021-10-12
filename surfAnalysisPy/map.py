#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import os
import sys
import re
import pathlib
import subprocess
import matplotlib.pyplot as plt
import scipy.stats as ss
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import nibabel as nib
import matplotlib as mpl
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import warnings

_base_dir = os.path.dirname(os.path.abspath(__file__))
_surf_dir = os.path.join(_base_dir, 'standard_mesh')

def affine_transform(x1,x2,x3,M):
    """
    Returns affine transform of x
    INPUT:
        x1 (np-array): X-coordinate of original
        x2 (np-array): Y-coordinate of original
        x3 (np-array): Z-coordinate of original
        M (2d-array): 4x4 transformation matrix

    OUTPUT:
        x1 (np-array): X-coordinate of transform
        x2 (np-array): Y-coordinate of transform
        x3 (np-array): Z-coordinate of transform
        transformed coordinates: same form as x1,x2,x3
    """
    y1 = np.multiply(M[0,0],x1) + np.multiply(M[0,1],x2) + np.multiply(M[0,2],x3) + M[0,3]
    y2 = np.multiply(M[1,0],x1) + np.multiply(M[1,1],x2) + np.multiply(M[1,2],x3) + M[1,3]
    y3 = np.multiply(M[2,0],x1) + np.multiply(M[2,1],x2) + np.multiply(M[2,2],x3) + M[2,3]
    return (y1,y2,y3)

def reslice_fs_to_wb(subj_name, subj_dir, out_dir,\
                 smoothing=1, surf_files=["white","pial","inflated"],\
                 curv_files=["curv","sulc","area"], hemisphere=[0,1],\
                 align_surf=[1,1,1], resolution="32k"):

    """
    Resamples a registered subject surface from freesurfer average to the new
    symmetric fs_LR_32 surface, standard in workbench.  
    This allows things to happen exactly in atlas space - each vertex number
    corresponds exactly to a anatomical location. 
    For more information, see: 
    https://wiki.humanconnectome.org/download/attachments/63078513/Resampling-FreeSurfer-HCP_5_8.pdf

    INPUT: 
        subj_name (string): 
            Subject name 
        subj_dir (string): 
            Path to freesurfer's SUBJECT_DIR (or location of freesurfer output)
        out_dir (string): 
            Path to where new resampled files will be written (out_dir/subj_name)
    
    OPTIONAL: 
        hemisphere (array): 
            Resample left, right, or both hemispheres? (0=left, 1=right) 
            DEFAULT = [0,1]
        align_surf (array):
            Shift the surface to correct for freesurfer convention? 
            DEFAULT = [1,1,1] 
        surf_files (list): 
            Surface files to be resampled. 
            DEFAULT = [".white",".pial",".inflated"]
        curv_files (list): 
            Curvature files to be resampled. 
            DEFAULT = [".curv",".sulc",".area"]
        resolution (string): 
            Resolution can be either set to '164k' or '32k'. 
            DEFAULT =  "32k"
        
    OUTPUT:
        Resampled surfaces (gifti files)
    """
    
    base_dir = pathlib.Path('surfAnalysisPy').resolve()
    atlas_dir = base_dir.joinpath('standard_mesh')
    
    hemisphere = np.array(hemisphere)
    align_surf = np.array(align_surf)
    
    structName = ["left","right"]
    hem = ["lh","rh"]
    Hem = ["L","R"]
    
    currentDirectory = os.getcwd()
    freesurferAverageSurfaceDirectory = os.path.join(os.getenv("FREESURFER_HOME"),"average","surf")
    
    if not subj_dir:
        subj_dir = os.getenv("SUBJECTS_DIR")
    
    # Read in freesurfer version
    fs_version_file = os.path.join(os.getenv("FREESURFER_HOME"),"build-stamp.txt") 
    f = open(fs_version_file, 'r')
    fs_version_str = f.readline()
    f.close()
    fs_version_str = fs_version_str.replace('-v',' ')
    fs_version_str = re.split('[ -F]+', fs_version_str)
    fs_version = fs_version_str[5]
    
    # Create new output directory for subject
    subjout_dir = os.path.join(out_dir,subj_name)
    if not subjout_dir:
        os.mkdir(subjout_dir)
     
    num_surf_files = len(surf_files)
    num_curv_files = len(curv_files)
    os.chdir(os.path.join(subj_dir,subj_name,"surf"))
    
    # Figure out the shifting of coordinate systems:
    # Freesurfer uses vertex coordinates in respect to
    # the center of the 256x256x256 image, independent
    # of the real zero point in the original image.
    # vox2surfTransformMatrix: Transform of voxels in 256x256 image to surface vertices
    # vox2spaceTransformMatrix: Transform of voxel to subject space

    anat_file = os.path.join(subj_dir,subj_name,"mri","brain.mgz")
    mriInfoVox2RasTkrProcess = subprocess.run(["mri_info", anat_file, "--vox2ras-tkr"],\
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE).\
                         stdout.decode('utf-8').split()
    mriInfoVox2RasTkrProcessOutput = np.array(list(map(float,mriInfoVox2RasTkrProcess)))
    vox2surfTransformMatrix = mriInfoVox2RasTkrProcessOutput.reshape(-1,4)


    mriInfoVox2RasProcess = subprocess.run(["mri_info", anat_file, "--vox2ras"],\
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE).\
                       stdout.decode('utf-8').split()
    mriInfoVox2RasProcessOutput = np.array(list(map(float,mriInfoVox2RasProcess)))
    vox2spaceTransformMatrix = mriInfoVox2RasProcessOutput.reshape(-1,4)
    surf2spaceTransformMatrix = np.matmul(vox2spaceTransformMatrix,np.linalg.inv(vox2surfTransformMatrix))
    
    # Transform the surfaces from the two hemispheres
    for h in hemisphere:
        #Convert reg_sphere
        reg_sphere = '.'.join((hem[h],"sphere.reg.surf.gii"))

        subprocess.call(["mris_convert", ('.'.join((hem[h],"sphere.reg"))),reg_sphere])
        
    # Transform all the surface files
        for i in range(num_surf_files):
            # Set up file names
            file_name = '.'.join((hem[h],surf_files[i],"surf.gii"))
            
            if len(subj_name) == 0:
                surf_gifti_name = os.path.join(subjout_dir,('.'.join((Hem[h],surf_files[i],\
                                                         resolution,'surf.gii'))))
            else:
                surf_gifti_name = os.path.join(subjout_dir,('.'.join((subj_name,Hem[h],\
                                                         surf_files[i],resolution,'surf.gii'))))
                
            atlas_name = os.path.join(atlas_dir,"resample_fsaverage",\
                                      (''.join(("fs_LR-deformed_to-fsaverage.",\
                                                Hem[h],".sphere.",resolution,"_fs_LR.surf.gii"))))
            

            subprocess.run(["mris_convert", ('.'.join((hem[h],surf_files[i]))),file_name])
            

            subprocess.run(["wb_command", "-surface-resample",\
                             file_name,reg_sphere,atlas_name,\
                             "BARYCENTRIC",surf_gifti_name])
            
            surf_gifti = nib.load(surf_gifti_name)
            
            if (align_surf[i]):
                [surf_gifti.darrays[0].coordsys.xform[:,0],surf_gifti.darrays[0].coordsys.xform[:,1],\
                 surf_gifti.darrays[0].coordsys.xform[:,2]]=\
                 affine_transform(surf_gifti.darrays[0].\
                 coordsys.xform[:,0],surf_gifti.darrays[0].coordsys.xform[:,1],\
                 surf_gifti.darrays[0].coordsys.xform[:,2],surf2spaceTransformMatrix)
                
            nib.save(surf_gifti,surf_gifti_name)
            
    # Transform all the curvature files
        for i in range(num_curv_files):
            # Set up file names
            file_name = '.'.join((hem[h],curv_files[i],"shape.gii"))
            if len(subj_name) == 0:
                curv_gifti_name = os.path.join(subjout_dir,('.'.join((Hem[h],curv_files[i],\
                                                         resolution,"shape.gii"))))
            else:
                curv_gifti_name = os.path.join(subjout_dir,('.'.join((subj_name,Hem[h],\
                                                         curv_files[i],resolution,"shape.gii"))))
            atlas_name = os.path.join(atlas_dir,"resample_fsaverage",\
                                     (''.join(("fs_LR-deformed_to-fsaverage.",\
                                               Hem[h],".sphere.",resolution,"_fs_LR.surf.gii"))))
                
            subprocess.run(["mris_convert", "-c", ('.'.join((hem[h],curv_files[i]))),\
                            ('.'.join((hem[h],surf_files[0]))), file_name])
            subprocess.run(["wb_command", "-metric-resample",\
                             file_name, reg_sphere, atlas_name, "BARYCENTRIC", curv_gifti_name])           

def coords_to_voxelidxs(coords, vol_def):
    """
    Maps coordinates to linear voxel indices

    INPUT:
        coords (3*N matrix or 3xPxQ array):
            (x,y,z) coordinates
        voldef (nibabel object):
            Nibabel object with attributes .affine (4x4 voxel to coordinate transformation matrix from the images to be sampled (1-based)) and shape (1x3 volume dimension in voxels)

    OUTPUT:
        linidxsrs (N-array or PxQ matrix):
            Linear voxel indices
    """
    mat = np.array(vol_def.affine)

    # Check that coordinate transformation matrix is 4x4
    if (mat.shape != (4,4)):
        sys.exit('Error: Matrix should be 4x4')

    rs = coords.shape
    if (rs[0] != 3):
        sys.exit('Error: First dimension of coords should be 3')

    if (np.size(rs) == 2):
        nCoordsPerNode = 1
        nVerts = rs[1]
    elif (np.size(rs) == 3):
        nCoordsPerNode = rs[1]
        nVerts = rs[2]
    else:
        sys.exit('Error: Coordindates have %d dimensions, not supported'.format(np.size(rs)))

    # map to 3xP matrix (P coordinates)
    coords = np.reshape(coords,[3,-1])
    coords = np.vstack([coords,np.ones((1,rs[1]))])

    ijk = np.linalg.solve(mat,coords)
    ijk = np.rint(ijk)[0:3,:]
    # Now set the indices out of range to -1
    for i in range(3):
        ijk[i,ijk[i,:]>=vol_def.shape[i]]=-1
    return ijk

def vol_to_surf(
    volumes, 
    white_surf, 
    pial_surf,
    ignore_zeros=0, 
    exclude_threshold=0, 
    depths=[0,0.2,0.4,0.6,0.8,1.0],
    stats='nanmean'
    ):
    """
    Maps volume data onto a surface, defined by white and pial surface.
    Function enables mapping of volume-based data onto the vertices of a
    surface. For each vertex, the function samples the volume along the line
    connecting the white and gray matter surfaces. The points along the line
    are specified in the variable 'depths'. default is to sample at 5
    locations between white an gray matter surface. Set 'depths' to 0 to
    sample only along the white matter surface, and to 0.5 to sample along
    the mid-gray surface.

    The averaging across the sampled points for each vertex is dictated by
    the variable 'stats'. For functional activation, use 'mean' or
    'nanmean'. For discrete label data, use 'mode'.

    If 'exclude_thres' is set to a value >0, the function will exclude voxels that
    touch the surface at multiple locations - i.e. voxels within a sulcus
    that touch both banks. Set this option, if you strongly want to prevent
    spill-over of activation across sulci. Not recommended for voxels sizes
    larger than 3mm, as it leads to exclusion of much data.

    For alternative functionality see wb_command volumne-to-surface-mapping
    https://www.humanconnectome.org/software/workbench-command/-volume-to-surface-mapping

    @author joern.diedrichsen@googlemail.com, Feb 2019 (Python conversion: switt)

    INPUTS:
        volumes (list):
            List of filenames, or nibable.NiftiImage  to be mapped
        white_surf (string or nibabel.GiftiImage):
            White surface, filename or loaded gifti object
        pial_surf (string or nibabel.GiftiImage):
            Pial surface, filename or loaded gifti object
    OPTIONAL:
        ignore_zeros (bool):
            Should zeros be ignored in mapping? DEFAULT:  False
        depths (array-like):
            Depths of points along line at which to map (0=white/gray, 1=pial).
            DEFAULT: [0.0,0.2,0.4,0.6,0.8,1.0]
        stats (str or lambda function):
            function that calculates the Statistics to be evaluated.
            lambda X: np.nanmean(X,axis=0) default and used for activation data
            lambda X: scipy.stats.mode(X,axis=0) used when discrete labels are sampled. The most frequent label is assigned.
        exclude_threshold (float):
            Threshold enables the exclusion of voxels that touch the surface
            in two distinct places
            (e.g., voxels that lie in the middle of a sulcus). If a voxel projects to two separate place
            on the surface, the algorithm excludes it, if the proportion of the bigger cluster
            is smaller than the threshold. (i.e. threshold = 0.9 means that the voxel has to
            lie at least to 90% on one side of the sulcus).
            **** Currently not supported.  exclude_threshold is automatically reset to 0. ****
            DEFAULT: 0

    OUTPUT:
        mapped_data (numpy.array):
            A Data array for the mapped data
    """

    Vols = []
    firstGood = None
    depths = np.array(depths)

    if exclude_threshold != 0:
        print('Warning: exclude_threshold option currently not supported. Resetting exclude_threshold to 0.')
        exclude_threshold = 0

    numPoints = len(depths)

    white_surf_img = nib.load(white_surf)
    pial_surf_img = nib.load(pial_surf)

    c1 = white_surf_img.darrays[0].data
    c2 = pial_surf_img.darrays[0].data
    faces = white_surf_img.darrays[1].data

    numVerts = c1.shape[0]

    if c2.shape[0] != numVerts:
        sys.exit('Error: White and pial surfaces should have same number of vertices.')

    for i in range(len(volumes)):
        try:
            a = nib.load(volumes[i])
            Vols.append(a)
            firstGood = i
        except:
            print(f'File {volumes[i]} could not be opened')
            Vols.append(None)

    if firstGood is None:
        sys.exit('Error: None of the images could be opened.')

    # Get the indices for all the points being sampled
    indices = np.zeros((numPoints,numVerts,3),dtype=int)
    for i in range(numPoints):
        c = (1-depths[i])*c1.T+depths[i]*c2.T
        ijk = coords_to_voxelidxs(c,Vols[firstGood])
        indices[i] = ijk.T

   # Read the data and map it
    data = np.zeros((numPoints,numVerts))
    mapped_data = np.zeros((numVerts,len(Vols)))
    for v,vol in enumerate(Vols):
        if vol is None:
            pass
        else:
            X = vol.get_data()
            if (ignore_zeros>0):
                X[X==0] = np.nan
            for p in range(numPoints):
                data[p,:] = X[indices[p,:,0],indices[p,:,1],indices[p,:,2]]
                outside = (indices[p,:,:]<0).any(axis=1) # These are vertices outside the volume
                data[p,outside] = np.nan

            # Determine the right statistics - if function - call it
            if stats=='nanmean':
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    mapped_data[:,v] = np.nanmean(data,axis=0)
            elif stats=='mode':
                mapped_data[:,v],_ = ss.mode(data,axis=0)
            elif callable(stats):
                mapped_data[:,v] = stats(data)

    return mapped_data

def make_func_gifti_cortex(
    data, 
    anatomical_struct='CortexLeft', 
    column_names=None
    ):
    """
    Generates a function GiftiImage from a numpy array
       @author joern.diedrichsen@googlemail.com, Feb 2019 (Python conversion: switt)

    Args:
        data (np array): shape (vertices x columns) 
        anatomical_struct (str): Anatomical Structure for the Meta-data default='CortexLeft'
        column_names (list or None): List of strings for column names, default is None
    Returns:
        gifti (functional GiftiImage)
    """

    if anatomical_struct=='L':
        anatomical_struct = 'CortexLeft'
    elif anatomical_struct=='R':
        anatomical_struct = 'CortexRight'

    try:
        num_verts, num_cols = data.shape
    except: 
        data = np.reshape(data, (len(data),1))
        num_verts, num_cols  = data.shape
  
    # Make columnNames if empty
    if column_names is None:
        column_names = []
        for i in range(num_cols):
            column_names.append("col_{:02d}".format(i+1))

    C = nib.gifti.GiftiMetaData.from_dict({
    'AnatomicalStructurePrimary': anatomical_struct,
    'encoding': 'XML_BASE64_GZIP'})

    E = nib.gifti.gifti.GiftiLabel()
    E.key = 0
    E.label= '???'
    E.red = 1.0
    E.green = 1.0
    E.blue = 1.0
    E.alpha = 0.0

    D = list()
    for i in range(num_cols):
        d = nib.gifti.GiftiDataArray(
            data=np.float32(data[:, i]),
            intent='NIFTI_INTENT_NONE',
            datatype='NIFTI_TYPE_FLOAT32',
            meta=nib.gifti.GiftiMetaData.from_dict({'Name': column_names[i]})
        )
        D.append(d)

    gifti = nib.gifti.GiftiImage(meta=C, darrays=D)
    gifti.labeltable.labels.append(E)

    return gifti

def make_label_gifti_cortex(
    data, 
    anatomical_struct='CortexLeft', 
    label_names=None,
    column_names=None, 
    label_RGBA=None
    ):
    """
    Generates a label GiftiImage from a numpy array
       @author joern.diedrichsen@googlemail.com, Feb 2019 (Python conversion: switt)

    INPUTS:
        data (np.array):
             numVert x numCol data
        anatomical_struct (string):
            Anatomical Structure for the Meta-data default= 'CortexLeft' or 'L'; 'CortexRight' or 'R'
        label_names (list): 
            List of strings for label names
        column_names (list):
            List of strings for names for columns
        label_RGBA (list):
            List of rgba vectors
    OUTPUTS:
        gifti (label GiftiImage)

    """

    if anatomical_struct=='L':
        anatomical_struct = 'CortexLeft'
    elif anatomical_struct=='R':
        anatomical_struct = 'CortexRight'

    try:
        num_verts, num_cols = data.shape
    except: 
        data = np.reshape(data, (len(data),1))
        num_verts, num_cols  = data.shape

    num_labels = len(np.unique(data))

    # check for 0 labels
    zero_label = 0 in data

    # Create naming and coloring if not specified in varargin
    # Make columnNames if empty
    if column_names is None:
        column_names = []
        for i in range(num_cols):
            column_names.append("col_{:02d}".format(i+1))

    # Determine color scale if empty
    if label_RGBA is None:
        hsv = plt.cm.get_cmap('hsv',num_labels)
        color = hsv(np.linspace(0,1,num_labels))
        # Shuffle the order so that colors are more visible
        color = color[np.random.permutation(num_labels)]
        label_RGBA = np.zeros([num_labels,4])
        for i in range(num_labels):
            label_RGBA[i] = color[i]
        if zero_label:
            label_RGBA = np.vstack([[0,0,0,1], label_RGBA[1:,]])

    # Create label names
    if label_names is None:
        idx = 0
        if not zero_label:
            idx = 1
        for i in range(num_labels):
            label_names.append("label-{:02d}".format(i + idx))

    # Create label.gii structure
    C = nib.gifti.GiftiMetaData.from_dict({
        'AnatomicalStructurePrimary': anatomical_struct,
        'encoding': 'XML_BASE64_GZIP'})

    E_all = []
    for (label,rgba,name) in zip(np.arange(num_labels),label_RGBA,label_names):
        E = nib.gifti.gifti.GiftiLabel()
        E.key = label 
        E.label= name
        E.red = rgba[0]
        E.green = rgba[1]
        E.blue = rgba[2]
        E.alpha = rgba[3]
        E.rgba = rgba[:]
        E_all.append(E)

    D = list()
    for i in range(num_cols):
        d = nib.gifti.GiftiDataArray(
            data=np.float32(data[:, i]),
            intent='NIFTI_INTENT_LABEL', 
            datatype='NIFTI_TYPE_FLOAT32', # was NIFTI_TYPE_INT32
            meta=nib.gifti.GiftiMetaData.from_dict({'Name': column_names[i]})
        )
        D.append(d)

    # Make and return the gifti file
    gifti = nib.gifti.GiftiImage(meta=C, darrays=D)
    gifti.labeltable.labels.extend(E_all)
    return gifti

def get_gifti_columns(
    gifti
    ):
    """get column names from gifti

    Args: 
        gifti (str or nib obj): full path to atlas (*.label.gii or *.func.gii) or nib gifti obj
    Returns: 
        column_names (list): list of column names
    """
    if isinstance(gifti, str):
        img = nib.load(gifti)
    else:
        img = gifti

    column_names = []
    for col in img.darrays:
        col_name =  list(col.metadata.values())[0]
        column_names.append(col_name)

    return column_names

def get_gifti_labels(
    gifti
    ):
    """get gifti labels for fpath (should be *.label.gii)

    Args: 
        gifti (str or nib obj): full path to atlas (*.label.gii) or nib obj
    Returns: 
        labels (list): list of label names
    """
    if isinstance(gifti, str):
        img = nib.load(gifti)
    else:
        img = gifti

    # labels = img.labeltable.get_labels_as_dict().values()
    label_dict = img.labeltable.get_labels_as_dict()

    return list(label_dict.values())

def get_gifti_anatomical_struct(
    gifti
    ):
    """
    Returns the primary anatomical structure for a gifti object.

    INPUT:
    gifti:				Nibabel gifti object

    OUTPUT:
    anat_struct:		AnatomicalStructurePrimary attribute from gifti object

    @author: jdiedrichsen (Python conversion: switt)
    """
    N = len(gifti._meta.data)
    anat_struct = []
    for i in range(N):
        if 'AnatomicalStructurePrimary' in gifti._meta.data[i].name:
            anat_struct.append(gifti._meta.data[i].value)
    return anat_struct

def get_gifti_colors(
    gifti,
    ignore_0=True
    ):
    """get gifti labels(should be *.label.gii)

    Args: 
        gifti (str or nibabel gifti obj): full path to atlas or gifti object
        ignore_0 (bool): default is True. ignores 0 index
    Returns: 
        rgba (np array): shape num_labels x num_rgba
        cpal (matplotlib color palette)
        cmap (matplotlib colormap)
    """
    if isinstance(gifti, str):
        img = nib.load(gifti)
    else:
        img = gifti

    labels = img.labeltable.labels

    rgba = np.zeros((len(labels),4))
    for i,label in enumerate(labels):
        rgba[i,] = labels[i].rgba
    
    if ignore_0:
        rgba = rgba[1:]
        labels = labels[1:]

    cmap = LinearSegmentedColormap.from_list('mylist', rgba, N=len(rgba))
    mpl.cm.register_cmap("mycolormap", cmap)
    cpal = sns.color_palette("mycolormap", n_colors=len(rgba))

    return rgba, cpal, cmap

def plotmap(
    data, 
    surf, 
    underlay=None,
    undermap='Greys', 
    underscale=None, 
    overlay_type='func', 
    threshold=None,
    cmap=None,
    cscale=None, 
    borders=None, 
    alpha=1.0,
    outputfile=None,
    render='matplotlib'
    ):
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
    if type(surf) is nib.gifti.gifti.GiftiImage:
        flatsurf = surf
    elif type(surf) is str:
        flatsurf = nib.load(surf)
    else: 
        raise NameError('surf needs to be a file name or a GiftiImage')
    vertices = flatsurf.darrays[0].data
    faces    = flatsurf.darrays[1].data

    # Underlay 
    if underlay is not None:
        # Load underlay and assign color
        if type(underlay) is not np.ndarray:
            underlay = nib.load(underlay).darrays[0].data
        underlay_color = _map_color(faces, underlay, underscale, undermap)
    else: 
        # set the underlay to white
        underlay_color = np.ones((faces.shape[0],4)) 

    # Load the overlay if it's a string
    if type(data) is str:
        data = nib.load(data)

    # If it is a giftiImage, figure out colormap
    if type(data) is nib.gifti.gifti.GiftiImage:
        if overlay_type == 'label':
            _, _, cmap = get_gifti_colors(data)
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
    ax = _render_matplotlib(vertices, faces, face_color, borders)
    return ax

def _map_color(
    faces, 
    data, 
    scale, 
    cmap=None, 
    threshold=None
    ):
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

def _render_matplotlib(
    vertices,
    faces,
    face_color, 
    borders
    ):
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

def plot_colorbar(
    gifti, 
    ):
    """plots colorbar for *.label.gii file
        
    Args:
        gifti (str or nib gifti obj): full path to *.label.gii or nibabel gifti obj
    """
    if isinstance(gifti, str):
        img = nib.load(gifti)
    else:
        img = gifti

    plt.figure()
    fig, ax = plt.subplots(figsize=(1,10)) # figsize=(1, 10)
    # fig, ax = plt.figure()

    rgba, cpal, cmap = get_gifti_colors(img)
    labels = get_gifti_labels(img)

    bounds = np.arange(cmap.N + 1)

    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb3 = mpl.colorbar.ColorbarBase(ax, cmap=cmap.reversed(cmap), 
                                    norm=norm,
                                    ticks=bounds,
                                    format='%s',
                                    orientation='vertical',
                                    )
    cb3.set_ticklabels(labels[::-1])  
    cb3.ax.tick_params(size=0)
    cb3.set_ticks(bounds+.5)
    cb3.ax.tick_params(axis='y', which='major', labelsize=30)

    return cb3

