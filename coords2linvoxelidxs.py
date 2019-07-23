#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 13:56:51 2019

Maps coordinates to linear voxel indices

INPUT:
coords:         A 3xN matrix or 3xPxQ array with (x,y,z) coordinates
voldef:         Nibabel gifti object with attributes .affine (4x4 voxel to coordinate transformation matrix
                from the images to be sampled (1-based)) and .shape (1x3 volume dimension in voxels)

OUTPUT:
linidxsrs:      Nx1 vector or PxQ matrix with linear voxel indices


@author: TW,NNO May 2010 (Python conversion: switt)

"""

import os
import sys
import numpy as np
import nibabel as nb
from . import subs2inds

def coords2linvoxelidxs(coords,volDef):
    
    mat = np.array(volDef.affine)
    dim = np.array(volDef.shape)
    
    # Check that coordinate transformation matrix is 4x4
    if (mat.shape != (4,4)):
        sys.exit('Error: Matrix should be 4x4')
    # Check that volume dimension is 1x3
    if (dim.shape != (3,)):
        sys.exit('Error: Dimensions should be a 1x3 vector')
        
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
    
    ijk = np.linalg.lstsq(mat,coords,rcond=-1)[0]
    ijk = np.rint(ijk)[0:3,:]
    
    allinidxs = subs2inds.subs2inds(dim,np.transpose(ijk))
    linidxs = allinidxs
    
    linidxsrs = np.transpose(np.reshape(linidxs,[nCoordsPerNode,nVerts]))
    
    return linidxsrs
    
    
    
    
    
    