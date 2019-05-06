#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 13:56:51 2019

@author: switt
"""

import os
import sys
import numpy as np
import nibabel as nb
from . import subs2inds

def coords2linvoxelidxs(coords,volDef):
    
    mat = np.array(volDef.affine)
    dim = np.array(volDef.shape)
    
#    mask = np.nan
    
#    if hasattr(volDef, 'centermask'):
#        mask = volDef.centermask
#        print('using centermask')
#    elif hasattr(volDef, 'mask'):
#        mask = volDef.mask
    
    if (mat.shape != (4,4)):
        sys.exit('matrix should be 4x4')
    if (dim.shape != (3,)):
        sys.exit('dim should be a 1x3 vector')
        
    rs = coords.shape
    if (rs[0] != 3):
        sys.exit('First dimension of coords should be 3')
    
    if (np.size(rs) == 2):
        nCoordsPerNode = 1
        nVerts = rs[1]
    elif (np.size(rs) == 3):
        nCoordsPerNode = rs[1]
        nVerts = rs[2]
    else:
        sys.exit('coords have %d dimensions, not support'.format(np.size(rs)))
        
    # map to 3xP matrix (P coordinates)
    coords = np.reshape(coords,[3,-1])
    coords = np.vstack([coords,np.ones((1,rs[1]))])
    
    ijk = np.linalg.lstsq(mat,coords,rcond=-1)[0]
    ijk = np.rint(ijk)
    
    allinidxs = subs2inds(dim,np.transpose(ijk))
    linidxs = allinidxs
    
    linidxsrs = np.transpose(np.reshape(linidxs,[nCoordsPerNode,nVerts]))
    
    return linidxsrs
    
    
    
    
    
    