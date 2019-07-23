#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 17:03:22 2019

Linear indices from multiple subscripts; a generalization of sub2ind

Function returns the linear indices for each row of subscripts in 
pos based on a matrix with dimensions siz.

INPUT:
siz:        Nibabel gifti object attribute .shape (1x3 volume dimension in voxels)        
pos:        Matrix of subindices referring to a single element

OUTPUT:
ids:        Vector of linear indices

If siz=[s1, ..., sn] refers to a matrix of size s1 x ... x sN, and pos is
an M x N where each row contains the subindices referring to a single
element, then ids is a M x 1 vector with linear indices.

For indices in pos that are out of bounds, the corresponding element in ids
is NaN. If pos is of type int32, however, the corresponding element is set
to zero.

If pos=[p1 ... pN] (i.e. M==1) and [p1,...,pN]=sub2ind(SIZ,IS), then
    [p1,...,pN]=SURFING_SUBS2INDS(siz,ids)

@author: NNO May 2009, updated June 2010 (Python conversion: switt)

"""

import os
import sys
import numpy as np
import nibabel as nb

def subs2inds(siz,pos):
    
    dimCount = np.shape(siz)[0]
    [posCount,dimCount2] = np.shape(pos)
    
    if (dimCount != dimCount2):
        sys.exit('Error: Number of dimensions do not match.')
    
    # Make sure pos and siz have same data type class    
    posClass = pos.dtype
    sizClass = siz.dtype 
    if (posClass != sizClass):
        siz = siz.astype(posClass)
    
    # Multiplication factors for the different positions    
    mply = np.zeros([1,dimCount],dtype=clpos)
    mply[0,0] = 1
    for k in range(dimCount-1):
        mply[0][k+1] = mply[0][k]*siz[k]
        
    alot = 1e6
    if posCount>alot:
        ids = np.zeros([posCount,1],dtype=clpos)
        for k in range(1,posCount,alot):
            idxs = range(k,np.min([k+alot-1,posCount]))
            ids[idxs] = subs2inds(siz,pos[idxs,:])
        return
    
    ids = np.sum(np.multiply((pos-1),np.tile(mply,(posCount,1))),axis=1)+1
    
    beyond = np.tile(siz,(posCount,1))-pos
    outOfBounds = np.sum(np.logical_or(pos<1,beyond<0),axis=1)>0
    if len(np.where(outOfBounds)[0]):
        ids[np.where(outOfBounds)]=np.nan
    
    return ids
    