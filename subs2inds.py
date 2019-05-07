#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 17:03:22 2019

@author: switt
"""

import os
import sys
import numpy as np
import nibabel as nb

def subs2inds(siz,pos):
    
    dimCount = np.shape(siz)[0]
    [posCount,dimCount2] = np.shape(pos)
    
    if (dimCount != dimCount2):
        sys.exit('Number of dimensions do not match')
        
    clpos = pos.dtype
    clsiz = siz.dtype 
    if (clpos != clsiz):
        siz = siz.astype(clpos)
        
    mply = np.zeros([1,dimCount],dtype=clpos)
    mply[0,0] = 1
    for k in range(dimCount-1):
        mply[k+1] = mply[k]*siz[k]
        
    alot = 1e6
    if posCount>alot:
        ids = np.zeros([posCount,1],dtype=clpos)
        for k in range(1,posCount,alot):
            idxs = range(k,np.min([k+alot-1,posCount]))
            ids[idxs] = subs2inds(siz,pos[idxs,:])
        return
    
    ids = np.sum((pos-1)*np.tile(mply,(posCount,1),2,'native'))+1
    
    beyond = np.tile(siz,(posCount,1),1)-pos
    outOfBounds = np.sum(((pos<1 or beyond<0),2)>0)
    ids[outOfBounds]=np.nan
    
    return ids
    