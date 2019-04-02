#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 21:32:02 2019

@author: switt
"""
import os
import numpy as np
import nibabel as nb

def makeFuncGifti(data,anatomicalStruct='CortexLeft',columnNames=[]):
    
    [N,Q] = [data.shape[0], data.shape[1]]
    #Make columnNames if empty
    if data.shape[0] == 0:
        for i in range(Q):
            columnNames.append("col_{:02d}".format(i+1))
    
    
    S = nb.gifti.gifti.GiftiImage()        
    
    C = nb.gifti.gifti.GiftiMetaData()
    C.data = [nb.gifti.gifti.GiftiNVPairs(name='AnatomicalStructurePrimary',value=anatomicalStruct),\
              nb.gifti.gifti.GiftiNVPairs(name='encoding',value='XML_BASE64_GZIP')]
    
    E = nb.gifti.gifti.GiftiLabel()
    E.key = 0
    E.label= '???'
    E.red = 1.0
    E.green = 1.0
    E.blue = 1.0
    E.alpha = 0.0
    
    D = list()
    for i in range(Q):
        d = nb.gifti.gifti.GiftiDataArray()
        d.data = np.float32(data[:,i])
        d.meta.data = nb.gifti.gifti.GiftiNVPairs(name='Name',value=columnNames[i])
        d.dims = N
        d.datatype = 'NIFTI_TYPE_FLOAT32'
        d.Intent = 'NIFTI_INTENT_NONE'
        d.coordsys.dataspace = 0
        D.append(d)
        
    S._meta = C
    S._labeltable = E
    S.darrays = D
       
    return S