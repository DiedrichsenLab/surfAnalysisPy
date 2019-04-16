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
    
    C = nb.gifti.GiftiMetaData.from_dict({
    'AnatomicalStructurePrimary': anatomicalStruct,
    'encoding': 'XML_BASE64_GZIP'})
    
    E = nb.gifti.gifti.GiftiLabel()
    E.key = 0
    E.label= '???'
    E.red = 1.0
    E.green = 1.0
    E.blue = 1.0
    E.alpha = 0.0
    
    D = list()
    for i in range(Q):
        d = nb.gifti.GiftiDataArray(
        data=np.float32(data[:, i]),
        intent='NIFTI_INTENT_NONE',
        datatype='NIFTI_TYPE_FLOAT32',
        meta=nb.gifti.GiftiMetaData.from_dict({'Name': columnNames[i]})
        )
        D.append(d)
        
    S = nb.gifti.GiftiImage(meta=C, darrays=D)
    S.labeltable.labels.append(E)
       
    return S