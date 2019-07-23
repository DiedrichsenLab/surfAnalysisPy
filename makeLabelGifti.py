#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 21:32:02 2019

Builds a label GIFTI structure from scratch

INPUT:
data:                 NxQ data matrix (N = number of nodes; Q = number of data series)

VARARGIN
anatomicalStruct:     Anatomical Structure for GIFTI header ('CortexLeft','CortexRight','Cerebellum')
                      DEFAULT: 'CortexLeft'
labelNames:           List of names for the labels
                      DEFAULT: empty list
columnNames:          List of names of the different data columns
                      DEFAULT: empty list
labelRGBA:            Number of labels x 4 array of RGB + alpha (opacacity) values
                      DEFAULT: empty array

@author: jdiedrichsen (Python conversion: switt)
"""
import os
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt

def makeLabelGifti(data,anatomicalStruct='CortexLeft',labelNames=[],columnNames=[],labelRGBA=[]):

  [N,Q] = [data.shape[0], data.shape[1]]
  numLabels = len(np.unique(data))

  # Create naming and coloring if not specified in varargin
  # Make columnNames if empty
  if len(columnNames) == 0:
      for i in range(numLabels):
          columnNames.append("col_{:02d}".format(i+1))

  # Determine color scale if empty
  if len(labelRGBA) == 0:
    hsv = plt.cm.get_cmap('hsv',numLabels)
    color = hsv(np.linspace(0,1,numLabels))
    # Shuffle the order so that colors are more visible
    color = color[np.random.permutation(numLabels)]
    labelRGBA = np.zeros([numLabels,4])
    for i in range(numLabels):
        labelRGBA[i] = color[i]

  # Create label names
  if len(labelNames) == 0:
    for i in range(numLabels):
      labelNames.append("label-{:02d}".format(i+1))

  # Create label.gii structure    
  C = nb.gifti.GiftiMetaData.from_dict({
    'AnatomicalStructurePrimary': anatomicalStruct,
    'encoding': 'XML_BASE64_GZIP'})
    
  E = nb.gifti.gifti.GiftiLabel()
  E.key = np.arange(labelNames)
  E.label= labelNames
  E.red = labelRGBA[:,0]
  E.green = labelRGBA[:,1]
  E.blue = labelRGBA[:,2]
  E.alpha = labelRGBA[:,3]
    
  D = list()
  for i in range(Q):
    d = nb.gifti.GiftiDataArray(
    data=np.float32(data[:, i]),
    intent='NIFTI_INTENT_LABEL',
    datatype='NIFTI_TYPE_INT32',
    meta=nb.gifti.GiftiMetaData.from_dict({'Name': columnNames[i]})
    )
    D.append(d)
        
  S = nb.gifti.GiftiImage(meta=C, darrays=D)
  S.labeltable.labels.append(E)
       
  return S