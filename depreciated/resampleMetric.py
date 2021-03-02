#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 21:32:02 2019

Builds a metric GIFTI structure and resamples it to target atlas space

INPUT: 
hemisphere:		Options: left [0] or right [1]. Currently no support to do both hemispheres in one call.
inMetric:		Full file path of metric file to be resampled (<inMetric>.metric)
outGifti:  		Full file path of output gifti file that has been resampled (<outGifti>.gii)

VARARGIN: 
resolution:  	Number of vertices in surface which we use to resample metric file ('32k' or '164k').
				DEFAULT: '32k'			

@author: SArbuckle 05/2019 (Python conversion: switt)
"""
import os
import sys
import numpy as np
import nibabel as nb
from . import makeFuncGifti

def resampleMetric(hemisphere,inMetric,outGifti,resolution='32k'):

	BASE_DIR = pathlib.Path('surfAnalysisPy').resolve()
    atlasDir = BASE_DIR.joinpath('standard_mesh')

    anatStruct = ['CortexLeft','CortexRight']
    Hem = ['L','R']

    # Check whether inMetric is a metric file
    inMetricFileName = os.path.split(inMetric)[1]
    inMetricFileExtension = os.path.splitext(inMetricFileName)[0]
    if inMetricFileExtension != '.metric':
    	sys.exit('Error: {} is not a .metric file'.format(inMetric))

   	# Check whether outGifti is a gifti file
   	outGiftiFileName = os.path.split(outGifti)[1]
   	outGiftiFileExtension = os.path.splitext(outGiftiFileName)[1]
   	if outGiftiFileExtension != '.gii':
   		sys.exit('Error: {} is a not a .gii file'.format(outGifti))

   	# Load .metric file and convert to .gii
   	metricFile = nb.load(inMetric)
   	giftiFile = makeFuncGifti.makeFuncGifti(metricFile.darrays.data,\
   		columnNames=metricFile.darrays.meta.data.value,\
   		anatomicalStruct=anatStruct[hemisphere])
   	nb.save(giftiFile,outGifti)

   	# Create filenames for inflated surface speheres.  
   	# fsaverageSphere is aligned to atlasSphere, and then outGifti is resampled.
   	fsaverageSphere = os.path.join(atlasDir,["fs_",Hem[hemisphere]],\
   		('',join(("fsaverage.",Hem[hemisphere],".sphere.164k_fs_",\
   			Hem[hemisphere],".surf.gii"))))
   	atlasSphere = os.path.join(atlasDir,"resample_fsaverage",\
   		(''.join(("fs_LR-deformed_to-fsaverage.",Hem[hemisphere],\
   			".sphere.",resolution,"_fs_LR.surf.gii"))))

   	subprocess.run(["wb_command","-metric-resample",outGifti,fsaverageSphere,atlasSphere,\
   		"ADAP_BARY_AREA",outGifti,"-area-surfs",fsaverageSphere,atlasSphere])









