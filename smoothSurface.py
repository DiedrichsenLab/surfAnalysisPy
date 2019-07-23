#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 21:32:02 2019

Smooth metric files

INPUT:
fileList:             Text file of full paths to metric files

VARARGIN
surfaceFile:     	  Full path to surface file on which to smooth
					  DEFAULT: empty
kernelSigma:		  Sigma (in mm) for gaussian kernel function
					  DEFAULT: 2.0
filePrefix			  Prefix for smoothed filenames (string)
					  DEFAULT: 'smooth'

OUTPUT:

Default algorithm (GEO_GAUSS_AREA) used, no other methods implemented or supported.
More information: https://www.humanconnectome.org/software/workbench-command/-metric-smoothing

@author: EBerlot 2019/05 (Python conversion: switt)
"""
import os
import sys
import numpy as np
import nibabel as nb

def smoothSurface(fileList,surfaceFile=[],kernelSigma=2.0,filePrefix='smooth'):

	if not fileList:
		sys.exit('Error: No fileList of .func.gii files provided.')

	if not surfaceFile:
		sys.exit('Error: No surfaceFile provided.')

	with open(fileList, 'r') as f:
        funcFileList = f.read()
     
    funcFileList = funcFileList.split()

    numFiles = len(funcFileList)

    smoothedFileList = list()

    for i in range(numFiles):
    	funcFilePath = os.path.split(funcFileList[i])
    	smoothedFileList.append(os.path.join(funcFilePath[0],'_'.join((filePrefix,funcFilePath[1]))))
    	subprocess.run(["wb_command","-metric-smoothing",surfaceFile,funcFileList[i],kernelSigma,smoothedFileList[i]])

