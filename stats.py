#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import nibabel as nb
import subprocess

def smooth_cifti(cifti_input, 
                 cifti_output,
                 left_surface, 
                 right_surface, 
                 surface_sigma = 2.0, 
                 volume_sigma = 0.0, 
                 direction = "COLUMN", 
                 ):
    """
    smoothes a cifti file on the direction specified by "direction"
    """
    # make up the command
    smooth_cmd = f"wb_command -cifti-smoothing {cifti_input} {surface_sigma} {volume_sigma} {direction} {cifti_output} -left-surface {left_surface} -right-surface {right_surface}"
    subprocess.run(smooth_cmd, shell=True)
    return
    
def smooth_metric():
    return

def gradient_cifti(cifti_input, cifti_output, left_surface, right_surface,
                   cerebellum=None, direction="COLUMN", surface_presmooth=2.0,
                   volume_presmooth=None, average_out=False, vectors=False):
    """Calculating the gradient of a cifti file using wb_command -cifti-gradient

    Args:
        cifti_input: path to input cifti file
        cifti_output: path to output cifti file
        left_surface: path to left surface
        right_surface: path to right surface
        cerebellum: path to cerebellum surface
        direction: direction to calculate gradient, default is COLUMN
        surface_presmooth: surface presmooth if applied, default is 2.0
                           don't apply to pre-smoothed functional connectivity
                           file, set to None
        volume_presmooth: volume presmooth if applied, default is None
        average_out: average out the gradient, default is False
        vectors: output vectors, default is False

    Returns:
        None - writes out cifti file of functional gradient

    Notes (from wb_command -cifti-gradient):
        Performs gradient calculation on each component of the cifti file, and
        optionally averages the resulting gradients.  The -vectors and
        -average-output options may not be used together.  You must specify a
        surface for each surface structure in the cifti file.  The COLUMN
        direction should be faster, and is the direction that works on dtseries.
        For dconn, you probably want ROW, unless you are using -average-output.
    """
    # make base command (minimum required arguments)
    cmd = f"wb_command -cifti-gradient {cifti_input} {direction} {cifti_output} " \
          f"-left-surface {left_surface} -right-surface {right_surface}"

    # add optional arguments if needed
    if cerebellum is not None:
        cmd += f" -cerebellum-surface {cerebellum}"
    if surface_presmooth is not None:
        cmd += f" -surface-presmooth {surface_presmooth}"
    if volume_presmooth is not None:
        cmd += f" -volume-presmooth {volume_presmooth}"
    if average_out:
        cmd += " -average-out"
    if vectors:
        cmd += " -vectors"

    subprocess.run(cmd, shell=True)
    return

# def group_giftis(fileList,inputCol=[],outColNames=[],\
#                 outFileNames=[],outFileNamePattern='{}.func.gii',\
#                 groupSummary=[],groupStats=True,replaceNaNs=True):
#     """
#     Created on Mon Mar 25 08:36:34 2019

#     Takes N func.gii files with P columns each (i.e., from each subject)
#     and creates P group func.gii files, with N columns (i.e., for each condition).
#     Also saves group summary file (mean across columns), if desired.

#     INPUTS:
#         fileList: 					   Text file with full paths to func.gii files

#     VARARGIN_OPTIONS:
#         inputCol:             Only use the column(s) specified by [colidx] from the input files
#                                             DEFAULT: use all columns from input files
#         outFileNames:  			 Name the outputfiles ['f1','f2']
#                                             DEFAULT: use outFileNamePattern to generate them
#         outColNames:					 Name the columns in output files
#                                             DEFAULT: names of input files
#         outFileNamePattern: 	 Specify pattern for naming of group files
#                                             DEFAULT: '{}.func.gii'
#         groupSummary: 				 Name of the group summary file
#                                             DEFAULT: empty (do not create group summary file)
#         groupStats: 					 Run group summary file using nanmean
#                                             DEFAULT: True (False will use mean)
#         replaceNaNs: 				 NaNs will be replaced with 0
#                                             DEFAULT: True (varagin deprecated)

#     @author: joern.diedrichsen@googlemail.com (Python conversion: switt)
#     """
# #    replaceNaNs = 0
#     funcGiftis = list()
#     numRows = []
#     numCols = []

#     with open(fileList, 'r') as f:
#         giftiFileList = f.read()

#     giftiFileList = giftiFileList.split()

#     for i in range(len(giftiFileList)):
#         fileName = (giftiFileList[i]).strip()
#         if os.path.isfile(fileName)==False:
#             print("Warning: File {} does not exist; replacing with NaNs".\
#                   format(fileName))
#             numCols.append(np.nan)
#             numRows.array(np.nan)
#         else:
#             funcGiftisTemp = nb.load(fileName)
#             funcGiftis.append(funcGiftisTemp)
#             numRows.append(funcGiftis[i].darrays[0].data.shape[0])
#             numCols.append(len(funcGiftis[i].darrays))

#     numRows = np.array(numRows)
#     numCols = np.array(numCols)

#     # Check if input func.gii files are same format
#     if np.var(numCols[~np.isnan(numCols)]) != 0:
#         sys.exit("Error: Number of columns is not the same.")
#     if np.var(numRows[~np.isnan(numRows)]) != 0:
#         sys.exit("Error: Number of vertices is not the same.")

#     # Determine the number of columns
#     if len(inputCol) == 0:
#         inputCol = [i for i in range(min(numCols))]

#     numrows = np.nanmean(numRows)

#     # Get the input column names
#     inputColumnNames = getGiftiColumnNames.getGiftiColumnNames(funcGiftis[0])

#     # Get main anatomical structure
#     anatStruct = getGiftiAnatomicalStruct.getGiftiAnatomicalStruct(funcGiftis[0])

#     # Determine output column names
#     if len(outColNames) == 0:
#         for col in range(len(giftiFileList)):
#         	giftiFileNameExt = os.path.split(giftiFileList[col])[1]
#           giftiFileName = os.path.splitext(giftiFileNameExt)[0]
#           outColNames.append(str.format("{}".format(giftiFileName)))
# #           for i in range(len(inputCol)):
# #             outColNames.append(str.format("col{:02d}".format(i+1)))

#     # Reorder into new metric files
#     for file in range(len(inputCol)):
#         OutData = np.empty([int(numrows),int(len(funcGiftis))])*np.nan
#         SumData = np.empty([int(numrows),int(len(inputCol))])*np.nan
#         for col in range(len(funcGiftis)):
#             if len(funcGiftis[col].darrays[file].data) != 0:
#                 OutData[:,col] = funcGiftis[col].darrays[inputCol[file]].data
#         if (len(outFileNames) == 0 or len(outFileNames)<file):
#             outFileNames.append(str.format("{}.func.gii".\
#                         format(inputColumnNames[inputCol[file]])))
#         funcGiftiOut = makeFuncGifti.makeFuncGifti(OutData,\
#                                              columnNames=outColNames,\
#                                              anatomicalStruct=anatStruct)
#         nb.save(funcGiftiOut,outFileNames[file])
#         if groupStats == True:
#             SumData[:,file] = np.nanmean(OutData,axis=1)
#         elif groupStats == False:
#         	SumData[:,file] = np.mean(OutData,axis=1)

#     # If Summary file name is given
#     if len(groupSummary) != 0:
#         groupSummaryGifti = makeFuncGifti.makeFuncGifti(SumData,columNames=outColNames,\
#                                anatomicalStruct=anatStruct)
#         nb.save(groupSummaryGifti,groupSummary)

# def smooth_surface(fileList,surfaceFile=[],kernelSigma=2.0,filePrefix='smooth'):
#     """
#     uses workbench to smooth files
#     """
#     if not fileList:
#         sys.exit('Error: No fileList of .func.gii files provided.')
#     if not surfaceFile:
#         sys.exit('Error: No surfaceFile provided.')
    
#     with open(fileList, 'r') as f:
#         funcFileList = f.read()

#     funcFileList = funcFileList.split()
#     numFiles = len(funcFileList)
#     smoothedFileList = list()
    
#     for i in range(numFiles):
#         funcFilePath = os.path.split(funcFileList[i])
#         smoothedFileList.append(os.path.join(funcFilePath[0],'_'.join((filePrefix,funcFilePath[1]))))
#         subprocess.run(["wb_command","-metric-smoothing",surfaceFile,funcFileList[i],kernelSigma,smoothedFileList[i]])
