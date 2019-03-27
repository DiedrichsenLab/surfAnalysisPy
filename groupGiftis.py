#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 08:36:34 2019

@author: switt
"""
import os
import numpy as np
import nibabel as nb
from . import makeFuncGifti
from . import getGiftiColumnNames
from . import getGiftiAnatomicalStruct


def groupGiftis(filelist,inputcol=[],outcolnames=[],\
                outfilenames=[],outfilenamePattern='{}.func.gii',\
                groupsummary=[],groupstats=True):
    
    replaceNans = 0
    M = list()
    numRows = []
    numCols = []
    
    inputcol = np.array(inputcol)
    
    with open(filelist, 'r') as f:
        P = f.read()
     
    P = P.split()
    
#    [d,f] = os.path.split(P[0])
#    b = os.path.splitext(f)[0]    
    
    for i in range(len(P)):
        fileName = (P[i]).strip()
        if os.path.isfile(fileName)==False:
            print("File {} does not exist; replacing with NaNs".\
                  format(fileName))
            numCols.append(np.nan)
            numRow.array(np.nan)
        else:
            A=nb.load(fileName)
            M.append(A)
#            da = M[i].darrays
            numRows.append(M[i].darrays[0].data.shape[0])
            numCols.append(len(M[i].darrays))
       
    numRows = np.array(numRows) 
    numCols = np.array(numCols)
    
    #Check if metric files are same format
    if np.var(numCols[~np.isnan(numCols)]) != 0:
        print("Number of columns is not the same")
    if np.var(numRows[~np.isnan(numRows)]) != 0:
        print("Number of vertices is not the same")
        
    #Determine the number of columns
    if inputcol.size == 0:
        inputcol = np.arange(np.amin(numCols))
    
    inputcol = np.array(inputcol)
        
    numrows = np.nanmean(numRows)
    
    #Get the input column names
    inputColumnNames = getGiftiColumnNames.getGiftiColumnNames(M[0])
    
    #Get main anatomical structure
    anatStruct = getGiftiAnatomicalStruct.getGiftiAnatomicalStruct(M[0])
    
    #Determine output column names -- NOT SURE WHAT IS GOING ON WITH THIS LOOP?
    if len(outcolnames) == 0:
        for col in range(len(P)):
            for i in range(len(inputcol)):
                a = os.path.split(P[col])[1]
                b = os.path.splitext(a)[0]
                outcolnames.append(b)
            
    #Reorder into new metric files
    for file in range(len(inputcol)):
        OutData = np.empty([int(numrows),int(len(M))])*np.nan
        for col in range(len(M)):
            if len(M[col].darrays[file].data) != 0:
                Outdata[:,col] = M[col].darrays[:,inputcol[file]]
        if (len(outfileNames) == 0 or len(outfilenames)<file):
            outfilenames[file] = print("{}.func.gii".\
                        format(inputColumnNames[inputcol[file]]))
        O = makeFuncGifti.makeFuncGifti(OutData,\
                                             columnNames=outcolnames,\
                                             anatomicalStruct=anatStruct)
        nb.save(O,outfilenames[file])
        if groupstats == True:
            SumData[:,file] = np.nanmean(OutData,axis=1)
        
    #If Summary file name is given
    if len(groupsummary) != 0:
        O = makeFuncGifti.makeFuncGifti(SumData,columNames=outcolnames,\
                               anatomicalStruct=anatStruct)
        nb.save(O,groupsummary)
        
    