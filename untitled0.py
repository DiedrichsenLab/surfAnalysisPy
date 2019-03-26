#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 15:14:05 2019

@author: switt
"""
def function():
    #Reorder into new metric files
    for file in range(inputcol):
        OutData = np.empty((numrows,len(M))*np.nan
        for col in range(len(M)):
            if len(M[col]) ~= 0:
                Outdata[:,col] = M[col].darrays(:,inputcol[file])
        if (len(outfilenames) == 0 or len(outfilenames)<file):
            outfilenames[file] = print("{}.func.gii".\
                        format(inputColumnNames[inputcol[file]]))
        O = surfAnalysisPy.makeFuncGifti(OutData,\
                                             columnNames=outcolnames,\
                                             anatomicalStruct=anatStruct)
        nb.save(O,outfilenames[file])
        SumData[:,file] = groupstats(OutData)