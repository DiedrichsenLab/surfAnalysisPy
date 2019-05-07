#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 09:20:41 2019

@author: switt
"""

import os
import sys
import numpy as np
import nibabel as nb
from scipy import sparse
from casadi import sparsify
from . import coords2linvoxelidxs

def vol2surf(whiteSurfGifti,pialSurfGifti,funcFileList,ignoreZeros=0,excludeThres=0,columnNames=[],\
             depths=[0,0.2,0.4,0.6,0.8,1.0],interp=0,anatomicalStruct='CortexLeft',\
             stats='nanmean'):
    
    
    A = []
    firstGood = []
    indices = []
    S = struct()
    M = struct()
    depths = np.array(depths)
    

    numPoints = len(depths)
    
    whiteSurfGiftiImage = nb.load(whiteSurfGifti)
    pialSurfGiftiImage = nb.load(pialSurfGifti)
    
    c1 = whiteSurfGiftiImage.darrays[0].data
    c2 = pialSurfGiftiImage.darrays[0].data
    faces = whiteSurfGiftiImage.darrays[1].data
    
    numVerts = len(c1[1])
    
    
    class struct:
        Tiles = []
        Edges = []
        numNodes = []
        data = []

    
    if ([len(c1[1]),len(c2[1])] != [numVerts,numVerts]):
        sys.exit('Vertex matrices should be of shape vertices x 3')
        
    if (len(c1[0]) != len(c2[0])):
        sys.exit('White and pial surfaces should have same number of vertices')
        
    with open(funcFileList, 'r') as f:
        V = f.read()
     
    V = V.split()
    
    for i in range(len(V)):
        fileName = (V[i]).strip()
        try:
            a = nb.load(fileName)
            A.append(a)
            if not firstGood:
                firstGood = i  
        except:
           print('{} could not be opened'.format(fileName))
           A.append('')
           
    V = A
    
    if not firstGood:
        sys.exit('None of the images could be opened')
        
    if (len(ignoreZeros) == 1):
        ignoreZeros = np.ones(len(V))*ignoreZeros
        
    for i in range(numPoints):
        c = (1-depths[i])*np.transpose(c1)+depths[i]*np.transpose(c2)
        indices[:,i] = coords2linvoxelidxs.coords2linvoxelidxs(c,V[firstgood])
    
    # Case: excludeThres > 0
    # If necessary, now ensure that voxels are mapped on to continuous location
    # only on the flat map - exclude voxels that are assigned to two sides of
    # the sulcus    
    if (excludeThres>0):
        exclude = np.zeros([np.prod(V[1].shape),1])
        if not faces:
            sys.exit('provide faces (topology data), so that projections should be avoided')
        S.Tiles = faces
        
        # Precalculate edges for fast cluster finding
        print('Calculating edges.')
        S.numNodes = np.max(np.max(S.Tiles))
        for i in range (3):
            i1 = S.Tiles[:,i]
            i2 = S.Tiles[:,np.remainder(i,3)+1]
            S.Edges.append(i1,i2)
        S.Edges = np.unique(S.Edges,axis=0)
        
        # Generate connectivity matrix
        # csr_matrix((data, (row_ind, col_ind)), [shape=(M, N)])
        data = np.ones(len(S.Edges))
        rowInd = S.Edges[:,0]
        colInd = S.Edges[:,1]
        M = S.numNodes
        G = sparse.csr_matrix((data,(rowInd,colInd)),shape=(M,M))
        
        # Cluster the projections to the surface and exclude all voxels
        # in which the weight of the biggest cluster is not > thres
        print('Checking projections.')
        I = np.unique(indices[np.isfinite(indices)])
        for i in I:
            
            # Calculate the weight of voxel on node
            weight = np.sum(indices==1,1)
            indx = np.where(weight>0)[0]
            
            # Check whether the nodes cluster
            H = G[indx,indx]
            H = H+np.transpose(H)+sparse.identity(len(H))
            # Matlab vairable translation: p=rowPerm,q=colPerm,r=rowBlock,s=colBlock
            nb,rowPerm,colPerm,rowBlock,colBlock,coarseRowBlock,coarseColBlock = H.sparsify().btf()
            CL = np.zeros(np.shape(indx))
            for c in range(len(rowBlock)-1):
                CL[rowPerm[rowBlock[c]]:rowBlock[(c+1)]-1,0]=c
            
            if (np.max(CL)>1):
                weight_cl=np.zeros([np.max(CL),1])
                for cl in range(np.max(CL)):
                    weight_cl[cl,0] = np.sum(weight[indx[CL==cl]])
                [m,cl] = np.max(weight_cl)
                if (m/np.sum(weight_cl)>excludeThres):
                    A = indices[indx[CL!=cl],:]
                    A[A==i] = np.nan
                    indices[indx[CL!=cl],:] = A
                    exclude[i] = 1
                    print('assigned: %2.3f'.format(m/np.sum(weight_cl)))
                else:
                    A[A==i] = np.nan
                    indices[indx,:] = A
                    exclude[i] = 2
                    print('excluded: %2.3f %d'.format(m/np.sum(weight_cl),np.max(CL)))
                    
        # For debugging: save the volume showing the exluded voxels in current
        # directory
        Vexcl = V[1]
        Vexcl.set_filename = 'excl.nii'
        nb.save(np.reshape(exclude,np.array(Vexcl.shape)),'excl.nii')
        
    # Case: excludeThres = 0
    i = np.where(np.isfinite(indices))
    data = np.zeros(len(indices))*np.nan
    
    for v in range(len(V)):
        if not V[v]:
            M.data[:,v] = np.nan([len(c1),1])
        else:
            X = V[v].get_data()
            if (ignoreZeros>0):
                X[X==0] = np.nan
            data[i] = X[indices[i]]
            M.data[:,v] = np.nanmean(data)
            
    # Determine the column names based on the filenames of the volumes
    if not columnNames:
        for i in range(len(V)):
            if not V[i]:
                columnNames[i] = 'void'
            else:
                columnNames[i] = os.path.splitext(V[i].get_filename)[0]
                
    # Build a functional GIFTI structure
    C = nb.gifti.GiftiMetaData.from_dict({'AnatomicalStructurePrimary': anatomicalStruct,
                                          'encoding': 'XML_BASE64_GZIP'})
    E = nb.gifti.gifti.GiftiLabel(
    key = 0,
    label= '???',
    red = 1.0,
    green = 1.0,
    blue = 1.0,
    alpha = 0.0,
    )
    
    for i in range(V):
        D = nb.gifti.GiftiDataArray(
        data=np.single(M.data[:,i]),
        shape=len(M.data[:,i]),
        intent='NIFTI_INTENT_NONE',
        datatype='NIFTI_TYPE_FLOAT32',
        meta=nb.gifti.GiftiMetaData.from_dict({'Name': columnNames[i]})
        )
        
    S = nb.gifti.GiftiImage(meta=C, darrays=D)
    S.labeltable.labels.append(E)
    nb.save(S)
            
        
        
            
    
