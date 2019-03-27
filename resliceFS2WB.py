#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 14:10:56 2019

@author: switt
"""
import os
import re
import numpy as np
import subprocess
import nibabel as nb
from . import affineTransform


def resliceFS2WB(subjName,subjDir,atlasDir,outDir,\
                 smoothing=1,surfFiles=["white","pial","inflated"],\
                 curvFiles=["curv","sulc","area"],hemisphere=[0,1],\
                 alignSurf=[1,1,1]):
    
    hemisphere = np.array(hemisphere)
    alignSurf = np.array(alignSurf)
    
    structName = ["left","right"]
    hem = ["lh","rh"]
    Hem = ["L","R"]
    
    currentDir = os.getcwd()
    direct = os.path.join(os.getenv("FREESURFER_HOME"),"average","surf")
    
    if not subjDir:
        subjDir = os.getenv("SUBJECTS_DIR")
    
    # read freesurfer version
    verFile = os.path.join(os.getenv("FREESURFER_HOME"),"build-stamp.txt") 
    f = open(verFile, 'r')
    verStr = f.readline()
    f.close()
    verStr = verStr.replace('-v',' ')
    verStr = re.split('[ -F]+', verStr)
    fslVer = verStr[5]
    
    # create new output directory
    newDir = os.path.join(outDir,subjName)
    if not newDir:
        os.mkdir(newDir)
     
        #---------------------
        #Transform surface files to standard mesh
    numSurf = len(surfFiles)
    numCurv = len(curvFiles)
    os.chdir(os.path.join(subjDir,subjName,"surf"))
    
    #---------------------
    #Figure out the shifting of coordinate systems:
    # Freesurfer uses vertex coordinates in respect to
    # the center of the 256x256x256 image.
    # Independent of the real zero point in the original image
    # So to find a transform of the
    # Mvox2surf: Transform of voxels in 256x256 image to surface vertices
    # Mvox2space: Transform of voxel to subject space
    anaFile = os.path.join(subjDir,subjName,"mri","brain.mgz")
    p = subprocess.run(["mri_info", anaFile, "--vox2ras-tkr"],\
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE).\
                         stdout.decode('utf-8').split()
    outp = np.array(list(map(float,p)))
    Mvox2surf = outp.reshape(-1,4)


    q = subprocess.run(["mri_info", anaFile, "--vox2ras"],\
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE).\
                       stdout.decode('utf-8').split()
    outq = np.array(list(map(float,q)))
    Mvox2space = outq.reshape(-1,4)
    Msurf2space = np.matmul(Mvox2space,np.linalg.inv(Mvox2surf))
    
    #---------------------
    #Transform the surfaces from the two hemispheres
    for h in hemisphere:
        #Convert regSphere
        regSphere = '.'.join((hem[h],"sphere.reg.surf.gii"))

        subprocess.call(["mris_convert", ('.'.join((hem[h],"sphere.reg"))),regSphere])
        
    #-----------------------
    #do all the surface files
        for i in range(numSurf):
            #Set up file names
            fileName = '.'.join((hem[h],surfFiles[i],"surf.gii"))
            
            if len(subjName) == 0:
                outName = os.path.join(newDir,('.'.join((Hem[h],surfFiles[i],\
                                                         '164k.surf.gii'))))
            else:
                outName = os.path.join(newDir,('.'.join((subjName,Hem[h],\
                                                         surfFiles[i],'164k.surf.gii'))))
                
            atlasName = os.path.join(atlasDir,"resample_fsaverage",\
                                      ('.'.join(("fs_LR-deformed_to-fsaverage",\
                                                Hem[h],"sphere.164k_fs_LR.surf.gii"))))
            

            subprocess.run(["mris_convert", ('.'.join((hem[h],surfFiles[i]))),fileName])
            

            subprocess.run(["wb_command", "-surface-resample",\
                             fileName,regSphere,atlasName,\
                             "BARYCENTRIC",outName])
            

    #Not sure I have the correct CoordinateSystemTransformMatrix...
            A = nb.load(outName)
            if (alignSurf[i]):
                [A.darrays[0].coordsys.xform[:,0],A.darrays[0].coordsys.xform[:,1],\
                 A.darrays[0].coordsys.xform[:,2]]=\
                 affineTransform.affineTransform(A.darrays[0].\
                 coordsys.xform[:,0],A.darrays[0].coordsys.xform[:,1],\
                 A.darrays[0].coordsys.xform[:,2],Msurf2space)
                
            nb.save(A,outName)
            
    #-------------------------
    #Do all the curvature files
        for i in range(numCurv):
            #Set up file names
            fileName = '.'.join((hem[h],curvFiles[i],"shape.gii"))
            
            if len(subjName) == 0:
                outName = os.path.join(newDir,('.'.join((Hem[h],curvFiles[i],\
                                                         "164k.shape.gii"))))
            else:
                outName = os.path.join(newDir,('.'.join((subjName,Hem[h],\
                                                         curvFiles[i],"164k.shape.gii"))))
                
            atlasName = os.path.join(atlasDir,"resample_fsaverage",\
                                     ('.'.join(("fs_LR-deformed_to-fsaverage",\
                                               Hem[h],"sphere.164k_fs_LR.surf.gii"))))
            
            #Convert surface to Gifti
            subprocess.run(["mris_convert", "-c", ('.'.join((hem[h],curvFiles[i]))),\
                            ('.'.join((hem[h],surfFiles[0]))), fileName])
           
            subprocess.run(["wb_command", "-metric-resample",\
                             fileName, regSphere, atlasName, "BARYCENTRIC", outName])           
            
            
            
        
    

        
    
    
