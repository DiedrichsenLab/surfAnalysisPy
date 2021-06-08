#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 14:10:56 2019

Resamples a registered subject surface from freesurfer average to the new
symmetric fs_LR_32 surface, standard in workbench.  
This allows things to happen exactly in atlas space - each vertex number
corresponds exactly to a anatomical location. 
For more information, see: 
https://wiki.humanconnectome.org/download/attachments/63078513/Resampling-FreeSurfer-HCP_5_8.pdf

INPUT: 
   subjName: subject name
   subjDir: freesurfer's SUBJECT_DIR (or location of freesurfer output)
   outDir: Location to where new resampled files will be written (outDir/subjName)
 VARARGIN: 
   'hemisphere'        : left / right or both hemispheres 
                       Default = [0,1]
   'alignSurf':        : Shift the surface to correct for freesurfer convention?
                       Default = [1,1,1] 
   'surfFiles'         : Surface files to be resampled 
                       Default = [".white",".pial",".inflated"]
   'curvFiles'         : Curvature files to be resampled 
                       Default = [".curv",".sulc",".area"]
   'resolution'        : Resolution can be either set to '164k' or '32k'.
                       Default =  "32k"

@author: jdiedrichsen (Python conversion: switt)
"""
import os
import re
import pathlib
import numpy as np
import subprocess
import nibabel as nb
from . import affineTransform


def reslice_fs_2_wb(subjName,subjDir,outDir,\
                 smoothing=1,surfFiles=["white","pial","inflated"],\
                 curvFiles=["curv","sulc","area"],hemisphere=[0,1],\
                 alignSurf=[1,1,1],resolution="32k"):
    
    BASE_DIR = pathlib.Path('surfAnalysisPy').resolve()
    atlasDir = BASE_DIR.joinpath('standard_mesh')
    
    
    hemisphere = np.array(hemisphere)
    alignSurf = np.array(alignSurf)
    
    
    structName = ["left","right"]
    hem = ["lh","rh"]
    Hem = ["L","R"]
    
    currentDirectory = os.getcwd()
    freesurferAverageSurfaceDirectory = os.path.join(os.getenv("FREESURFER_HOME"),"average","surf")
    
    if not subjDir:
        subjDir = os.getenv("SUBJECTS_DIR")
    
    # Read in freesurfer version
    freesurferVersionFile = os.path.join(os.getenv("FREESURFER_HOME"),"build-stamp.txt") 
    f = open(freesurferVersionFile, 'r')
    freesurferVersionString = f.readline()
    f.close()
    freesurferVersionString = freesurferVersionString.replace('-v',' ')
    freesurferVersionString = re.split('[ -F]+', freesurferVersionString)
    freesurferVersion = freesurferVersionString[5]
    
    # Create new output directory for subject
    subjOutDir = os.path.join(outDir,subjName)
    if not subjOutDir:
        os.mkdir(subjOutDir)
     
    numSurfFiles = len(surfFiles)
    numCurvFiles = len(curvFiles)
    os.chdir(os.path.join(subjDir,subjName,"surf"))
    
    # Figure out the shifting of coordinate systems:
    # Freesurfer uses vertex coordinates in respect to
    # the center of the 256x256x256 image, independent
    # of the real zero point in the original image.
    # vox2surfTransformMatrix: Transform of voxels in 256x256 image to surface vertices
    # vox2spaceTransformMatrix: Transform of voxel to subject space

    anatFile = os.path.join(subjDir,subjName,"mri","brain.mgz")
    mriInfoVox2RasTkrProcess = subprocess.run(["mri_info", anatFile, "--vox2ras-tkr"],\
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE).\
                         stdout.decode('utf-8').split()
    mriInfoVox2RasTkrProcessOutput = np.array(list(map(float,mriInfoVox2RasTkrProcess)))
    vox2surfTransformMatrix = mriInfoVox2RasTkrProcessOutput.reshape(-1,4)


    mriInfoVox2RasProcess = subprocess.run(["mri_info", anatFile, "--vox2ras"],\
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE).\
                       stdout.decode('utf-8').split()
    mriInfoVox2RasProcessOutput = np.array(list(map(float,mriInfoVox2RasProcess)))
    vox2spaceTransformMatrix = mriInfoVox2RasProcessOutput.reshape(-1,4)
    surf2spaceTransformMatrix = np.matmul(vox2spaceTransformMatrix,np.linalg.inv(vox2surfTransformMatrix))
    
    # Transform the surfaces from the two hemispheres
    for h in hemisphere:
        #Convert regSphere
        regSphere = '.'.join((hem[h],"sphere.reg.surf.gii"))

        subprocess.call(["mris_convert", ('.'.join((hem[h],"sphere.reg"))),regSphere])
        
    # Transform all the surface files
        for i in range(numSurfFiles):
            # Set up file names
            fileName = '.'.join((hem[h],surfFiles[i],"surf.gii"))
            
            if len(subjName) == 0:
                surfGiftiFileName = os.path.join(subjOutDir,('.'.join((Hem[h],surfFiles[i],\
                                                         resolution,'surf.gii'))))
            else:
                surfGiftiFileName = os.path.join(subjOutDir,('.'.join((subjName,Hem[h],\
                                                         surfFiles[i],resolution,'surf.gii'))))
                
            atlasName = os.path.join(atlasDir,"resample_fsaverage",\
                                      (''.join(("fs_LR-deformed_to-fsaverage.",\
                                                Hem[h],".sphere.",resolution,"_fs_LR.surf.gii"))))
            

            subprocess.run(["mris_convert", ('.'.join((hem[h],surfFiles[i]))),fileName])
            

            subprocess.run(["wb_command", "-surface-resample",\
                             fileName,regSphere,atlasName,\
                             "BARYCENTRIC",surfGiftiFileName])
            
            surfGifti = nb.load(surfGiftiFileName)
            
            if (alignSurf[i]):
                [surfGifti.darrays[0].coordsys.xform[:,0],surfGifti.darrays[0].coordsys.xform[:,1],\
                 surfGifti.darrays[0].coordsys.xform[:,2]]=\
                 affine_transform(surfGifti.darrays[0].\
                 coordsys.xform[:,0],surfGifti.darrays[0].coordsys.xform[:,1],\
                 surfGifti.darrays[0].coordsys.xform[:,2],surf2spaceTransformMatrix)
                
            nb.save(surfGifti,surfGiftiFileName)
            
    # Transform all the curvature files
        for i in range(numCurvFiles):
            # Set up file names
            fileName = '.'.join((hem[h],curvFiles[i],"shape.gii"))
            if len(subjName) == 0:
                curvGiftiFileName = os.path.join(subjOutDir,('.'.join((Hem[h],curvFiles[i],\
                                                         resolution,"shape.gii"))))
            else:
                curvGiftiFileName = os.path.join(subjOutDir,('.'.join((subjName,Hem[h],\
                                                         curvFiles[i],resolution,"shape.gii"))))
            atlasName = os.path.join(atlasDir,"resample_fsaverage",\
                                     (''.join(("fs_LR-deformed_to-fsaverage.",\
                                               Hem[h],".sphere.",resolution,"_fs_LR.surf.gii"))))
                
            subprocess.run(["mris_convert", "-c", ('.'.join((hem[h],curvFiles[i]))),\
                            ('.'.join((hem[h],surfFiles[0]))), fileName])
            subprocess.run(["wb_command", "-metric-resample",\
                             fileName, regSphere, atlasName, "BARYCENTRIC", curvGiftiFileName])           
