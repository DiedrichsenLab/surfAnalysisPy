import numpy as np
import os
import sys
import re
from pathlib import Path
import nibabel as nib
import subprocess
import warnings

# Authors: Maedbh King, Suzanne Witt

def _affine_transform(x1, x2, x3, M):
    """
    Returns affine transform of x
    
    Args:
        x1 (np-array): X-coordinate of original
        x2 (np-array): Y-coordinate of original
        x3 (np-array): Z-coordinate of original
        M (2d-array): 4x4 transformation matrix

    Returns:
        x1 (np-array): X-coordinate of transform
        x2 (np-array): Y-coordinate of transform
        x3 (np-array): Z-coordinate of transform
        transformed coordinates: same form as x1,x2,x3
    """
    y1 = np.multiply(M[0,0],x1) + np.multiply(M[0,1],x2) + np.multiply(M[0,2],x3) + M[0,3]
    y2 = np.multiply(M[1,0],x1) + np.multiply(M[1,1],x2) + np.multiply(M[1,2],x3) + M[1,3]
    y3 = np.multiply(M[2,0],x1) + np.multiply(M[2,1],x2) + np.multiply(M[2,2],x3) + M[2,3]
    
    return (y1,y2,y3)

def _coords_to_voxelidxs(
    coords, 
    vol_def
    ):
    """
    Maps coordinates to linear voxel indices
    
    Args:
        coords (3*N matrix or 3xPxQ array): (x,y,z) coordinates
        vol_def (nibabel object): nibabel object with attributes .affine (4x4 voxel to coordinate transformation matrix from the images to be sampled (1-based)) and shape (1x3 volume dimension in voxels)
    Retuerns:
        linidxsrs (N-array or PxQ matrix): linear voxel indices
    """
    mat = np.array(vol_def.affine)

    # Check that coordinate transformation matrix is 4x4
    if (mat.shape != (4,4)):
        sys.exit('Error: Matrix should be 4x4')

    rs = coords.shape
    if (rs[0] != 3):
        sys.exit('Error: First dimension of coords should be 3')

    if (np.size(rs) == 2):
        nCoordsPerNode = 1
        nVerts = rs[1]
    elif (np.size(rs) == 3):
        nCoordsPerNode = rs[1]
        nVerts = rs[2]
    else:
        sys.exit('Error: Coordindates have %d dimensions, not supported'.format(np.size(rs)))

    # map to 3xP matrix (P coordinates)
    coords = np.reshape(coords,[3,-1])
    coords = np.vstack([coords,np.ones((1,rs[1]))])

    ijk = np.linalg.solve(mat,coords)
    ijk = np.rint(ijk)[0:3,:]
    # Now set the indices out of range to -1
    for i in range(3):
        ijk[i,ijk[i,:]>=vol_def.shape[i]]=-1
    return ijk

def _transform_coordinates(
    subj_id, 
    subject_dir
    ):
    """Figure out the shifting of coordinate systems. 

    Freesurfer uses vertex coordinates in respect to the center of the 256x256x256 image, 
    independent of the real zero point in the original image.
    surf2space_mat: Transform of voxel to subject space

    Args: 
        subj_id (str): ex. 'sub-01'
        subject_dir (str): path to freesurfer's SUBJECT_DIR

    Returns: 
        surf2space_mat (np array)
    """
    anat_fpath = os.path.join(subject_dir, subj_id, "mri", "brain.mgz")
    vox2rastkr_process = subprocess.run(["mri_info", anat_fpath, "--vox2ras-tkr"],\
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE).\
                            stdout.decode('utf-8').split()
    vox2rastkr_output = np.array(list(map(float, vox2rastkr_process)))
    vox2surf_mat = vox2rastkr_output.reshape(-1,4)


    vox2ras_process = subprocess.run(["mri_info", anat_fpath, "--vox2ras"],\
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE).\
                        stdout.decode('utf-8').split()
    vox2ras_output = np.array(list(map(float,vox2ras_process)))
    vox2space_mat = vox2ras_output.reshape(-1,4)
    surf2space_mat = np.matmul(vox2space_mat,np.linalg.inv(vox2surf_mat))

    return surf2space_mat

def freesurfer_reconall(
    subj_id, 
    fpath, 
    subject_dir=None
    ):
    """Run freesurfer routine on `subj_id` and `fpath`

    Calls freesurfer `recon-all`, which performs, all of the freeSurfer cortical reconstruction process

    Args: 
        subj_id (str): ex: 'sub-01'
        fpath (str): full path to anatomical nifti. ex: '*/<subj_id>_desc-preproc_T1w.nii'
        subject_dir (str or None): If None, path to freesurfer's SUBJECT_DIR (or location of freesurfer output)
    """
    old_dir = os.getenv("SUBJECTS_DIR")

    # set environment variables
    if not subject_dir:
        subject_dir = os.getenv("SUBJECTS_DIR")
    else:
        os.environ["SUBJECTS_DIR"] = str(subject_dir)

    subprocess.run(["recon-all", '-s', subj_id, '-i', fpath, '-all', '-cw256'])

    os.environ["SUBJECTS_DIR"] = str(old_dir)

def freesurfer_register_xhem(
    subj_id, 
    subject_dir=None, 
    hemispheres=['lh', 'rh']
    ):
    """Does xhemi-alignment of the subject to the fsaverage_sym template 

    Calls freesurfer `surfreg`

    Args: 
        subj_id (str): subj_id
        subject_dir (str or None): If None, path to freesurfer's SUBJECT_DIR (or location of freesurfer output)
        hemispheres (list of str): default is ['lh', 'rh']
    """
    old_dir = os.getenv("SUBJECTS_DIR")

    # set environment variables
    if not subject_dir:
        subject_dir = os.getenv("SUBJECTS_DIR")
    else:
        os.environ["SUBJECTS_DIR"] = str(subject_dir)

    for hem in hemispheres:

        if hem=='lh':
            subprocess.run(["surfreg", '--s', subj_id, '--t', "fsaverage_sym", '--lh']) 
        
        elif hem=='rh':
            subprocess.run(["xhemireg", '--s', subj_id])
            subprocess.run(["surfreg", '--s', subj_id, '--t', "fsaverage_sym", '--lh', '--xhemi']) 

    os.environ["SUBJECTS_DIR"] = str(old_dir)

def freesurfer_mapicosahedron_xhem(
    subj_id,   
    subject_dir=None, 
    file_type='surf',
    files=['white','pial','inflated'],
    hemispheres=['lh', 'rh'],
    smoothing=1
    ):
    """Resampels a registered subject surface to a regular isocahedron. 
    
    Calls freesurfer `mri_surf2surf`. This allows things to happen exactly in atlas space - each vertex number
    corresponds exactly to a anatomical location. Makes a new folder, called x<subj> that contains the remapped subject

    Args: 
        subj_id (str): subject name
        subject_dir (str or None): If None, path to freesurfer's SUBJECT_DIR (or location of freesurfer output)
        file_type (str): 'surf' or 'curv'
        files (list of str): surf files are: ['white','pial','inflated']; curv files are: ['curv','sulc','area']
        hemispheres (list of str): default is ['lh', 'rh']
        smoothing (int): default is 1. 
    """

    direct = os.path.join(os.getenv("FREESURFER_HOME"), "average",  "surf")
    old_dir = os.getenv("SUBJECTS_DIR")

    # set environment variables
    if not subject_dir:
        subject_dir = os.getenv("SUBJECTS_DIR")
    else:
        os.environ["SUBJECTS_DIR"] = str(subject_dir)

    # create new directory
    new_dir = os.path.join(subject_dir, f'x{subj_id}', 'surf')
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    # loop over hemispheres
    for hem in hemispheres:
        if hem=='lh':
            orig_dir = os.path.join(subject_dir, subj_id, 'surf')
        elif hem=='rh':
            orig_dir = os.path.join(subject_dir, subj_id, 'xhemi', 'surf')
        os.chdir(orig_dir)
    
        # read .reg file
        _, vertices_sphere = nib.freesurfer.io.read_geometry(os.path.join(direct, 'lh.sphere.reg'))

        # read surf/curv files
        for file in files:

            if file_type=='surf':
                vertices, _ = nib.freesurfer.io.read_geometry(os.path.join(orig_dir, f'lh.{file}'))
            elif file_type=='curv':
                vertices = nib.freesurfer.io.read_morph_data(os.path.join(orig_dir, f'lh.{file}'))
                try: 
                    col, row = vertices.shape
                except:
                    vertices = np.reshape(vertices, (len(vertices), 1))

            # loop over xyz
            xyz_all = []
            for xyz in np.arange(vertices.shape[1]):

                # write temp curv data to disk
                nib.freesurfer.io.write_morph_data(os.path.join(orig_dir, 'lh.temp.curv'), vertices[:, xyz], 10000)
                
                if hem=='lh':
                    subprocess.run(["mri_surf2surf", '--srcsubject', subj_id, '--hemi', 'lh',
                                    '--trgsubject', 'ico', '--trgicoorder', '7', 
                                    '--surfreg', 'fsaverage_sym.sphere.reg', 
                                    '--srcsurfval', 'temp.curv', '--src_type', 'curv',
                                    '--trgsurfval', 'tempN.curv', '--trg_type', 'curv', 
                                    '--mapmethod', 'nnf', '--nsmooth-out', str(smoothing)
                                    ])
                elif hem=='rh':
                    subprocess.run(["mri_surf2surf", '--srcsubject', f'{subj_id}/xhemi', '--hemi', 'lh',
                                    '--trgsubject', 'ico', '--trgicoorder', '7', 
                                    '--surfreg', 'fsaverage_sym.sphere.reg', 
                                    '--srcsurfval', 'temp.curv', '--src_type', 'curv',
                                    '--trgsurfval', 'tempN.curv', '--trg_type', 'curv', 
                                    '--mapmethod', 'nnf', '--nsmooth-out', str(smoothing)
                                    ])

                # read temp curv data
                xyz_temp = nib.freesurfer.io.read_morph_data(os.path.join(orig_dir, 'lh.tempN.curv')) 
                xyz_all.append(xyz_temp)

            # stack the temp files
            xyz_all = np.stack(xyz_all)
            xyz_all = np.reshape(xyz_all, (xyz_all.shape[1], xyz_all.shape[0]))

            # write out new files
            if file_type=='surf':
                nib.freesurfer.io.write_geometry(os.path.join(new_dir, f'{hem}.{file}'), xyz_all, vertices_sphere)
            elif file_type=='curv':
                nib.freesurfer.io.write_morph_data(os.path.join(new_dir, f'{hem}.{file}'), xyz_all[:,0], vertices_sphere.shape[0])

        # delete temporary files
        file = os.path.join(orig_dir, 'lh.temp.curv')
        if os.path.exists(file):
            os.remove(file)

        file = os.path.join(orig_dir, 'lh.tempN.curv')
        if os.path.exists(file):
            os.remove(file)

    # reset subjects dir
    os.environ["SUBJECTS_DIR"] = str(old_dir)

def reslice_surf_fs_to_wb(
    subj_id,
    mesh_dir,
    subject_dir=None,
    out_dir=None, 
    surfaces=['white', 'pial', 'inflated'],
    hemispheres=['lh', 'rh'],
    align_surfaces=[True, True, True],
    resolution='32k'
    ):
    """Resamples a registered subject surface from freesurfer average to the new symmetric fs_LR_32 surface, standard in workbench.  
    
    This allows things to happen exactly in atlas space - each vertex number corresponds exactly to a anatomical location. 
    For more information, see: https://wiki.humanconnectome.org/download/attachments/63078513/Resampling-FreeSurfer-HCP_5_8.pdf

    Args: 
        subj_id (str): subject name 
        mesh_dir (str): path to where standard meshes are saved
        subject_dir (str or None): freesurfer's SUBJECT_DIR
        out_dir (str or None): path to where new resampled files will be written 
        surfaces (list of str): surf files to be resampled: ['white', 'pial', 'inflated']
        hemispheres (list of str): default is ['lh', 'rh']
        align_surf (list of int): shift the surface to correct for freesurfer convention. default is [True, True, True]
        resolution (str): default is 32k. other option is '164k'
    
    Returns:
        Resampled surfaces (gifti files)
    """
    old_dir = os.getenv("SUBJECTS_DIR")

    # set environment variables
    if not subject_dir:
        subject_dir = os.getenv("SUBJECTS_DIR")
    else:
        os.environ["SUBJECTS_DIR"] = str(subject_dir)

    # Create new output directory for subject
    subj_out_dir = os.path.join(out_dir, subj_id)
    if not os.path.exists(subj_out_dir):
        os.makedirs(subj_out_dir)

    # change dir
    os.chdir(os.path.join(subject_dir, subj_id, 'surf'))

    # transform coordinate system
    surf2space_mat = _transform_coordinates(subj_id, subject_dir)

    # convert reg
    for hem in hemispheres:

        if hem=='lh':
            hemi = 'L'
        elif hem=='rh':
            hemi = 'R'
            
        reg_sphere = f'{hem}.sphere.reg.surf.gii'
        subprocess.run(["mris_convert", f'{hem}.sphere.reg', reg_sphere])

        # Transform all the surface files
        for (surf, align) in zip(surfaces, align_surfaces):
            # Set up file names
            surf_name = f'{hem}.{surf}.surf.gii'
            surf_fpath = os.path.join(subj_out_dir, f'{subj_id}.{hemi}.{surf}.{resolution}.surf.gii')
            atlas_name = os.path.join(mesh_dir, 'resample_fsaverage',\
                                    f'fs_LR-deformed_to-fsaverage.{hemi}.sphere.{resolution}_fs_LR.surf.gii')
            
            # run subproceses
            subprocess.run(["mris_convert", f'{hem}.{surf}', surf_name])
            subprocess.run(["wb_command", "-surface-resample",\
                            surf_name, reg_sphere, atlas_name,\
                            "BARYCENTRIC", surf_fpath])
            # load gifti
            surf_gii = nib.load(surf_fpath)
            
            if align:
                [surf_gii.darrays[0].coordsys.xform[:,0], surf_gii.darrays[0].coordsys.xform[:,1],\
                surf_gii.darrays[0].coordsys.xform[:,2]]=_affine_transform(surf_gii.darrays[0].\
                coordsys.xform[:,0], surf_gii.darrays[0].coordsys.xform[:,1],\
                surf_gii.darrays[0].coordsys.xform[:,2], surf2space_mat)
                
            nib.save(surf_gii, surf_fpath)
    
    os.environ["SUBJECTS_DIR"] = str(old_dir)

def reslice_curv_fs_to_wb(
    subj_id,
    mesh_dir,
    subject_dir=None,
    out_dir=None, 
    surfaces=['curv', 'sulc', 'area'],
    hemispheres=['lh', 'rh'],
    resolution='32k'
    ):
    """Resamples a registered subject curv file from freesurfer average to the new symmetric fs_LR_32 surface, standard in workbench.  
    
    This allows things to happen exactly in atlas space - each vertex number corresponds exactly to a anatomical location. 
    For more information, see: https://wiki.humanconnectome.org/download/attachments/63078513/Resampling-FreeSurfer-HCP_5_8.pdf

    Args: 
        subj_id (str): subject name 
        mesh_dir (str): path to where standard meshes are saved
        subject_dir (str): path to freesurfer's SUBJECT_DIR
        out_dir (str): path to where new resampled files will be written
        surfaces (list of str): surf files to be resampled: ['white', 'pial', 'inflated']
        hemispheres (list of str): default is ['lh', 'rh']
        resolution (str): default is 32k. other option is '164k'
    
    Returns:
        Resampled curv files (gifti files)
    """
    old_dir = os.getenv("SUBJECTS_DIR")

    # set environment variables
    if not subject_dir:
        subject_dir = os.getenv("SUBJECTS_DIR")
    else:
        os.environ["SUBJECTS_DIR"] = str(subject_dir)

    if not out_dir:
        out_dir = subject_dir

    # Create new output directory for subject
    subj_out_dir = os.path.join(out_dir, subj_id)
    if not os.path.exists(subj_out_dir):
        os.makedirs(subj_out_dir)

    # change dir
    os.chdir(os.path.join(subject_dir, subj_id, 'surf'))

    for hem in hemispheres:
    
        if hem=='lh':
            hemi = 'L'
        elif hem=='rh':
            hemi = 'R'

        # convert reg
        reg_sphere = f'{hem}.sphere.reg.surf.gii'
        subprocess.run(["mris_convert", f'{hem}.sphere.reg', reg_sphere])

        # Transform all the curv files
        for curv in surfaces:
            # Set up file names
            curv_fname = f'{hem}.{curv}.shape.gii'
            curv_fpath = os.path.join(subj_out_dir, f'{subj_id}.{hemi}.{curv}.{resolution}.shape.gii')
            atlas_fpath = os.path.join(mesh_dir,"resample_fsaverage",\
                                    f'fs_LR-deformed_to-fsaverage.{hemi}.sphere.{resolution}_fs_LR.surf.gii')         

            # run subprocesses
            subprocess.run(["mris_convert", "-c", f'{hem}.{curv}', f'{hem}.white', curv_fname])
            subprocess.run(["wb_command", "-metric-resample", curv_fname, reg_sphere, atlas_fpath, "BARYCENTRIC", curv_fpath])

    os.environ["SUBJECTS_DIR"] = str(old_dir)
