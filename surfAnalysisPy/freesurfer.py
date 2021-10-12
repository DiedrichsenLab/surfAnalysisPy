import numpy as np
import os
import sys
import re
from pathlib import Path
import nibabel as nib
import subprocess
import warnings

def freesurfer_reconall(
    subj_id, 
    fpath, 
    subject_dir=None
    ):
    """Run freesurfer routine on `subj_id`

    Calls freesurfer `recon-all`, which performs, all of the freeSurfer cortical reconstruction process

    Args: 
        subj_id (str): ex: 'sub-01'
        fpath (str): full path to T1-weighted anatomical nifti
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