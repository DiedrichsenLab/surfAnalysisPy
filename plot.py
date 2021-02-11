from nilearn import plotting
import nibabel as nib
from nilearn import surface
from nilearn import datasets

import os
from pathlib import Path
from scipy.stats import norm
import matplotlib.pyplot as plt

import glob

class Defaults: 
    BASE_DIR = os.path.join(Path(__file__).absolute().parent.parent.parent, "data")
    # BASE_DIR = "/global/scratch/maedbhking/projects/cerebellum_language/data"
    GLM_DIR = os.path.join(BASE_DIR, "glm_firstlevel")
    SUIT_FUNCTIONAL_DIR = os.path.join(BASE_DIR, "suit", "functional")
    SUIT_ANATOMICAL_DIR = os.path.join(BASE_DIR, "suit", "anatomical")

    # subjs
    subjs = ['sub-01', 'sub-02', 'sub-03']

    task_labels = ['covertnoun', 'covertverb', 'overtverb', 'covertadjective']
    
class VisualizeCortex:

    def __init__(self):
        self.contrast_type = "effect_size" # "z_score" and "effect_size" are the options
        self.z_map_threshold = 2 # 0.01
        self.surface_threshold = 1.0
        self.vmax = 30
        self.interactive_threshold = '10%'
        self.display_mode = 'ortho' # or ortho, l
        self.plot_abs = False
        self.glm = "glm2"
        self.interactive = False
        self.hemi = ["right", "left"]

    def plot_glass_brain(self):
        plotting.plot_glass_brain(self.img, colorbar=True, 
                                threshold=norm.isf(self.z_map_threshold),
                                title=self.title, plot_abs=self.plot_abs, 
                                display_mode=self.display_mode)
        plt.show()

    def plot_surface(self): 
        plotting.plot_surf_stat_map(self.inflated, self.texture,
                                    hemi=self.hem, threshold=self.surface_threshold,
                                    title=self.title, colorbar=True, 
                                    bg_map=self.bg_map, vmax=self.vmax)
        
        plt.show()

    def plot_interactive_surface(self):
        view = plotting.view_surf(self.inflated, self.texture, 
                                threshold=self.interactive_threshold, 
                                bg_map=self.bg_map)

        view.open_in_browser()

    def _get_surfaces(self):
        if self.hem=="right":
            inflated = self.fsaverage.infl_right
            bg_map = self.fsaverage.sulc_right
            texture = surface.vol_to_surf(self.img, self.fsaverage.pial_right)
        elif self.hem=="left":
            inflated = self.fsaverage.infl_left
            bg_map = self.fsaverage.sulc_left
            texture = surface.vol_to_surf(self.img, self.fsaverage.pial_left)
        else:
            print('hemisphere not provided')

        return inflated, bg_map, texture
    
    def visualize_subj_glass_brain(self):

        # loop over subjects
        for subj in Defaults.subjs:

            GLM_SUBJ_DIR = os.path.join(Defaults.GLM_DIR, self.glm, subj)
            os.chdir(GLM_SUBJ_DIR)
            fpaths = glob.glob(f'*{self.contrast_type}-{self.glm}.nii')
            
            for fpath in fpaths: 
                self.img = nib.load(os.path.join(GLM_SUBJ_DIR, fpath))
                self.title = f'{Path(fpath).stem}-(unc p<{self.z_map_threshold})'
                self.plot_glass_brain()

    def visualize_group_glass_brain(self):

        GLM_GROUP_DIR = os.path.join(Defaults.GLM_DIR, self.glm, "group")
        os.chdir(GLM_GROUP_DIR)
        fpaths = glob.glob(f'*{self.contrast_type}-{self.glm}.nii')

        for fpath in fpaths:
            # define contrast names
            self.img = nib.load(os.path.join(GLM_GROUP_DIR, fpath))
            self.title = f'{Path(fpath).stem}-(unc p<{self.z_map_threshold})'
            self.plot_glass_brain()

    def visualize_group_surface(self):
        self.fsaverage = datasets.fetch_surf_fsaverage()

        GLM_GROUP_DIR = os.path.join(Defaults.GLM_DIR, self.glm, "group")
        os.chdir(GLM_GROUP_DIR)
        fpaths = glob.glob(f'*{self.contrast_type}-{self.glm}.nii')

        for fpath in fpaths:
            self.img = nib.load(os.path.join(GLM_GROUP_DIR, fpath))

            # loop over hemispheres and visualize
            for self.hem in self.hemi:

                # set plot title
                self.title = f'{Path(fpath).stem}-{self.hem}'

                # return surfaces and texture
                self.inflated, self.bg_map, self.texture = self._get_surfaces()

                # plot interactive or static surface(s)
                if self.interactive:  
                    self.plot_interactive_surface()
                else:
                    self.plot_surface()

        return self.fsaverage

class VisualizeCerebellum:

    def __init__(self):
        self.contrast_type = "effect_size"
        self.glm = "glm2" 
        self.surface_threshold = 1
        self.vmax = 10

    def plot_surface(self):
        view = plotting.view_surf(surf_mesh=self.surf_mesh, 
                                surf_map=self.surf_map, 
                                colorbar=True,
                                threshold=self.surface_threshold,
                                vmax=self.vmax,
                                title=self.title) 
        # view.resize(500,500)

        view.open_in_browser()
   
    def visualize_group_surface(self):

        # get functional group dir
        SUIT_FUNCTIONAL_GROUP_DIR = os.path.join(Defaults.SUIT_FUNCTIONAL_DIR, self.glm, "group")

        os.chdir(SUIT_FUNCTIONAL_GROUP_DIR)

        # get all contrast images in gifti format
        fpaths = glob.glob(f'*{self.contrast_type}-{self.glm}.gii')

        # get surface mesh for SUIT
        self.surf_mesh = os.path.join(Defaults.SUIT_ANATOMICAL_DIR, "FLAT.surf.gii")

        # loop over all gifti files
        for fpath in fpaths:

            self.surf_map = surface.load_surf_data(fpath).astype(int)
            self.title = Path(fpath).stem

            self.plot_surface() 

# example code to visualize group contrast(s) on flat map
vis = VisualizeCerebellum()
vis.visualize_group_surface()