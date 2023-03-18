import numpy as np
import surfAnalysisPy as surf
import nibabel as nb
import os

_base_dir = os.path.dirname(os.path.abspath(__file__))
_surf_dir = os.path.join(_base_dir, 'standard_mesh')
_individ_dir = os.path.join(_base_dir, 'example_individual')

s02white = os.path.join(_individ_dir,'s02.L.white.32k.surf.gii')
s02pial = os.path.join(_individ_dir,'s02.L.pial.32k.surf.gii')
s02func = os.path.join(_individ_dir,'con_motorImagery-average_4.nii')


D = surf.map.vol_to_surf([s02func],s02pial,s02white)
D.shape
fig=surf.plot.plotmap(D,'fs32k_R',cscale = [-5,5],render='matplotlib')
pass