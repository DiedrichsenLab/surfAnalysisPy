import surfAnalysisPy as surf
import nibabel as nb
import os

os.chdir('/Users/jdiedrichsen/Python/surfAnalysisPy/standard_mesh/suit')

G = surf.map.vol_to_surf('WHITE_SUIT.surf.gii','PIAL_SUIT.surf.gii',['SUIT.nii'])

