## ENA; conda activate cortex; ipython

import cortex
import nibabel as nib
import numpy as np
from scipy import ndimage

# Load your statistical map
stat_data=nib.load("/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/SPM/lang-iq-pcs/task-resids/spm_results/p_uncorr_0001.nii").get_fdata()
stat_data=nib.load("/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/SPM/lang-iq-pcs/task-resids/spm_results/p_uncorr_0003.nii").get_fdata()
stat_data=nib.load("/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/SPM/lang-iq-pcs/task-resids/spm_results/p_uncorr_0005.nii").get_fdata()
stat_data=nib.load("/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/SPM/lang-iq-pcs/task-resids/spm_results/p_uncorr_0007.nii").get_fdata()



# 2. Resample to 1mm (256x256x256)
# This is what pycortex's identity transform expects
target_shape = (256, 256, 256)
# Calculate exact zoom factors
zoom_factors = [
    target_shape[0] / stat_data.shape[0],
    target_shape[1] / stat_data.shape[1],
    target_shape[2] / stat_data.shape[2]
]
# Use order=1 for linear interpolation (preserves values better than nearest)
data_1mm = ndimage.zoom(stat_data, zoom_factors, order=1)
print(f"Resampled to 1mm: {data_1mm.shape}")
# 3. Create volume with identity transform
vol = cortex.Volume(data_1mm, "2mm", xfmname="identity")


cortex.webshow(vol, '2mm', 'flat')
# got to port http://localhost:XXXX as specified



######################## flattening problems
# # CRITICAL STEP: Check if flatmap surfaces exist
# try:
#     # Try to get a flat surface
#     flat_surf = cortex.db.get_surf("2mm", "flat")
#     print("Flat surfaces already exist!")
# except:
#     print("Flat surfaces don't exist. Creating them...")
#     
#     # Create flatmap surfaces
#     cortex.segment.cut_surface(
#         cx_subject="2mm", 
#         hemi="lh",
#         name="flat",
#         flatten_with="freesurfer",
#         method='couldnotgetthiswork',
#         blender_path="/snap/bin/blender"
#     )
#     cortex.segment.cut_surface(
#         cx_subject="2mm", 
#         hemi="rh",
#         name="flat",
#         flatten_with="freesurfer",
#         blender_path="/snap/bin/blender"
#     )
#     print("Flat surfaces created!")
# 


