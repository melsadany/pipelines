import numpy as np
import numpy as np
import nibabel as nib
from nilearn import plotting, datasets

stat_img='/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/task-resids/spm_results/spmT_0001_shared-pos.nii'
display = plotting.plot_glass_brain(None, plot_abs=False, display_mode="lzry")
hex_colors=['#ff4600','#4782b4']
rgb_colors = [tuple(int(h.lstrip('#')[i:i+2], 16)/255 for i in (0, 2, 4)) for h in hex_colors]
display.add_contours(
    stat_img, 
    levels=[3], 
    colors=rgb_colors, filled=True
    )
plotting.show()
