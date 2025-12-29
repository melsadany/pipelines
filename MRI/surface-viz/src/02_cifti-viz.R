library(ciftiTools)
library(neurobase)

map <- readnii("/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/USE_THIS/2E_044/PS-VC/sla/major/sub-2E_044_word_association__number_Tstat.nii.gz")
ciftiTools.setOption("wb_path", "/wdata/msmuhammad/workbench/connectome-wb/workbench")


# Read your CIFTI file
cii <- read_cifti("map.dscalar.nii")
# Add the correct surfaces that match your vertex counts
cii <- add_surf(cii, 
                surfL = "../../../../refs/mni-freesurfer/2mm/SUMA/lh.smoothwm.gii",
                surfR = "../../../../refs/mni-freesurfer/2mm/SUMA/rh.smoothwm.gii")

# Now visualize
view_xifti_surface(cii, color_mode = "diverging",colors  = viridis::inferno(n=500),
                   legend_embed = F,widget = T)
