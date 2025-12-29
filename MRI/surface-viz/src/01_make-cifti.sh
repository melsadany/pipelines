#!/bin/bash


## needed to prepare surface files for the MNI freesurfer output
export SUBJECTS_DIR=/wdata/msmuhammad/refs/mni-freesurfer/
sid="2mm"
TEMP_DIR=$SUBJECTS_DIR/$sid

# Generate SUMA files if missing, then launch SUMA
if [ ! -f $TEMP_DIR/SUMA/fsaverage_both.spec ]; then
    cd $TEMP_DIR
    @SUMA_Make_Spec_FS -sid $sid -NIFTI
fi


## get your statistic map
STAT_MAP="/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/USE_THIS/2E_044/PS-VC/sla/major/sub-2E_044_word_association__word_Tstat.nii.gz"
P_DIR=/wdata/msmuhammad/pipelines/MRI/surface-viz
cd $P_DIR/data

###

export wb_command=/wdata/msmuhammad/workbench/connectome-wb/workbench/bin_linux64/wb_command
$wb_command -volume-to-surface-mapping \
    $STAT_MAP \
    $TEMP_DIR/SUMA/lh.smoothwm.gii \
    lh_stats.func.gii \
    -trilinear

$wb_command -volume-to-surface-mapping \
    $STAT_MAP \
    $TEMP_DIR/SUMA/rh.smoothwm.gii \
    rh_stats.func.gii \
    -trilinear

$wb_command -cifti-create-dense-scalar \
    map.dscalar.nii \
    -left-metric lh_stats.func.gii \
    -right-metric rh_stats.func.gii


