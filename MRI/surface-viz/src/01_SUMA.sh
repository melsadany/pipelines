#!/bin/bash


## needed to prepare surface files for the MNI freesurfer output
export SUBJECTS_DIR=/wdata/msmuhammad/refs/mni-freesurfer
SUBJECT_ID="2mm"
TEMP_DIR=$SUBJECTS_DIR/$SUBJECT_ID

# Generate SUMA files if missing, then launch SUMA
if [ ! -f $TEMP_DIR/SUMA/fsaverage_both.spec ]; then
    cd $TEMP_DIR
    @SUMA_Make_Spec_FS -sid $SUBJECT_ID -NIFTI
fi


## get your statistic map
STAT_MAP="/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/USE_THIS/2E_044/PS-VC/sla/major/sub-2E_044_word_association__word_Tstat.nii.gz"

P_DIR=/wdata/msmuhammad/pipelines/MRI/surface-viz
cd $P_DIR/data

###
@SUMA_Make_Spec_FS -sid ${SUBJECT_ID} -fspath ${TEMP_DIR}/surf


3dVol2Surf \
  -spec $TEMP_DIR/SUMA/std.141.fsaverage_both.spec \
  -surf_A lh.white \
  -surf_B lh.pial \
  -sv $TEMP_DIR/SUMA/fsaverage_SurfVol.nii \
  -grid_parent statmap_resampled.nii.gz \
  -map_func max \
  -f_steps 10 -f_index nodes \
  -out_niml statmap_lh.niml.dset

3dVol2Surf \
  -spec $TEMP_DIR/SUMA/std.141.fsaverage_both.spec \
  -surf_A rh.white \
  -surf_B rh.pial \
  -sv $TEMP_DIR/SUMA/fsaverage_SurfVol.nii \
  -grid_parent statmap_resampled.nii.gz \
  -map_func max \
  -f_steps 10 -f_index nodes \
  -out_niml statmap_rh.niml.dset
  
suma -spec $TEMP_DIR/SUMA/std.141.fsaverage_both.spec \
  -sv $TEMP_DIR/SUMA/fsaverage_SurfVol.nii




## convert lang network
lang_MAP="/wdata/msmuhammad/data/neurosynth/language_network_association-test_z_FDR_0.01.nii.gz"
3dVol2Surf \
  -spec $TEMP_DIR/SUMA/std.141.fsaverage_both.spec \
  -surf_A lh.white \
  -surf_B lh.pial \
  -sv $TEMP_DIR/SUMA/fsaverage_SurfVol.nii \
  -grid_parent $lang_MAP \
  -map_func max \
  -f_steps 10 -f_index nodes \
  -out_niml langmap_lh.niml.dset

3dVol2Surf \
  -spec $TEMP_DIR/SUMA/std.141.fsaverage_both.spec \
  -surf_A rh.white \
  -surf_B rh.pial \
  -sv $TEMP_DIR/SUMA/fsaverage_SurfVol.nii \
  -grid_parent $lang_MAP \
  -map_func max \
  -f_steps 10 -f_index nodes \
  -out_niml langmap_rh.niml.dset





