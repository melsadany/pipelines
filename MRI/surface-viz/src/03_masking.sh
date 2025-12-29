#!/bin/bash

language_net="/wdata/msmuhammad/data/neurosynth/language_network_association-test_z_FDR_0.01.nii.gz"

P_DIR=/wdata/msmuhammad/pipelines/MRI/surface-viz
cd $P_DIR/data/masked-input

#### create masks 
# Create binary mask from language network (threshold > 0)
fslmaths "$language_net" -abs -thr 0.0001 -bin language_mask_binary.nii.gz
# Create complement mask and apply
fslmaths language_mask_binary.nii.gz -binv complement_mask_binary.nii.gz



#### apply masks
type=STAT
input="/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/USE_THIS/2E_044/PS-VC/sla/major/sub-2E_044_word_association__word_Tstat.nii.gz"
masked_output=${type}_language-network.nii.gz
complement_output=${type}_language-network-complement.nii.gz
fslmaths "$input" -mas language_mask_binary.nii.gz "$masked_output"
fslmaths "$input" -mas complement_mask_binary.nii.gz "$complement_output"


type=fALFF
input="/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/2E_044/run-3/REST1/L_sub-2E_044_REST1_MNI-reg-fALFF.nii.gz"
masked_output=${type}_language-network.nii.gz
complement_output=${type}_language-network-complement.nii.gz
fslmaths "$input" -mas language_mask_binary.nii.gz "$masked_output"
fslmaths "$input" -mas complement_mask_binary.nii.gz "$complement_output"


##### prep for SUMA

export SUBJECTS_DIR=/wdata/msmuhammad/refs/mni-freesurfer
SUBJECT_ID="2mm"
TEMP_DIR=$SUBJECTS_DIR/$SUBJECT_ID

types=("STAT" "fALFF")
langs=("language-network" "language-network-complement")
hemis=("rh" "lh")
for type in ${types[@]}; do
    for lang in ${langs[@]}; do
    	prefix=${type}_${lang}
    	for hemi in ${hemis[@]}; do
    	    in=${prefix}.nii.gz
    	    out=${prefix}_${hemi}.niml.dset
    	    3dVol2Surf -spec $TEMP_DIR/SUMA/std.141.fsaverage_both.spec \
    	        -surf_A $hemi.white -surf_B $hemi.pial \
    	        -sv $TEMP_DIR/SUMA/fsaverage_SurfVol.nii \
    	        -grid_parent $in \
    	        -map_func max -f_steps 10 -f_index nodes \
    	        -out_niml ${out}
    	    echo $out
    	done
    done
done


suma -spec $TEMP_DIR/SUMA/std.141.fsaverage_both.spec -sv $TEMP_DIR/SUMA/fsaverage_SurfVol.nii

