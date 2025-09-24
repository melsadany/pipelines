#!/bin/bash

module load stack/2022.2
module load ants
module load FSL 
#module load FreeSurfer
MRTRIXDIR=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/mrtrix3/3.0.4/bin
export PATH=${MRTRIXDIR}:$PATH
export PATH=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07:/Users/msmuhammad/packages/FreeSurfer/7.4.1:$PATH

participant_id=$1

###### trial run
: << 'trialrun'
participant_id=2E_124
trialrun
######
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/dwi
sub_DIR=${P_DIR}/data/derivatives/dwi/${participant_id}
task_NAME="DTI_32_DIR"

mkdir -p ${sub_DIR}
cd ${sub_DIR}
IMAGE_prefix=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_${task_NAME}
OUT_prefix=sub-${participant_id}_dwi_${task_NAME}

MNI_W=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm.nii.gz

##########
# STEP8 #
##########
# generate tensor data
dtifit \
    --data=${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased.nii \
    --mask=${sub_DIR}/04_${OUT_prefix}_hd-bet-brain_mask.nii.gz \
    --bvecs=${IMAGE_prefix}.bvec \
    --bvals=${IMAGE_prefix}.bval \
    --out=07_${OUT_prefix}-dtifit \
    --save_tensor
    


##########
# STEP9 #
##########
# register FA to T1/T2
T2_Scan=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/anat/${participant_id}/run-2/T2/F_sub-${participant_id}_anat-T2_T1-rigid-reg.nii.gz
T1_Scan=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/anat/${participant_id}/run-2/E_sub-${participant_id}_anat-T1_ANTs-bfcorr.nii.gz
T1_BM=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/anat/${participant_id}/run-2/F_sub-${participant_id}_anat-T1_hd-bet-brain_mask.nii.gz

antsRegistration \
    --dimensionality 3 \
    --output 08_${OUT_prefix}-dtifit_FA_T1-rigid-reg \
    --initial-moving-transform [${T1_Scan},07_${OUT_prefix}-dtifit_FA.nii.gz,1] \
    --transform Rigid[0.1] \
    --metric Mattes[${T1_Scan},07_${OUT_prefix}-dtifit_FA.nii.gz,1,32,Regular,0.25] \
    --convergence [2000x2000x2000x2000x2000,1e-6,10] \
    --smoothing-sigmas 4x3x2x1x0vox \
    --shrink-factors 8x8x4x2x1 \
    --use-histogram-matching 1 \
    --verbose 1 \
    --random-seed 13983981 \
    --winsorize-image-intensities [0.005,0.995]

# apply transformation matrix
antsApplyTransforms \
    -d 3 \
    -n BSpline[3] \
    -i 07_${OUT_prefix}-dtifit_FA.nii.gz \
    -o 08_${OUT_prefix}-dtifit_FA_T1-rigid-reg.nii.gz \
    -t 08_${OUT_prefix}-dtifit_FA_T1-rigid-reg0GenericAffine.mat \
    -r ${T1_Scan}


slicer ${T1_Scan} 08_${OUT_prefix}-dtifit_FA_T1-rigid-reg.nii.gz -a FA-T1-reg-overlay2


##########
# STEP10 #
##########
# Calculate RD
# Radial Diffusivity (RD) = (lambda2 + lambda3)/2
fslmaths 07_${OUT_prefix}-dtifit_L2.nii.gz \
  -add 07_${OUT_prefix}-dtifit_L3.nii.gz \
  -div 2 07_${OUT_prefix}-dtifit_RD.nii.gz



##########
# STEP11 #
##########
# apply transformation to tensor data to native T1 space
mkdir -p T1-reg_tensor

# copy FA
cp 08_${OUT_prefix}-dtifit_FA_T1-rigid-reg.nii.gz T1-reg_tensor/FA_T1-rigid-reg.nii.gz

# transform the rest
modes=("MD" "RD" "L1" "L2" "L3" "V1" "V2" "V3")
for mode in "${modes[@]}"; do
    echo "processing ${mode}"
    
    antsApplyTransforms \
        -d 3 \
        -n BSpline[3] \
        -i 07_${OUT_prefix}-dtifit_${mode}.nii.gz \
        -o T1-reg_tensor/${mode}_T1-rigid-reg.nii.gz \
        -t 08_${OUT_prefix}-dtifit_FA_T1-rigid-reg0GenericAffine.mat \
        -r ${T1_Scan}
done



##########
# STEP11 #
##########
# apply transformation to tensor data to atlas space



##########
# STEP12 #
##########
# get ROIs summary statistics from these different modes
# use the transformation matrix from the T1w to MNI
# apply this using three different atlases:
# 	Schaefer
#	JHU
#	labeled MNI/Desikan


# transformation files from SyN
transform_prefix=${P_DIR}/../mri/data/derivatives/anat-registration/${participant_id}/sub-${participant_id}-T1-to-MNI_
warp=${transform_prefix}1Warp.nii.gz
affine=${transform_prefix}0GenericAffine.mat

# define your atlases with labeled ROIs here
ROI_mask_MNI=""
ROI_mask_SCHF=/Dedicated/jmichaelson-wdata/msmuhammad/refs/Schaefer2018/Parcellations/MNI/Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz
ROI_mask_JHU=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/atlases/JHU/JHU-ICBM-labels-0.5mm_ME.nii.gz


# just do JHU for now
#atlases=("MNI" "SCHF" "JHU")
atlases=("JHU")
mkdir -p labeled
for atlas in "${atlases[@]}"; do
    mkdir -p labeled/${atlas}
    if [ atlas=="MNI" ]; then
    	ROI_mask=
    fi
    if [ atlas=="SCHF" ]; then
    	ROI_mask=/Dedicated/jmichaelson-wdata/msmuhammad/refs/Schaefer2018/Parcellations/MNI/Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz
    fi    
    if [ atlas=="JHU" ]; then
    	ROI_mask=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/atlases/JHU/JHU-ICBM-labels-0.5mm_ME.nii.gz
    fi
    ROI_out=labeled/${atlas}/${atlas}-in-T1.nii.gz
    
    
    # Apply inverse transformation to bring ROI mask to native T1w space
    antsApplyTransforms -d 3 \
        -i $ROI_mask \
        -r $T1_Scan \
        -o $ROI_out \
        -t [${affine},1] \
        -t ${warp}
    
    
    # Extract the list of ROI labels
    ROI_labels=$(fslstats $ROI_out -R | awk '{for (i=int($1); i<=int($2); i++) if (i != 0) print i}')

    
    modes=("FA" "MD" "RD" "L1" "L2" "L3" "V1" "V2" "V3")
    for mode in "${modes[@]}"; do
    	output_csv=labeled/${atlas}/${mode}-summary-stats.csv
    	echo "Participant,ROI,Max,Min,Mean,Std" > "$output_csv"
    	# Loop over each ROI in the mask
	for ROI in $ROI_labels; do
	    echo "  Analyzing ROI: $ROI"

            # Create a binary mask for the current ROI
            ROI_mask_current=labeled/${atlas}/ROI_${participant_id}_${ROI}_mask.nii.gz
            fslmaths $ROI_out -thr $ROI -uthr $ROI -bin $ROI_mask_current

            i_image=T1-reg_tensor/${mode}_T1-rigid-reg.nii.gz
            # Extract statistics for the current ROI
            max_val=$(fslstats "$i_image" -k "$ROI_mask_current" -R | awk '{print $2}')
            min_val=$(fslstats "$i_image" -k "$ROI_mask_current" -R | awk '{print $1}')
            mean_val=$(fslstats "$i_image" -k "$ROI_mask_current" -M)
            std_val=$(fslstats "$i_image" -k "$ROI_mask_current" -S)
    
            # Append the results to the CSV
            echo "$participant_id,$ROI,$max_val,$min_val,$mean_val,$std_val" >> "$output_csv"
    
            # Remove temporary ROI mask
            rm "$ROI_mask_current"
        done
    done
done




#########
echo "------------DONE------------"

    






