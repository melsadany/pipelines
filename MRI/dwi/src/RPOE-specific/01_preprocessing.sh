#!/bin/bash
# this script is for preprocessing the dwi


# load all neuroimaging needed stuff
module load stack/2022.2
module load ants
module load FSL 
#module load FreeSurfer
MRTRIXDIR=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/mrtrix3/3.0.4/bin
export PATH=${MRTRIXDIR}:$PATH
export PATH=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07:/Users/msmuhammad/packages/FreeSurfer/7.4.1:$PATH


participant_id=$1
P_DIR=$2
sub_DIR=$3
task_NAME=$4

###### trial run
: << 'trialrun'
participant_id=2E_126
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/dwi
sub_DIR=${P_DIR}/data/derivatives/dwi/${participant_id}
task_NAME="DTI_32_DIR"
trialrun
######


mkdir -p ${sub_DIR}
cd ${sub_DIR}
IMAGE_prefix=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_${task_NAME}
#IMAGE_rev_prefix=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_DTI_Rev_PE
IMAGE_nii=${IMAGE_prefix}.nii.gz
IMAGE_json=${IMAGE_prefix}.json
IMAGE_bval=${IMAGE_prefix}.bval
IMAGE_bvec=${IMAGE_prefix}.bvec

OUT_prefix=sub-${participant_id}_dwi_${task_NAME}
#OUT_rev_prefix=sub-${participant_id}_dwi_DTI_Rev_PE


#########
# STEP1 #
#########
# convert to mrtrix format
mrconvert ${IMAGE_nii} ${sub_DIR}/00_${OUT_prefix}.mif \
    -fslgrad ${IMAGE_bvec} ${IMAGE_bval}

# for the reverse coding
#mrconvert ${IMAGE_rev_prefix}.nii.gz \
#    ${sub_DIR}/00_${OUT_rev_prefix}.mif
#mrconvert ${sub_DIR}/00_${OUT_rev_prefix}.mif \
#    -fslgrad ${IMAGE_rev_prefix}.bvec ${IMAGE_rev_prefix}.bval \
#    ${sub_DIR}/00_02_${OUT_rev_prefix}.mif
    
#mrmath ${sub_DIR}/00_02_${OUT_rev_prefix}.mif mean ${sub_DIR}/00_02_${OUT_rev_prefix}-mean-b0.mif -axis 3


#########
# STEP2 #
#########
# denoise
dwidenoise ${sub_DIR}/00_${OUT_prefix}.mif \
    ${sub_DIR}/01_${OUT_prefix}_denoised.mif \
    -noise ${sub_DIR}/01_${OUT_prefix}_noise.mif

# if you don't think it's denoised enought, add the extent parameter and change it to 7 or other (default is 5)
#   -extent 7


# if there's Gibbs ringing, run this
#mrdegibbs ${sub_DIR}/01_${OUT_prefix}_denoised.mif \
#    ${sub_DIR}/01_${OUT_prefix}_denoised-unringed.mif



#########
# STEP3 #
#########
# combine main and reverse decoding
#dwiextract ${sub_DIR}/01_${OUT_prefix}_denoised.mif - -bzero | mrmath - mean ${sub_DIR}/00_02_${OUT_prefix}-mean-b0.mif -axis 3
#mrcat ${sub_DIR}/00_02_${OUT_prefix}-mean-b0.mif ${sub_DIR}/00_02_${OUT_rev_prefix}-mean-b0.mif -axis 3 ${sub_DIR}/01_${OUT_prefix}-pair.mif



#########
# STEP4 #
#########
# preprocessing: this includes topup and eddy

#dwifslpreproc ${sub_DIR}/01_${OUT_prefix}_denoised.mif ${sub_DIR}/01_${OUT_prefix}_denoised-preproc.mif \
#    -nocleanup \
#    -pe_dir AP \
#    -rpe_pair \
#    -se_epi ${sub_DIR}/01_${OUT_prefix}-pair.mif \
#    -eddy_options " --slm=linear --data_is_shelled" \
#    -nthreads 45

# got rid of this paramater
# -pe_dir AP 
# -json_import ${IMAGE_json} \    

# without rev
dwifslpreproc ${sub_DIR}/01_${OUT_prefix}_denoised.mif ${sub_DIR}/02_${OUT_prefix}_denoised-preproc.mif \
    -nocleanup \
    -pe_dir AP \
    -rpe_none \
    -eddy_options " --slm=linear --data_is_shelled" \
    -nthreads 55



#########
# STEP5 #
#########
# Bias Correction
dwibiascorrect ants \
    ${sub_DIR}/02_${OUT_prefix}_denoised-preproc.mif \
    ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased.mif \
    -bias ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-bias.mif



#########
# STEP6 #
#########
# mask generation
#dwi2mask ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased.mif \
#    ${sub_DIR}/04_${OUT_prefix}_mask.mif

# try other methods in case the first generated mask wasn't good
# convert input to nii
mrconvert ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased.mif \
    ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased.nii

# try using bet
#bet2 ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased.nii \
#    ${sub_DIR}/04_${OUT_prefix}_denoised-preproc-unbiased_masked \
#    -m -f 0.6
#mrconvert ${sub_DIR}/04_${OUT_prefix}_denoised-preproc-unbiased_masked_mask.nii.gz \
#    ${sub_DIR}/04_${OUT_prefix}_mask-bet.mif

# hd-bet
# it only works with 3D images
fslroi ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased.nii \
       ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased_TP0.nii.gz \
       0 1
conda activate ENA
hd-bet -i ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased_TP0.nii -o ${sub_DIR}/04_${OUT_prefix}_hd-bet-brain.nii.gz
slicer ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased_TP0.nii.gz ${sub_DIR}/04_${OUT_prefix}_hd-bet-brain.nii.gz -a hd-bet-mask-overlay

# mask the input
fslmaths ${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased.nii \
         -mas ${sub_DIR}/04_${OUT_prefix}_hd-bet-brain_mask.nii.gz \
         ${sub_DIR}/04_${OUT_prefix}_denoised-preproc-unbiased_masked.nii.gz



mrconvert ${sub_DIR}/04_${OUT_prefix}_hd-bet-brain_mask.nii.gz \
    ${sub_DIR}/04_${OUT_prefix}_mask-hd-bet.mif -force

#########
echo "------------DONE------------"


