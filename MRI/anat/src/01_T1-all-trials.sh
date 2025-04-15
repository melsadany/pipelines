module load stack/2022.2
module load ants
# the following modules/tools are installed locally and not available among argon software list
module load FSL 
module load FreeSurfer 
export PATH=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07:$PATH
conda activate ENA # mainly for HD-bet


## reference data
MNI_W=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm.nii.gz
MNI_B=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm_brain.nii.gz
MNI_BM=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm_brain_mask.nii.gz


## project setup
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/pipelines/MRI/anat
export SUBJECTS_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/pipelines/MRI/anat/data/derivatives/freesurfer
mkdir -p $SUBJECTS_DIR

## project/subject data
participant_id='2E_109'
inputScan=${P_DIR}/data/raw/sub-${participant_id}_anat-T1.nii.gz
sub_DIR=${P_DIR}/data/derivatives/anat/${participant_id}
mkdir -p ${sub_DIR}
cd ${sub_DIR}
outputNameBase=sub-${participant_id}_anat-T1

if [ ! -d $sub_DIR ]; then
  mkdir $sub_DIR
elif [ -f ${sub_DIR}/F_${outputNameBase}_ANTs-brain.nii.gz ]; then
  echo "Preprocessing has already been run"
  exit 0
fi


##################################################################################################
##################################################################################################
## order of processing
# 1. reoriting to RPI orientation
# 2. Denoise
# 3. intensity non-uniformaty correction
# 4. brain extraction
# 5. brain segmentation
# 6: STL file creation
##################################################################################################
##################################################################################################

## Step 1: reorient to RPI
3dresample \
    -orient rpi \
    -overwrite \
    -input $inputScan \
    -prefix A_${outputNameBase}_reorient_RPI.nii.gz
if [ -f A_${outputNameBase}_reorient_RPI.nii.gz ]; then
    echo "------Done with RPI reorientation"
fi


#### FreeSrfer start #1
subj=RPI_start
mkdir -p $SUBJECTS_DIR/$subj/mri/orig
export logfile=$SUBJECTS_DIR/$subj/logfile.log
cd $SUBJECTS_DIR/$subj
mri_convert ${sub_DIR}/A_${outputNameBase}_reorient_RPI.nii.gz $SUBJECTS_DIR/$subj/mri/orig/001.mgz
recon-all -all -subjid $subj -sd $SUBJECTS_DIR -threads 35
####



##################################################################################################
## Step 2: denoising
# not including mask in denoising since it hasn't been generated yet
# Denoising done by ANTs
DenoiseImage \
    -d 3 \
    -n Rician \
    -s 1 -p 1 -r 2 -v 1 \
    -i A_${outputNameBase}_reorient_RPI.nii.gz \
    -o [B_${outputNameBase}_reorient_RPI_denoise.nii.gz,B_${outputNameBase}_reorient_RPI_noise.nii.gz]

: <<'gibbsfix'
# some participants had gibbs ringing, so I tried doing some smoothing hoping it can fix that
fslmaths B_${outputNameBase}_reorient_RPI_denoise.nii.gz \
    -s 0.5 \
    B_${outputNameBase}_reorient_RPI_denoise.nii.gz
gibbsfix

if [ -f B_${outputNameBase}_reorient_RPI_denoise.nii.gz ]; then
    echo "------Done with denoising step"
fi

##################################################################################################
## Step 3: bias correction
# could be done using different tools: ANTs N$ or AFNI

# ANTs intensity correction
N4BiasFieldCorrection \
    -d 3 \
    -i B_${outputNameBase}_reorient_RPI_denoise.nii.gz \
    -r 1 -s 4 \
    -c [50x50x50x50,0.0] \
    -b [200,3] \
    -t [0.15,0.01,200] \
    -o [C_${outputNameBase}_ANTs-bfcorr.nii.gz,C_${outputNameBase}_ANTs-bf.nii.gz]
if [ -f C_${outputNameBase}_ANTs-bfcorr.nii.gz ]; then
    echo "------Done with N4Bias correction"
fi

#### FreeSrfer start #2
subj=N4ANTs_start
mkdir -p $SUBJECTS_DIR/$subj/mri/orig
export logfile=$SUBJECTS_DIR/$subj/logfile.log
cd $SUBJECTS_DIR/$subj
mri_convert ${sub_DIR}/C_${outputNameBase}_ANTs-bfcorr.nii.gz $SUBJECTS_DIR/$subj/mri/orig/001.mgz
recon-all -all -subjid $subj -sd $SUBJECTS_DIR -threads 35
####


# AFNI intensity correction; takes a long time but good output
3dUnifize \
    -prefix C_${outputNameBase}_AFNI-bfcorr.nii.gz \
    -input B_${outputNameBase}_reorient_RPI_denoise.nii.gz \
    -GM \
    -Urad 30
if [ -f C_${outputNameBase}_AFNI-bfcorr.nii.gz ]; then
    echo "------Done with AFNI bias correction"
fi

#### FreeSrfer start #3
subj=UnifizeAFNI_start
mkdir -p $SUBJECTS_DIR/$subj/mri/orig
export logfile=$SUBJECTS_DIR/$subj/logfile.log
cd $SUBJECTS_DIR/$subj
mri_convert ${sub_DIR}/C_${outputNameBase}_AFNI-bfcorr.nii.gz $SUBJECTS_DIR/$subj/mri/orig/001.mgz
recon-all -all -subjid $subj -sd $SUBJECTS_DIR -threads 35
####

##################################################################################################
## Step 4: brain extraction
# this is the main brain mask generation step
# tools include BET from FSL, ANTs, and AFNI
BC=("AFNI" "ANTs")
for BCM in "${BC[@]}"; do
    echo ${BCM}
    
    ## ANTs brain extraction
    antsBrainExtraction.sh \
        -d 3 \
        -a C_${outputNameBase}_${BCM}-bfcorr.nii.gz \
        -e /Dedicated/jmichaelson-wdata/msmuhammad/refs/OASIS/T_template0.nii.gz \
        -m /Dedicated/jmichaelson-wdata/msmuhammad/refs/OASIS/T_template0_BrainCerebellumProbabilityMask.nii.gz \
        -f /Dedicated/jmichaelson-wdata/msmuhammad/refs/OASIS/T_template0_BrainCerebellumRegistrationMask.nii.gz \
        -o D_${outputNameBase}_${BCM}-bfcorr_ants-brain_
    rm -rf D_${outputNameBase}_${BCM}-bfcorr_ants-brain_

    CopyImageHeaderInformation C_${outputNameBase}_${BCM}-bfcorr.nii.gz \
        D_${outputNameBase}_${BCM}-bfcorr_ants-brain_BrainExtractionMask.nii.gz \
        D_${outputNameBase}_${BCM}-bfcorr_ants-brain_BrainExtractionMask.nii.gz 1 1 1
    mv D_${outputNameBase}_${BCM}-bfcorr_ants-brain_BrainExtractionMask.nii.gz \
        D_${outputNameBase}_${BCM}-bfcorr_ants-mask.nii.gz
    mv D_${outputNameBase}_${BCM}-bfcorr_ants-brain_BrainExtractionBrain.nii.gz \
        D_${outputNameBase}_${BCM}-bfcorr_ants-brain.nii.gz
    mv D_${outputNameBase}_${BCM}-bfcorr_ants-brain_BrainExtractionPrior0GenericAffine.mat \
        D_${outputNameBase}_${BCM}-bfcorr_ants-Affine.mat

    
    ## BET 
    bet C_${outputNameBase}_${BCM}-bfcorr.nii.gz \
        D_${outputNameBase}_${BCM}-bfcorr_BET.nii.gz \
        -R -m -f 0.4

    
    ## AFNI
    3dSkullStrip \
        -input C_${outputNameBase}_${BCM}-bfcorr.nii.gz \
        -prefix D_${outputNameBase}_${BCM}-bfcorr_AFNI-mask
    3dcalc -a D_${outputNameBase}_${BCM}-bfcorr_AFNI-mask+orig'[0]' \
        -expr 'a' -prefix D_${outputNameBase}_${BCM}-bfcorr_AFNI-mask.nii.gz
    rm D_${outputNameBase}_${BCM}-bfcorr_AFNI-mask+orig*
    
    
    ## HD_BET
    hd-bet -i C_${outputNameBase}_${BCM}-bfcorr.nii.gz -o D_${outputNameBase}_${BCM}-bfcorr_HD-bet.nii.gz
    
    #### FreeSrfer start #3
    subj=HDBEt_${BCM}_start
    mkdir -p $SUBJECTS_DIR/$subj/mri/orig
    export logfile=$SUBJECTS_DIR/$subj/logfile.log
    cd $SUBJECTS_DIR/$subj
    mri_convert ${sub_DIR}/D_${outputNameBase}_${BCM}-bfcorr_HD-bet.nii.gz $SUBJECTS_DIR/$subj/mri/orig/001.mgz
    recon-all -all -subjid $subj -sd $SUBJECTS_DIR -threads 35
    ####

done

##################################################################################################
## Step 5: brain segmentation
# once you get a good brain mask, extract roi segmentation on it
# segmentation could be on the tissue type or region
# the best brain mask method is HD-BET 

BC=("AFNI" "ANTs")
for BCM in "${BC[@]}"; do
    echo ${BCM}
    
    # FAST
    fast -v D_${outputNameBase}_${BCM}-bfcorr_HD-bet.nii.gz
    
    
    ## mri_synthseg
    # it needs Freesurfer with a new version. at least 7.4.1
    export FREESURFER_HOME=/Users/msmuhammad/packages/FreeSurfer/7.4.1
    export FSFAST_HOME=/Users/msmuhammad/packages/FreeSurfer/7.4.1/fsfast
    export SUBJECTS_DIR=/Users/msmuhammad/packages/FreeSurfer/7.4.1/subjects
    export MNI_DIR=/Users/msmuhammad/packages/FreeSurfer/7.4.1/mni
    export PATH=${FREESURFER_HOME}/bin:$PATH
    
    # using the bias-field corrected
    mkdir -p synthseg/${BCM}
    mri_synthseg \
        --i C_${outputNameBase}_${BCM}-bfcorr.nii.gz \
        --o synthseg/${BCM}/ \
        --vol synthseg/${BCM}/vols.csv \
        --threads 4
    # using the HD-BET brain-extracted
    mkdir -p synthseg/${BCM}_HD-BET
    mri_synthseg \
        --i D_${outputNameBase}_${BCM}-bfcorr_HD-bet.nii.gz \
        --o synthseg/${BCM}_HD-BET/ \
        --vol synthseg/${BCM}_HD-BET/vols.csv \
        --threads 4
    
done
##################################################################################################
## Step 6: STL file generation for 3d-printing
# choose one good output from syntheg
BC=("AFNI" "ANTs")
for BCM in "${BC[@]}"; do
    echo ${BCM}
    
    ## extract cortex regions
    # bias-field corrected
    mri_binarize --i synthseg/${BCM}/D_${outputNameBase}_${BCM}-bfcorr_synthseg.nii.gz \
        --match 2 3 31 41 42 63 72 77 51 52 13 12 43 50 4 11 26 58 49 10 17 18 53 54 44 5 80 14 30 62 \
        --o synthseg/${BCM}/cortex-mask.nii.gz
    # HD-BET brains corrected
    mri_binarize --i synthseg/${BCM}_HD-BET/D_${outputNameBase}_${BCM}-bfcorr_HD-bet_synthseg.nii.gz \
        --match 2 3 31 41 42 63 72 77 51 52 13 12 43 50 4 11 26 58 49 10 17 18 53 54 44 5 80 14 30 62 \
        --o synthseg/${BCM}_HD-BET/cortex-mask.nii.gz
        
        
    ## mesh generation    
    flirt \
        -in D_${outputNameBase}_${BCM}-bfcorr_HD-bet.nii.gz \
        -ref synthseg/${BCM}_HD-BET/cortex-mask.nii.gz \
        -out synthseg/${BCM}_HD-BET/resampled_brain.nii.gz \
        -applyxfm -usesqform
    fslmaths synthseg/${BCM}_HD-BET/resampled_brain.nii.gz \
        -mas synthseg/${BCM}_HD-BET/cortex-mask.nii.gz \
        synthseg/${BCM}_HD-BET/brain.nii.gz
    # needs the binary of nii2mesh
    /Dedicated/jmichaelson-wdata/msmuhammad/workbench/nii2mesh/src/nii2mesh \
        synthseg/${BCM}_HD-BET/brain.nii.gz \
        -r 0.5 \
        synthseg/${BCM}_HD-BET/final.stl

done



################################################################################
################################################################################
### Muhammad's approach
################################################################################
################################################################################
## Step 3: bias correction
# use a generated HD-BET mask to do the bias correction

mkdir -p iterative_BC
IN_IMG=B_${outputNameBase}_reorient_RPI_denoise.nii.gz
OUT_MASK=iterative_BC/C_${outputNameBase}_ANTs-bfcorr_HD-bet.nii.gz
OUT_IMG_PRE=iterative_BC/D_${outputNameBase}_ANTs
i=1
while [ $i -le 3 ]; do
    echo ${i}
    
    # generate mask
    hd-bet -i ${IN_IMG} \
        -o ${OUT_MASK}
    # ANTs intensity correction
    N4BiasFieldCorrection \
        -d 3 \
        -i ${IN_IMG} \
        -x ${OUT_MASK} \
        -r 1 -s 4 \
        -c [50x50x50x50,0.0] \
        -b [200,3] \
        -t [0.15,0.01,200] \
        -o [${OUT_IMG_PRE}-bfcorr.nii.gz,${OUT_IMG_PRE}-bf.nii.gz]
    
    IN_IMG=${OUT_IMG_PRE}-bfcorr.nii.gz
    i=$(($i + 1 ))
done
# last
hd-bet -i ${OUT_IMG_PRE}-bfcorr.nii.gz -o ${OUT_MASK}

#### FreeSrfer start #6
subj=N4ANTs_HDBET_iterative_start
mkdir -p $SUBJECTS_DIR/$subj/mri/orig
export logfile=$SUBJECTS_DIR/$subj/logfile.log
cd $SUBJECTS_DIR/$subj
mri_convert ${sub_DIR}/iterative_BC/C_${outputNameBase}_ANTs-bfcorr_HD-bet.nii.gz $SUBJECTS_DIR/$subj/mri/orig/001.mgz
recon-all -all -subjid $subj -sd $SUBJECTS_DIR -threads 35
####

################################################################################
## Step 5: brain segmentation

mkdir -p iterative_BC/synthseg
## mri_synthseg
export FREESURFER_HOME=/Users/msmuhammad/packages/FreeSurfer/7.4.1
export FSFAST_HOME=/Users/msmuhammad/packages/FreeSurfer/7.4.1/fsfast
export SUBJECTS_DIR=/Users/msmuhammad/packages/FreeSurfer/7.4.1/subjects
export MNI_DIR=/Users/msmuhammad/packages/FreeSurfer/7.4.1/mni
export PATH=${FREESURFER_HOME}/bin:$PATH
mri_synthseg \
    --i iterative_BC/C_${outputNameBase}_ANTs-bfcorr_HD-bet.nii.gz \
    --o iterative_BC/synthseg \
    --vol iterative_BC/synthseg/vols.csv \
    --threads 4
################################################################################
## Step 6: MNI-registration
### ANTs SyN registration
antsRegistrationSyN.sh \
     -d 3 \
     -f ${MNI_W} \
     -m iterative_BC/D_${outputNameBase}_ANTs-bfcorr.nii.gz \
     -o ${sub_DIR}/iterative_BC/E_${outputNameBase}_T1-to-MNI-SyN_ \
     -x ${MNI_BM} \
     -n 4
cp ${sub_DIR}/iterative_BC/E_${outputNameBase}_T1-to-MNI-SyN_Warped.nii.gz ${sub_DIR}/iterative_BC/E_${outputNameBase}_T1-to-MNI-SyN.nii.gz

# make a PNG of the registration result and save it
slicer ${MNI_W} ${sub_DIR}/iterative_BC/E_${outputNameBase}_T1-to-MNI-SyN.nii.gz \
    -a ${sub_DIR}/iterative_BC/E_${outputNameBase}_T1-to-MNI-SyN.png

################################################################################
################################################################################
################################################################################
################################################################################






