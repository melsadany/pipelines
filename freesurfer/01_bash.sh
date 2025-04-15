

conda activate ENA
module load FreeSurfer
module load FSL


export SUBJECTS_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni-freesurfer

subj=2mm

mkdir -p $SUBJECTS_DIR/$subj/mri/orig
export logfile=$SUBJECTS_DIR/$subj/logfile.log
cd $SUBJECTS_DIR/$subj

mri_convert /Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_$subj.nii.gz $SUBJECTS_DIR/$subj/mri/orig/001.mgz

recon-all -all -subjid $subj -sd $SUBJECTS_DIR -threads 127




## freesurfer does the surface for cortical parts only
# do it for the subcortical parts here
mkdir -p $SUBJECTS_DIR/$subj/surf2
mri_convert $SUBJECTS_DIR/$subj/mri/aseg.mgz $SUBJECTS_DIR/$subj/surf2/start.nii

## right hemisphere
mri_binarize --i $SUBJECTS_DIR/$subj/surf2/start.nii \
    --match 1 2 3 4 5 10 11 12 13 17 18 19 20 26 27 28 31 \
    --o $SUBJECTS_DIR/$subj/surf2/lh-bin.nii
fslmaths $SUBJECTS_DIR/$subj/surf2/start.nii \
    -mul $SUBJECTS_DIR/$subj/surf2/lh-bin.nii \
    $SUBJECTS_DIR/$subj/surf2/lh.nii.gz
cp $SUBJECTS_DIR/$subj/surf2/lh.nii.gz $SUBJECTS_DIR/$subj/surf2/lh_tmp.nii.gz
gunzip -f $SUBJECTS_DIR/$subj/surf2/lh_tmp.nii.gz
for i in 1 2 3 4 5 10 11 12 13 17 18 19 20 26 27 28 31; do 
    mri_pretess $SUBJECTS_DIR/$subj/surf2/lh_tmp.nii \
    $i $SUBJECTS_DIR/$subj/mri/norm.mgz\
    $SUBJECTS_DIR/$subj/surf2/lh_tmp.nii 
done
fslmaths $SUBJECTS_DIR/$subj/surf2/lh_tmp.nii -bin \
    $SUBJECTS_DIR/$subj/surf2/lh_bin.nii
mri_tessellate $SUBJECTS_DIR/$subj/surf2/lh_bin.nii.gz 1 \
    $SUBJECTS_DIR/$subj/surf2/lh.surf
mris_convert $SUBJECTS_DIR/$subj/surf2/lh.surf $SUBJECTS_DIR/$subj/surf2/lh.stl
mris_convert $SUBJECTS_DIR/$subj/surf2/lh.surf $SUBJECTS_DIR/$subj/surf2/lh.surf.gii

## left hemispher
mri_binarize --i $SUBJECTS_DIR/$subj/surf2/start.nii \
    --match 40 41 42 43 44 49 50 51 52 53 54 55 56 58 59 60 63 \
    --o $SUBJECTS_DIR/$subj/surf2/rh-bin.nii
fslmaths $SUBJECTS_DIR/$subj/surf2/start.nii \
    -mul $SUBJECTS_DIR/$subj/surf2/rh-bin.nii \
    $SUBJECTS_DIR/$subj/surf2/rh.nii.gz
cp $SUBJECTS_DIR/$subj/surf2/rh.nii.gz $SUBJECTS_DIR/$subj/surf2/rh_tmp.nii.gz
gunzip -f $SUBJECTS_DIR/$subj/surf2/rh_tmp.nii.gz
for i in 40 41 42 43 44 49 50 51 52 53 54 55 56 58 59 60 63; do 
    mri_pretess $SUBJECTS_DIR/$subj/surf2/rh_tmp.nii \
    $i $SUBJECTS_DIR/$subj/mri/norm.mgz\
    $SUBJECTS_DIR/$subj/surf2/rh_tmp.nii 
done
fslmaths $SUBJECTS_DIR/$subj/surf2/rh_tmp.nii -bin \
    $SUBJECTS_DIR/$subj/surf2/rh_bin.nii
mri_tessellate $SUBJECTS_DIR/$subj/surf2/rh_bin.nii.gz 1 \
    $SUBJECTS_DIR/$subj/surf2/rh.surf
mris_convert $SUBJECTS_DIR/$subj/surf2/rh.surf $SUBJECTS_DIR/$subj/surf2/rh.stl
mris_convert $SUBJECTS_DIR/$subj/surf2/rh.surf $SUBJECTS_DIR/$subj/surf2/rh.surf.gii

## subcortical only
mri_binarize --i $SUBJECTS_DIR/$subj/surf2/start.nii \
    --match 4 5 10 11 12 13 17 18 19 20 26 27 28 31 43 44 49 50 51 52 53 54 55 56 58 59 60 63 \
    --o $SUBJECTS_DIR/$subj/surf2/sc-bin.nii
fslmaths $SUBJECTS_DIR/$subj/surf2/start.nii \
    -mul $SUBJECTS_DIR/$subj/surf2/sc-bin.nii \
    $SUBJECTS_DIR/$subj/surf2/sc.nii.gz
cp $SUBJECTS_DIR/$subj/surf2/sc.nii.gz $SUBJECTS_DIR/$subj/surf2/sc_tmp.nii.gz
gunzip -f $SUBJECTS_DIR/$subj/surf2/sc_tmp.nii.gz
for i in 4 5 10 11 12 13 17 18 19 20 26 27 28 31 43 44 49 50 51 52 53 54 55 56 58 59 60 63; do 
    mri_pretess $SUBJECTS_DIR/$subj/surf2/sc_tmp.nii \
    $i $SUBJECTS_DIR/$subj/mri/norm.mgz\
    $SUBJECTS_DIR/$subj/surf2/sc_tmp.nii 
done
fslmaths $SUBJECTS_DIR/$subj/surf2/sc_tmp.nii -bin \
    $SUBJECTS_DIR/$subj/surf2/sc_bin.nii
mri_tessellate $SUBJECTS_DIR/$subj/surf2/sc_bin.nii.gz 1 \
    $SUBJECTS_DIR/$subj/surf2/sc.surf
mris_convert $SUBJECTS_DIR/$subj/surf2/sc.surf $SUBJECTS_DIR/$subj/surf2/sc.stl
mris_convert $SUBJECTS_DIR/$subj/surf2/sc.surf $SUBJECTS_DIR/$subj/surf2/sc.surf.gii







