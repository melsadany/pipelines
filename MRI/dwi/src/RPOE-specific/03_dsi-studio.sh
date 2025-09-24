#!/bin/bash

conda activate ENA

participant_id=$1

###### trial run
: << 'trialrun'
participant_id=2E_126
trialrun
######
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/dwi
sub_DIR=${P_DIR}/data/derivatives/dwi/${participant_id}
task_NAME="DTI_32_DIR"

cd ${sub_DIR}
IMAGE_prefix=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_${task_NAME}
OUT_prefix=sub-${participant_id}_dwi_${task_NAME}

MNI_W=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm.nii.gz



##########
# STEP13 #
##########
# tractography using DSI studio
# all info about the CLI usage and options are here
# http://dsi-studio.labsolver.org/Manual/command-line-for-dsi-studio


DSI_bin=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/dsi-studio/dsi_studio

mkdir -p dsi-studio
cd dsi-studio
prefix=sub-${participant_id}

### A
# Create the src file 
${DSI_bin} \
    --action=src \
    --source=${sub_DIR}/03_${OUT_prefix}_denoised-preproc-unbiased.nii \
    --bvec=${IMAGE_prefix}.bvec \
    --bval=${IMAGE_prefix}.bval \
    --output=${prefix}.src.gz

### B
# QC src file 
${DSI_bin} \
    --action=qc \
    --source=dsi-studio \
    --output ${prefix}


### C
# Image reconstruction 
# method 7 is QSDR reconstruction
# template 0 is ICBM152
${DSI_bin} \
    --action=rec \
    --source=${prefix}.src.gz \
    --mask=${sub_DIR}/04_${OUT_prefix}_hd-bet-brain_mask.nii.gz \
    --method=7 \
    --template=0 \
    --output=${prefix}



#FIB_FILE=($(ls $prefix*fib.gz))
#SRC_FILE=($(ls ${prefix}*src.gz))
FIB_FILE=$(find . -maxdepth 1 -name "${prefix}*fib.gz" -print -quit)
SRC_FILE=$(find . -maxdepth 1 -name "${prefix}*src.gz" -print -quit)


### D
# Fiber tracking
# basic tracking with cluster recognition
${DSI_bin} \
    --action=trk \
    --source=${FIB_FILE} \
    --output=${prefix}.tt.gz \
    --fa_threshold=0.2 \
    --export=stat,tdi \
    --recognize=cluster_info


# connectivity using different atlases
atlases=("FreeSurferDKT_Cortical" "HCP-MMP" "HCP842_tractography")
# every folder will have .tt.gz per tract/roi
# the folder will also have a .tt.gz for the subject
for atlas in "${atlases[@]}"; do
    mkdir -p conn_${atlas}
    ${DSI_bin} \
        --action=trk \
        --source=${FIB_FILE} \
        --fiber_count=1000000 \
        --output=conn_${atlas} \
        --connectivity=${atlas} \
        --connectivity_value=count,qa,trk
done



# automatic fiber tracking
${DSI_bin} \
    --action=atk \
    --source=${FIB_FILE} \
    --output=${prefix}.atk \
    --template=0
    
    
    

### E
# visualization
#${DSI_bin} \
#    --action=vis \
#    --source=${FIB_FILE} \
#    --tract=${prefix}.tt.gz \
#    --cmd="Full+save_h3view_image,${prefix}.png"



### F
# analysis
# for the whole file
${DSI_bin} \
    --action=ana \
    --source=${FIB_FILE} \
    --atlas=FreeSurferDKT_Cortical,HCP-MMP,HCP842_tractography \
    --output=${prefix}.statistics.txt




#########
echo "------------DONE------------"

    


