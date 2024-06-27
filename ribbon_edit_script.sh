#!/bin/bash

#################
## Help
#################
Help ()
{
builtin echo "
AUTHOR: Benoît Verreman

LAST UPDATE: 2024-06-27

DESCRIPTION: 
Use ribbon and subcortical NIFTI files to recompute pial surface,
based on previously created <subjid> folder using Freesurfer 7.4.1
Reuse Freesurfer 7.4.1 functions, but not always in the same order as scripts/recon-all.cmd (when output of the function is not immediately used after)
Create a log 'report.sh'.
Create a folder 'outputs' with all the output files.

PREREQUISITES:
*Freesurfer variables:
Export SUBJECTS_DIR and FREESURFER_HOME correctly

*Python script:
Test if you have access to python: 'which python'
Install two python libraries: 'pip install nibabel scipy'
OR use conda environment:
conda create --name env_ribbon_edit_script
conda activate env_ribbon_edit_script
conda install nibabel scipy -c conda-forge
conda list | grep -E 'nibabel|scipy'

*Freesurfer output folder:
If you want to launch Freesurfer 7.4.1 recon-all pipeline using the script:
	Add argument -i
Else:
	Launch Freesurfer 7.4.1 command before using script: 
$ recon-all -s <subjid> -i <subject_image> -autorecon1 -autorecon2 -hires -parallel -openmp 4 -expert expert_file.txt
Prepare ribbon and subcortical NIFTI files (for step 0)
Put the script inside <subjid> folder

*Specific labels and statical values:
Modify the value of some constants in the script if needed: LABELS_SUBCORTICAL, LABEL_RIBBON_WM_LH, LABEL_RIBBON_WM_RH, PIAL_BORDER_LOW

EXAMPLES:
$ bash ribbon_edit_script.sh -i 133019_T1w_acpc_dc_restore.nii.gz -s 133019 -b 133019_ribbon.nii.gz -c 133019_subcortical.nii.gz

$ bash ribbon_edit_script.sh -s <subjid> -t 8 -r #start from pial computation step (8), right hemisphere only (r)

PARAMETERS:

HELP
-h: Print this string, and exit

INPUT FILES
-i: Relative or absolute path to T1w image file
-s: Relative or absolute path to subjid folder (Necessary)
-b: Relative or absolute path to ribbon file
-c: Relative or absolute path to subcortical file

TAG
-t 0: (ribbons) Start with resizing RIBBON_EDIT and SUBCORTICAL
-t 1: (bmask) Start with BRAIN_MASK
-t 2: (maskT1) Start with T1_MASKED
-t 3: (brain.finalsurfs) Start with skull-stripping up to BRAIN_FINALSURFS
-t 4: (wm-bmask) Start the creation of WM_BMASK based on RIBBON_EDIT
-t 5: (wm) Start from computing WM based on WM_BMASK
-t 6: (orig) Start from computing orig surface based on wm from RIBBON_EDIT
-t 7: (brain.finalsurfs-edit) edit brain.finalsurfs with GM from RIBBON_EDIT
-t 8: (stats) Start from computing stats
-t 9: (pial) Start from computing pial surface
-t 10: (smooth) Start from smoothing pial surface
-t 11: (aseg+aparc) Compute stats and other files

VALUES
-p: Give value of PIAL_BORDER_LOW

HEMI
-r: Compute only right hemisphere surface
-l: Compute only left hemisphere surface

RESET
-d: Reset outputs folder and report.sh script

TROUBLESHOOTS:
-Missing argument: check if you put the necessary option flags, and for each flag, if it needs an argument or not.
"
}

#################
## Default global variables
#################
current_date_time=$(date)
TAG=-1 # Start from beginning
HEMI=-1 # Both hemispheres
FS=0 # Default: No recon-all
OUTPUT_FOLDER="outputs"
LABELS_SUBCORTICAL="5 15 29 30 32 31"
declare -a H=("lh" "rh") #Left then Right hemispheres
declare -a LABEL_RIBBON_WM=("2" "41")
declare -a LABEL_RIBBON_GM=("3" "42")

CHANGE_AUTODET=1 # Default: Change autodet with parameters bellow or given as parameter
PIAL_BORDER_LOW=5

#################
## Manage flags
#################
unset -v IMAGE
unset -v SUBJID
unset -v RIBBON
unset -v SUBCORTICAL

#If a character is followed by :, then it needs an argument
VALID_ARGS="i:s:b:c:t:p:hlrd"

while getopts ${VALID_ARGS} opt; do
  case ${opt} in
    i)
        IMAGE=${OPTARG}
        FS=1
        ;;
    s)
        SUBJID=${OPTARG}
        ;;
    b)
        RIBBON=${OPTARG}
        ;;
    c)
        SUBCORTICAL=${OPTARG}
        ;;     
    t)
	TAG=${OPTARG}
	;;
    p)
	PIAL_BORDER_LOW=${OPTARG}
	;;
    h)
	Help
	exit 1
	;;
    l)
	HEMI=0 #left hemisphere only
	;;
    r)
	HEMI=1 #right hemisphere only
	;;
    d)
	Delete #Delete report.sh and $OUTPUT_FOLDER
	;;
    :)
      	echo "Option -${OPTARG} requires an argument."
      	exit 1
      	;;
    ?)
      	echo "Invalid option: -${OPTARG}."
      	exit 1
      	;;
  esac
done

# Test if user provided SUBJID
: ${SUBJID:?Missing argument -s}

#################
## Remove slaches from $SUBJECTS_DIR and $SUBJID
#################
export var="${SUBJECTS_DIR: -1}"
if [[ "$var" == "/" ]]; then
export SUBJECTS_DIR="${SUBJECTS_DIR:0:-1}"
fi

export var="${SUBJID: -1}"
if [[ "$var" == "/" ]]; then
export SUBJID="${SUBJID:0:-1}"
fi

export var="${SUBJID:0:1}"
if [[ "$var" == "/" ]]; then
export SUBJID="${SUBJID:1}"
fi

#################
## Set default permission of working directory to a+rwx
#################
#cd $SUBJECTS_DIR
#umask 0000

#################
## Function to print both on terminal and on script report.sh
#################
Echo ()
{
    builtin echo "$@" | tee -a $SUBJECTS_DIR/report.sh
}

#################
## Function to reset report.sh, $OUTPUT_FOLDER and mri_convert_correction_by_translation.py
#################
CreateScripts()
{
if [ ! -f "$SUBJECTS_DIR/report.sh" ]
then
Echo "#!/bin/bash"
fi

script_nifti_padding
script_brain-finalsurfs-edit
script_expert_file

}

CreateFolders()
{
if [ ! -d "$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER" ]
then
	cmd "Create $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER" \
	"mkdir $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER;
	mkdir $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/scripts;
	mkdir $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf;
	mkdir $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri;
	mkdir $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/transforms;
	mkdir $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/label;
	mkdir $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/stats;"
fi
}

Delete()
{
cmd "Reset $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER" \
"rm -r $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER;"

rm $SUBJECTS_DIR/report.sh
}

#################
## Create a python script "nifti_padding.py"
#################
script_nifti_padding()
{
if [ ! -f "$SUBJECTS_DIR/nifti_padding.py" ]
then

cat > $SUBJECTS_DIR/nifti_padding.py <<EOF
import os
import nibabel as nib
import nibabel.processing #Used in nib.processing.conform
import scipy.ndimage #Used in nib.processing.conform
import sys #To add arguments

# SUBJID directory
img_in = sys.argv[1]
img_out = sys.argv[2]

#Load image to be treated
if not os.path.isfile(img_in):
    raise FileNotFoundError("The following path doesn't exist: " + img_in)
else:
    img = nib.load(img_in)
        
#Padding function: reshape the image to (max_dim, max_dim, max_dim) with same resolution and an orientation of 'LAS'
def padding(img, new_name):
    d = max(img.header.get_data_shape())
    new_img = nib.processing.conform(img, out_shape=(d, d, d), \
    voxel_size = img.header.get_zooms(), order=0, cval=0, orientation='LAS', out_class=None)
    nib.save(new_img, new_name)

padding(img, img_out)

EOF
fi
}

#################
## Create a python script "brain-finalsurfs-edit.py"
#################
script_brain-finalsurfs-edit()
{
if [ ! -f "$SUBJECTS_DIR/brain-finalsurfs-edit.py" ]
then

cat > $SUBJECTS_DIR/brain-finalsurfs-edit.py <<EOF
import os
import numpy as np #To compute motion
import nibabel as nib #To edit MRI images
import sys #To add arguments
import copy #For deepcopy

#Outside parameters
path_bf = sys.argv[1] #Brain.finalsurfs without the cerebellum
path_gmbm = sys.argv[2] #Gray Matter binary mask at 128 (by default) (based on ribbon-edit.mgz)
path_out = sys.argv[3] #Absolute path of the output brain.finalsurfs

#Load brain.finalsurfs without the cerebellum
if not os.path.isfile(path_bf):
    raise FileNotFoundError("The following path doesn't exist: " + path_bf)
else:
    img_bf = nib.load(path_bf)
data_bf = img_bf.get_fdata()
(a,b,c)=img_bf.header.get_data_shape()

#Load Gray Matter binary mask at 128
if not os.path.isfile(path_gmbm):
    raise FileNotFoundError("The following path doesn't exist: " + path_gmbm)
else:
    img_gmbm = nib.load(path_gmbm)
data_gmbm = img_gmbm.get_fdata()
(e,f,g)=img_gmbm.header.get_data_shape()

#Test if both images have same size
if (a!=e or b!=f or c!=g):
    sys.exit(0)

#List of coordinates of the 26-nearest-neighbors around (0,0,0) plus itself
motion = np.transpose(np.indices((3,3,3)) - 1).reshape(-1, 3)

#Deepcopy of data_bf
data_bf_new = copy.deepcopy(data_bf)

#Change Gray Matter in data_bf_new to: 80.0/mean_neigbours*val
for x in range(a):
    for y in range(b):
        for z in range(c):
            if int(data_gmbm[x,y,z]) == 128:
                val=data_bf[x,y,z]
                n_coordinates = motion + [[x, y, z]]
                mean=0
                nn=0
                for (k,n,m) in n_coordinates:
                    if int(data_gmbm[k,n,m]) == 128:
                        mean+=data_bf[k,n,m]
                        nn+=1
                data_bf_new[x,y,z]=80.0/(mean/nn)*val

#Create and save new image
img_bf_new = nib.Nifti1Image(data_bf_new, img_bf.affine.copy())
nib.save(img_bf_new, path_out)


EOF
fi
}


#################
## Create an expert_file.txt for FreeSurfer
#################
script_expert_file()
{
if [ ! -f "$SUBJECTS_DIR/expert_file.txt" ]
then
cat > $SUBJECTS_DIR/expert_file.txt <<EOF
mris_inflate -n 30
EOF
fi
}

#################
## Function to both print and launch commands ($2) with a description ($1)
#################
cmd () {
if [ -z "$1" ] #First variable is empty
then
	Echo "
$2"
else
	Echo "
#---------------------------------
#@# $1: $current_date_time

$2"
fi
eval $2
}

O="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER"

#################
## Input files
#################
T1="$SUBJECTS_DIR/$SUBJID/mri/T1.mgz"

RB_ALL_WITHSKULL="$FREESURFER_HOME/average/RB_all_withskull_2016-05-10.vc700.gca"
TALAIRACH_WITH_SKULL="$SUBJECTS_DIR/$SUBJID/mri/transforms/talairach_with_skull.lta"
NU="$SUBJECTS_DIR/$SUBJID/mri/nu.mgz"
RB_ALL="$FREESURFER_HOME/average/RB_all_2020-01-02.gca"
CTRL_PTS="$SUBJECTS_DIR/$SUBJID/mri/ctrl_pts.mgz"
CC_UP="$SUBJECTS_DIR/$SUBJID/mri/transforms/cc_up.lta"

WM="$SUBJECTS_DIR/$SUBJID/mri/wm.mgz"

SUBCORTICALMASSLUT="$FREESURFER_HOME/SubCorticalMassLUT.txt"

# Curv + Thickness + Stats + APARC + APEG ...
declare -a CD_APARC_ATLAS=("$FREESURFER_HOME/average/lh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs" "$FREESURFER_HOME/average/rh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs")
declare -a DKT_APARC_ATLAS=("$FREESURFER_HOME/average/lh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs" "$FREESURFER_HOME/average/rh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs")
WMPARC_STATS_LUT="$FREESURFER_HOME/WMParcStatsLUT.txt"
ASEG_STATS_LUT="$FREESURFER_HOME/ASegStatsLUT.txt"
FSAVERAGE="$FREESURFER_HOME/subjects/fsaverage"
COLORTABLE_VPNL_TXT="$FREESURFER_HOME/average/colortable_vpnl.txt"
COLORTABLE_BA_TXT="$FREESURFER_HOME/average/colortable_BA.txt"

declare -a FOLDING_ATLAS_ACFB40=("$FREESURFER_HOME/average/lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif" "$FREESURFER_HOME/average/rh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif")
declare -a DKAPARC_ATLAS_ACFB40=("$FREESURFER_HOME/average/lh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs" "$FREESURFER_HOME/average/rh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs")

RAWAVG="$SUBJECTS_DIR/$SUBJID/mri/rawavg.mgz"
ORIG="$SUBJECTS_DIR/$SUBJID/mri/orig.mgz"


#ADDED for mri_segstats
TALAIRACH_XFM="$SUBJECTS_DIR/$SUBJID/mri/transforms/talairach.xfm"
TALAIRACH_XFM_COPY="$O/mri/transforms/talairach.xfm"

#################
## Output files
#################
IMAGE_PADDED="$SUBJECTS_DIR/image-padded.mgz"

RIBBON_PADDED="$O/mri/ribbon-precorrection.mgz"
SUBCORTICAL_PADDED="$O/mri/subcortical-precorrection.mgz"
RIBBON_EDIT="$O/mri/ribbon-edit.mgz"
SUBCORTICAL_EDIT="$O/mri/subcortical-edit.mgz"

SUBCORTICAL_MASK="$O/mri/subcortical-mask.mgz"
BRAIN_MASK="$O/mri/brain-mask.mgz"

T1_MASKED="$O/mri/T1-masked.mgz"

TALAIRACH="$O/mri/transforms/talairach.lta"
NORM="$O/mri/norm.mgz"
TALAIRACH_M3Z="$O/mri/transforms/talairach.m3z"
ASEG_AUTO_NOCCSEG="$O/mri/aseg.auto_noCCseg.mgz"
ASEG_AUTO="$O/mri/aseg.auto.mgz"
ASEG_PRESURF="$O/mri/aseg.presurf.mgz"

BRAIN="$O/mri/brain.mgz"
BRAIN_FINALSURFS="$O/mri/brain.finalsurfs.mgz"

WM_BMASK="$O/mri/wm-bmask.mgz"
WM_MASK="$O/mri/wm-mask.mgz"
WM_CONCAT="$O/mri/wm-concat.mgz"
WM_BMASK_250="$O/mri/wm-bmask-250.mgz"
WM_ASEGEDIT="$O/mri/wm-asegedit.mgz"
WM_EDITED="$O/mri/wm.mgz" # Use this name for mri_fix_topology

declare -a FILLED_PRETRESS=("$O/mri/filled_pretress_lh.mgz" "$O/mri/filled_pretress_rh.mgz")
declare -a ORIG_NOFIX_PREDEC=("$O/surf/lh.orig.nofix.predec" "$O/surf/rh.orig.nofix.predec")
declare -a ORIG_NOFIX=("$O/surf/lh.orig.nofix" "$O/surf/rh.orig.nofix")

declare -a SMOOTHW_NOFIX=("$O/surf/lh.smoothwm.nofix" "$O/surf/rh.smoothwm.nofix")
declare -a INFLATED_NOFIX=("$O/surf/lh.inflated.nofix" "$O/surf/rh.inflated.nofix")
declare -a QSPHERE_NOFIX=("$O/surf/lh.qsphere.nofix" "$O/surf/rh.qsphere.nofix")
declare -a ORIG_PREMESH=("$O/surf/lh.orig.premesh" "$O/surf/rh.orig.premesh")
declare -a ORIG=("$O/surf/lh.orig" "$O/surf/rh.orig")

#RJR ADDED
declare -a INFLATED=("$O/surf/lh.inflated" "$O/surf/rh.inflated")
declare -a SMOOTHW=("$O/surf/lh.smoothwm" "$O/surf/rh.smoothwm")
declare -a CURV=("$O/surf/lh.curv" "$O/surf/rh.curv")

declare -a CORTEX_LABEL=("$O/label/lh.cortex.label" "$O/label/rh.cortex.label")
declare -a CORTEX_HIPAMYG_LABEL=("$O/label/lh.cortex+hipamyg.label" "$O/label/rh.cortex+hipamyg.label")

declare -a GM_BMASK=("$O/mri/gm-bmask_lh.mgz" "$O/mri/gm-bmask_rh.mgz")
declare -a BMASK=("$O/mri/bmask_lh.mgz" "$O/mri/bmask_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB=("$O/mri/brain.finalsurfs_no_cereb_lh.mgz" "$O/mri/brain.finalsurfs_no_cereb_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80=("$O/mri/brain.finalsurfs_no_cereb_uniform_gm_80_lh.mgz" "$O/mri/brain.finalsurfs_no_cereb_uniform_gm_80_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB_EDITED=("$O/mri/brain.finalsurfs_no_cereb_edited_lh.mgz" "$O/mri/brain.finalsurfs_no_cereb_edited_rh.mgz")

declare -a AUTODET_NEW_GW_STATS=("$O/surf/autodet-new.gw.stats.lh.dat" "$O/surf/autodet-new.gw.stats.rh.dat")

declare -a RIBBON_EDIT_PIAL=("$O/surf/lh.ribbon_edit.pial" "$O/surf/rh.ribbon_edit.pial")
declare -a RIBBON_EDIT_PIAL_SECOND_PASS_i=("$O/surf/lh.ribbon_edit-second-pass_i.pial" "$O/surf/rh.ribbon_edit-second-pass_i.pial")
declare -a RIBBON_EDIT_PIAL_SECOND_PASS_r=("$O/surf/lh.ribbon_edit-second-pass_r.pial" "$O/surf/rh.ribbon_edit-second-pass_r.pial")
declare -a RIBBON_EDIT_PIAL_SECOND_PASS_i_r=("$O/surf/lh.ribbon_edit-second-pass_i_r.pial" "$O/surf/rh.ribbon_edit-second-pass_i_r.pial")
declare -a RIBBON_EDIT_PIAL_SECOND_PASS_i_w=("$O/surf/lh.ribbon_edit-second-pass_i_w.pial" "$O/surf/rh.ribbon_edit-second-pass_i_w.pial")
declare -a RIBBON_EDIT_PIAL_SECOND_PASS_i_w_r=("$O/surf/lh.ribbon_edit-second-pass_i_w_r.pial" "$O/surf/rh.ribbon_edit-second-pass_i_w_r.pial")

declare -a RIBBON_EDIT_PIAL_THIRD_PASS_SMOOTH=("$O/surf/lh.ribbon_edit.smooth-third-pass.pial" "$O/surf/rh.ribbon_edit.smooth-third-pass.pial")

# Curv + Thickness + Stats + APARC + APEG ...
declare -a CURV=("$O/surf/lh.curv" "$O/surf/rh.curv")
declare -a AREA=("$O/surf/lh.area" "$O/surf/rh.area")
declare -a PIAL_CURV=("$O/surf/lh.curv.pial" "$O/surf/rh.curv.pial")
declare -a PIAL_AREA=("$O/surf/lh.area.pial" "$O/surf/rh.area.pial")
declare -a THICKNESS=("$O/surf/lh.thickness" "$O/surf/rh.thickness")
declare -a CURV_STATS=("$O/stats/lh.curv.stats" "$O/stats/rh.curv.stats")

declare -a WHITE=("$O/surf/lh.white" "$O/surf/rh.white") #Copy of the rh.orig surface
declare -a PIAL=("$O/surf/lh.pial" "$O/surf/rh.pial") #Copy of the pial surface


declare -a SPHERE=("$O/surf/lh.sphere" "$O/surf/rh.sphere")
declare -a SPHERE_REG=("$O/surf/lh.sphere.reg" "$O/surf/rh.sphere.reg")

RAWAVG_MASKED="$O/mri/rawavg.mgz"
ORIG_MASKED="$O/mri/orig.mgz"

declare -a CD_APARC_ANNOT=("$O/label/lh.aparc.a2009s.annot" "$O/label/rh.aparc.a2009s.annot")
declare -a DKT_APARC_ANNOT=("$O/label/lh.aparc.DKTatlas.annot" "$O/label/rh.aparc.DKTatlas.annot")

ASEG_PRESURF_HYPOS="$O/mri/aseg.presurf.hypos.mgz" 
ASEG="$O/mri/aseg.mgz"
APARC_PLUS_ASEG="$O/mri/aparc+aseg.mgz"

declare -a APARC_ANNOT=("$O/label/lh.aparc.annot" "$O/label/rh.aparc.annot")
APARC_A2009S_ASEG="$O/mri/aparc.a2009s+aseg.mgz"
declare -a APARC_A2009S_ANNOT=("$O/label/lh.aparc.a2009s.annot" "$O/label/rh.aparc.a2009s.annot")
APARC_DKT_ATLAS_ASEG="$O/mri/aparc.DKTatlas+aseg.mgz"

WMPARC="$O/mri/wmparc.mgz"
WMPARC_STATS="$O/stats/wmparc.stats"

APARC_ANNOT_CTAB="$O/label/aparc.annot.ctab"
declare -a APARC_STATS=("$O/stats/lh.aparc.stats" "$O/stats/rh.aparc.stats")

declare -a APARC_PIAL_STATS=("$O/stats/lh.aparc.pial.stats" "$O/stats/rh.aparc.pial.stats")

declare -a APARC_A2009S_STATS=("$O/stats/lh.aparc.a2009s.stats" "$O/stats/rh.aparc.a2009s.stats")
APARC_A2009S_CTAB="$O/label/aparc.annot.a2009s.ctab"
declare -a DKT_APARC_STATS=("$O/stats/lh.aparc.DKTatlas.stats" "$O/stats/rh.aparc.DKTatlas.stats")

DKT_APARC_CTAB="$O/label/aparc.annot.DKTatlas.ctab"

ASEG_STATS="$O/stats/aseg.stats"

declare -a BA1_EXVIVO_LABEL=("$O/label/lh.BA1_exvivo.label" "$O/label/rh.BA1_exvivo.label")
declare -a BA2_EXVIVO_LABEL=("$O/label/lh.BA2_exvivo.label" "$O/label/rh.BA2_exvivo.label")
declare -a BA3A_EXVIVO_LABEL=("$O/label/lh.BA3a_exvivo.label" "$O/label/rh.BA3a_exvivo.label")
declare -a BA3B_EXVIVO_LABEL=("$O/label/lh.BA3b_exvivo.label" "$O/label/rh.BA3b_exvivo.label")
declare -a BA4A_EXVIVO_LABEL=("$O/label/lh.BA4a_exvivo.label" "$O/label/rh.BA4a_exvivo.label")
declare -a BA4P_EXVIVO_LABEL=("$O/label/lh.BA4p_exvivo.label" "$O/label/rh.BA4p_exvivo.label")
declare -a BA6_EXVIVO_LABEL=("$O/label/lh.BA6_exvivo.label" "$O/label/rh.BA6_exvivo.label")
declare -a BA44_EXVIVO_LABEL=("$O/label/lh.BA44_exvivo.label" "$O/label/rh.BA44_exvivo.label")
declare -a BA45_EXVIVO_LABEL=("$O/label/lh.BA45_exvivo.label" "$O/label/rh.BA45_exvivo.label")
declare -a V1_EXVIVO_LABEL=("$O/label/lh.V1_exvivo.label" "$O/label/rh.V1_exvivo.label")
declare -a V2_EXVIVO_LABEL=("$O/label/lh.V2_exvivo.label" "$O/label/rh.V2_exvivo.label")
declare -a MT_EXVIVO_LABEL=("$O/label/lh.MT_exvivo.label" "$O/label/rh.MT_exvivo.label")
declare -a ENTORHINAL_EXVIVO_LABEL=("$O/label/lh.enthorinal_exvivo.label" "$O/label/rh.enthorinal_exvivo.label")
declare -a PERIRHINAL_EXVIVO_LABEL=("$O/label/lh.perirhinal_exvivo.label" "$O/label/rh.perirhinal_exvivo.label")

declare -a FG1_MPM_VPNL_LABEL=("$O/label/lh.FG1.mpm.vpnl.label" "$O/label/rh.FG1.mpm.vpnl.label")
declare -a FG2_MPM_VPNL_LABEL=("$O/label/lh.FG2.mpm.vpnl.label" "$O/label/rh.FG2.mpm.vpnl.label")
declare -a FG3_MPM_VPNL_LABEL=("$O/label/lh.FG3.mpm.vpnl.label" "$O/label/rh.FG3.mpm.vpnl.label")
declare -a FG4_MPM_VPNL_LABEL=("$O/label/lh.FG4.mpm.vpnl.label" "$O/label/rh.FG4.mpm.vpnl.label")
declare -a H0C1_MPM_VPNL_LABEL=("$O/label/lh.h0c1.mpm.vpnl.label" "$O/label/rh.h0c1.mpm.vpnl.label")
declare -a H0C2_MPM_VPNL_LABEL=("$O/label/lh.h0c2.mpm.vpnl.label" "$O/label/rh.h0c2.mpm.vpnl.label")
declare -a H0C3V_MPM_VPNL_LABEL=("$O/label/lh.h0c3v.mpm.vpnl.label" "$O/label/rh.h0c3v.mpm.vpnl.label")
declare -a H0C4V_MPM_VPNL_LABEL=("$O/label/lh.h0c4v.mpm.vpnl.label" "$O/label/rh.h0c4v.mpm.vpnl.label")

declare -a BA1_EXVIVO_THRESH_LABEL=("$O/label/lh.BA1_exvivo.thresh.label" "$O/label/rh.BA1_exvivo.thresh.label")
declare -a BA2_EXVIVO_THRESH_LABEL=("$O/label/lh.BA2_exvivo.thresh.label" "$O/label/rh.BA2_exvivo.thresh.label")
declare -a BA3A_EXVIVO_THRESH_LABEL=("$O/label/lh.BA3a_exvivo.thresh.label" "$O/label/rh.BA3a_exvivo.thresh.label")
declare -a BA3B_EXVIVO_THRESH_LABEL=("$O/label/lh.BA3b_exvivo.thresh.label" "$O/label/rh.BA3b_exvivo.thresh.label")
declare -a BA4A_EXVIVO_THRESH_LABEL=("$O/label/lh.BA4a_exvivo.thresh.label" "$O/label/rh.BA4a_exvivo.thresh.label")
declare -a BA4P_EXVIVO_THRESH_LABEL=("$O/label/lh.BA4p_exvivo.thresh.label" "$O/label/rh.BA4p_exvivo.thresh.label")
declare -a BA6_EXVIVO_THRESH_LABEL=("$O/label/lh.BA6_exvivo.thresh.label" "$O/label/rh.BA6_exvivo.thresh.label")
declare -a BA44_EXVIVO_THRESH_LABEL=("$O/label/lh.BA44_exvivo.thresh.label" "$O/label/rh.BA44_exvivo.thresh.label")
declare -a BA45_EXVIVO_THRESH_LABEL=("$O/label/lh.BA45_exvivo.thresh.label" "$O/label/rh.BA45_exvivo.thresh.label")
declare -a V1_EXVIVO_THRESH_LABEL=("$O/label/lh.V1_exvivo.thresh.label" "$O/label/rh.V1_exvivo.thresh.label")
declare -a V2_EXVIVO_THRESH_LABEL=("$O/label/lh.V2_exvivo.thresh.label" "$O/label/rh.V2_exvivo.thresh.label")
declare -a MT_EXVIVO_THRESH_LABEL=("$O/label/lh.MT_exvivo.thresh.label" "$O/label/rh.MT_exvivo.thresh.label")
declare -a ENTORHINAL_EXVIVO_THRESH_LABEL=("$O/label/lh.enthorinal_exvivo.thresh.label" "$O/label/rh.enthorinal_exvivo.thresh.label")
declare -a PERIRHINAL_EXVIVO_THRESH_LABEL=("$O/label/lh.perirhinal_exvivo.thresh.label" "$O/label/rh.perirhinal_exvivo.thresh.label")

declare -a BA_EXVIVO_STATS=("$O/stats/lh.BA_exvivo.stats" "$O/stats/rh.BA_exvivo.stats")
declare -a BA_EXVIVO_ANNOT=("$O/label/lh.BA_exvivo.annot" "$O/label/rh.BA_exvivo.annot")
BA_EXVIVO_CTAB="$O/label/BA_exvivo.ctab"

declare -a BA_EXVIVO_THRESH_STATS=("$O/stats/lh.BA_exvivo.thresh.stats" "$O/stats/rh.BA_exvivo.thresh.stats")
declare -a BA_EXVIVO_THRESH_ANNOT=("$O/label/lh.BA_exvivo.thresh.annot" "$O/label/rh.BA_exvivo.thresh.annot")
BA_EXVIVO_THRESH_CTAB="$O/label/BA_exvivo.thresh.ctab"

#################
## New invocation in report.sh and create
#################
#Test if $OUTPUT_FOLDER folder already exist, and if it does not, create one
CreateScripts

Echo "
#*******************
#*******************
#*******************
# New invocation of ribbon_edit_script.sh

# Given subjid: $SUBJID"

#################
## FreeSurfer 7.4.1 on $IMAGE creating $SUBJECTS_DIR/$SUBJID folder
#################
if ((FS == 1 && TAG < 0))
then
# Test if user provided $IMAGE
: ${IMAGE:?Missing argument -i}
Echo "# Given image: $IMAGE"

if [ -d "$SUBJECTS_DIR/$SUBJID" ]
then
	Echo "Do not re-run FreeSurfer on same SUBJID: $SUBJECTS_DIR/$SUBJID"
	exit 1
fi

cmd "Compensate for future translation in FreeSurfer" \
"python $SUBJECTS_DIR/nifti_padding.py $IMAGE $IMAGE_PADDED padding"

#-xopts-overwrite is used when expert file already used before
cmd "Apply recon-all -autorecon 1 and 2 on $IMAGE_PADDED" \
"recon-all -autorecon1 -autorecon2 -s $SUBJID -i $IMAGE_PADDED -hires -parallel -openmp 4 -expert expert_file.txt -xopts-overwrite" 
fi

#################
## Convert ribbon and subcortical
#################
if ((TAG<=0))
then
# Test if user provided RIBBON and SUBCORTICAL
: ${RIBBON:?Missing argument -b} ${SUBCORTICAL:?Missing argument -c}
Echo "# Given ribbon: $RIBBON"
Echo "# Given subcortical: $SUBCORTICAL"

CreateFolders

cmd "Use script $SUBJECTS_DIR/nifti_padding.py on $RIBBON" \
"python $SUBJECTS_DIR/nifti_padding.py $RIBBON $RIBBON_PADDED padding"

cmd "Use script $SUBJECTS_DIR/nifti_padding.py on $SUBCORTICAL" \
"python $SUBJECTS_DIR/nifti_padding.py $SUBCORTICAL $SUBCORTICAL_PADDED padding"

#Necessary for correcting the orientation of the image
cmd "Convert $RIBBON_PADDED" \
"mri_convert $RIBBON_PADDED $RIBBON_EDIT -rt nearest -ns 1 --conform_min"

cmd "Convert $SUBCORTICAL_PADDED" \
"mri_convert $SUBCORTICAL_PADDED $SUBCORTICAL_EDIT -rt nearest -ns 1 --conform_min"
fi

#################
## Exctract labels from ribbon and subcortical into brain_mask
#################
if ((TAG<=1))
then
cmd "Extract labels from $SUBCORTICAL_EDIT (Cerebellum, Medulla oblongata, Pons and Midbrain) into $SUBCORTICAL_MASK" \
"mri_extract_label $SUBCORTICAL_EDIT $LABELS_SUBCORTICAL $SUBCORTICAL_MASK"

cmd "Concatenate $RIBBON_EDIT with $SUBCORTICAL_MASK into $BRAIN_MASK" \
"mri_concat --i $RIBBON_EDIT --i $SUBCORTICAL_MASK --o $BRAIN_MASK --combine"
fi

#################
## Recompute brain.finalsurfs.mgz
#################
if ((TAG<=2))
then
cmd "Mask $T1 with $BRAIN_MASK into $T1_MASKED" \
"mri_mask $T1 $BRAIN_MASK $T1_MASKED"
fi

# Recompute recon-all steps starting at EM Register up to brain.finalsurfs.mgz
if ((TAG<=3))
then
cmd "EM Register: mri_em_register" \
"mri_em_register -uns 3 -mask $T1_MASKED $NU $RB_ALL $TALAIRACH"

cmd "CA Normalize: mri_ca_normalize" \
"mri_ca_normalize -c $CTRL_PTS -mask $T1_MASKED $NU $RB_ALL $TALAIRACH $NORM"

cmd "CA Register: mri_ca_register" \
"mri_ca_register -nobigventricles -T $TALAIRACH -align-after -mask $T1_MASKED $NORM $RB_ALL $TALAIRACH_M3Z"

cmd "Subcortical Segment: mri_ca_label" \
"mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align $NORM $TALAIRACH_M3Z $RB_ALL $ASEG_AUTO_NOCCSEG"

cmd "CC Segment: mri_cc" \
"mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta $CC_UP $SUBJID/$OUTPUT_FOLDER" # Function "mri_cc" add "mri/" to $ASEG_AUTO_NOCCSEG and $ASEG_AUTO

cmd "Merge ASeg" \
"cp $ASEG_AUTO $ASEG_PRESURF"

cmd "Intensity Normalize" \
"mri_normalize -seed 1234 -mprage -noconform -aseg $ASEG_PRESURF -mask $T1_MASKED $NORM $BRAIN"

#Not needed
cmd "Mask BFS" \
"mri_mask -T 5 $BRAIN $T1_MASKED $BRAIN_FINALSURFS"
#OR
#cmd "Copy" \
#"cp $BRAIN $BRAIN_FINALSURFS"
fi

#################
## Compute wm.mgz : wm-bmask AND if(wm == 250 & wm-bmask), then wm-mask = 250
#################
if ((TAG<=4))
then
# Extract white matter from ribbon-edit to create wm-bmask.mgz
cmd "Extract WM from $RIBBON_EDIT" \
"mri_extract_label $RIBBON_EDIT ${LABEL_RIBBON_WM[0]} ${LABEL_RIBBON_WM[1]} $WM_BMASK" #0/128 binary mask
fi

# Compute WM_EDIT based on BRAIN_FINALSURFS masked by WM_BMASK
if ((TAG<=5))
then
cmd "Concatenate $WM_BMASK with $WM into $WM_CONCAT" \
"mri_concat --i $WM_BMASK --i $WM --o $WM_CONCAT --sum" #ROI at 378 (128+250)

cmd "Binarize $WM_CONCAT at 251 into $WM_BMASK_250" \
"mri_binarize --i $WM_CONCAT --o $WM_BMASK_250 --match 378"

cmd "Replace 1 by 250 into $WM_BMASK_250" \
"mri_binarize --i $WM_BMASK_250 --o $WM_BMASK_250 --replace 1 250"

# May also use $BRAIN_FINALSURFS
cmd "Mask $BRAIN with $WM_BMASK into $WM_MASK" \
"mri_mask -T 5 $BRAIN $WM_BMASK $WM_MASK"

cmd "Concatenate $WM_MASK with $WM_BMASK_250 into $WM_ASEGEDIT" \
"mri_concat --i $WM_MASK --i $WM_BMASK_250 --o $WM_ASEGEDIT --max"

cmd "Pretess $WM_ASEGEDIT: Solve connectivity issue" \
"mri_pretess $WM_ASEGEDIT wm $NORM $WM_EDITED"
fi

#################
## Compute ORIG: Don't need mri_fill, use ribbon-edit wm
#################
if ((TAG<=6))
then
for (( i=0; i<2; i++ ));
do
	if ((HEMI>=0 && HEMI!=i)); #Cases when the current hemi ($i) is not to be processed
	then
		continue;
	else
	# Compute directly ORIG_NOFIX
	cmd "${H[$i]} Pretress WM from $RIBBON_EDIT" \
	"mri_pretess $RIBBON_EDIT ${LABEL_RIBBON_WM[$i]} $NORM ${FILLED_PRETRESS[$i]}"

	cmd "${H[$i]} Tessellate WM surf" \
	"mri_tessellate ${FILLED_PRETRESS[$i]} ${LABEL_RIBBON_WM[$i]} ${ORIG_NOFIX_PREDEC[$i]}"

	cmd "${H[$i]} Extract main component WM surf" \
	"mris_extract_main_component ${ORIG_NOFIX_PREDEC[$i]} ${ORIG_NOFIX_PREDEC[$i]}"

	cmd "${H[$i]} Remesh WM surf" \
	"mris_remesh --desired-face-area 0.5 --input ${ORIG_NOFIX_PREDEC[$i]} --output ${ORIG_NOFIX[$i]}"

	# Smooth 1
	cmd "${H[$i]} Smooth WM surf" \
	"mris_smooth -n 1 -nw -seed 1234 ${ORIG_NOFIX[$i]} ${SMOOTHW_NOFIX[$i]}"

	# Inflate 1
	cmd "${H[$i]} Inflate WM surf" \
	"mris_inflate -no-save-sulc -n 30 ${SMOOTHW_NOFIX[$i]} ${INFLATED_NOFIX[$i]}"

	# Sphere 1
	cmd "${H[$i]} Make spherical WM surf" \
	"mris_sphere -q -p 6 -a 128 -seed 1234 ${INFLATED_NOFIX[$i]} ${QSPHERE_NOFIX[$i]}"

	# Fix topology
	cmd "${H[$i]} Fix tolpology WM surf" \
	"mris_fix_topology -mgz -sphere qsphere.nofix -inflated inflated.nofix -orig orig.nofix -out orig.premesh -ga -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]}"
	#-threads 1 #7.4.1: no difference seen; 7.4.0: small differences with different threads ? (https://surfer.nmr.mgh.harvard.edu/fswiki/ReleaseNotes)
	#Takes $SUBJID/$OUTPUT_FOLDER/mri/wm.mgz and brain.mgz
	
	# Remesh
	cmd "${H[$i]} Remesh WM surf" \
	"mris_remesh --remesh --iters 3 --input ${ORIG_PREMESH[$i]} --output ${ORIG[$i]}"

	# Remove intersection
	cmd "${H[$i]} Remove intersection" \
	"mris_remove_intersection ${ORIG[$i]} ${ORIG[$i]}"
	
	# Smooth 2 : Used for mris_curvature_stats later
	cmd "${H[$i]} Smooth WM surf" \
	"mris_smooth -n 1 -nw -seed 1234 ${ORIG[$i]} ${SMOOTHW[$i]}"
	
	# INFLATE 2 
	cmd "${H[$i]} Inflate to produce sulc file" \
	"mris_inflate ${ORIG[$i]} ${INFLATED[$i]}"

	# Compute WM CURV
	cmd "${H[$i]} Generate curv file from orig" \
	"mris_place_surface --curv-map ${ORIG[$i]} 2 10 ${CURV[$i]}"
	
	# Compute WM AREA
	cmd "${H[$i]} White area" \
 	"mris_place_surface --area-map ${ORIG[$i]} ${AREA[$i]}"
	fi
done
fi

#################
## Edit brain.finalsurfs with GM ribbon
#################
if ((TAG<=7))
then
for (( i=0; i<2; i++ ));
do
	if ((HEMI>=0 && HEMI!=i)); #Cases when the current hemi ($i) is not to be processed
	then
		continue;
	else
	# Extract brain from ribbon-edit to create bmask.mgz (no cerebellum), and use the latter on brain.finalsurfs
	cmd "${H[$i]} Extract GM from $RIBBON_EDIT" \
	"mri_extract_label $RIBBON_EDIT ${LABEL_RIBBON_GM[$i]} ${LABEL_RIBBON_WM[$i]} ${BMASK[$i]}" #0/128 binary mask

	cmd "${H[$i]} Replace 128 by 1 into ${BMASK[$i]}" \
	"mri_binarize --i ${BMASK[$i]} --o ${BMASK[$i]} --replace 128 1"

	cmd "${H[$i]} Mask $BRAIN_FINALSURFS with ${BMASK[$i]} into ${BRAIN_FINALSURFS_NO_CEREB[$i]}" \
	"mri_mask $BRAIN_FINALSURFS ${BMASK[$i]} ${BRAIN_FINALSURFS_NO_CEREB[$i]}"

	# Extract gray matter from ribbon-edit to create gm-bmask.mgz, and create with it bf_80
	cmd "${H[$i]} Extract GM from $RIBBON_EDIT" \
	"mri_extract_label $RIBBON_EDIT ${LABEL_RIBBON_GM[$i]} ${GM_BMASK[$i]}" #0/128 binary mask

	cmd "${H[$i]} Concatenate ${GM_BMASK[$i]} with ${BRAIN_FINALSURFS_NO_CEREB[$i]} into ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]}" \
	"mri_concat --i ${GM_BMASK[$i]} --i ${BRAIN_FINALSURFS_NO_CEREB[$i]} --o ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --max"

	cmd "${H[$i]} Replace 128 by 80 into ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_80[$i]}" \
	"mri_binarize --i ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --o ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --replace 128 80"

	# Use script brain-finalsurfs-edit.py to edit brain.finalsurfs.mgz
	cmd "${H[$i]} Use script $SUBJECTS_DIR/brain-finalsurfs-edit.py on ${BRAIN_FINALSURFS_NO_CEREB[$i]} with ${GM_BMASK[$i]}" \
	"python $SUBJECTS_DIR/brain-finalsurfs-edit.py ${BRAIN_FINALSURFS_NO_CEREB[$i]} ${GM_BMASK[$i]} ${BRAIN_FINALSURFS_NO_CEREB_EDITED[$i]}"
	fi
done
fi

#################
## Compute stats for surface
#################
if ((TAG<=8))
then
for (( i=0; i<2; i++ ));
do
	if ((HEMI>=0 && HEMI!=i)); #Cases when the current hemi ($i) is not to be processed
	then
		continue;
	else
	# Compute stats
	cmd "${H[$i]} Computes stats for pial surface" \
	"mris_autodet_gwstats --o ${AUTODET_NEW_GW_STATS[$i]} --i ${BRAIN_FINALSURFS_NO_CEREB_EDITED[$i]} --wm $WM_EDITED --surf ${ORIG[$i]}"
	
	if ((CHANGE_AUTODET==1))
	then
		# In order to improve pial surface, you can lower 'pial_border_low' to 20 
		# Change stats
		cmd "${H[$i]} Change stats" \
		"sed -i'' -e 's/^pial_border_low[^/n]*/pial_border_low $PIAL_BORDER_LOW/' ${AUTODET_NEW_GW_STATS[$i]}"
		#ex -s -c '%s/^pial_border_low.*/pial_border_low   $PIAL_BORDER_LOW/g|x' $AUTODET_NEW_GW_STATS_LH
	fi
	
	# Compute labels for pin-medial-wall
	cmd "${H[$i]} Label2label for cortex" \
	"mri_label2label --label-cortex ${ORIG[$i]} $ASEG_PRESURF 0 ${CORTEX_LABEL[$i]}"
	
	# Compute labels to remove HIPOCAMPUS AND AMYGDALA from pial surface in mris_place_surface
	cmd "${H[$i]} Label2label for cortex" \
	"mri_label2label --label-cortex ${ORIG[$i]} $ASEG_PRESURF 1 ${CORTEX_HIPAMYG_LABEL[$i]}"
	fi
done
fi

#################
## Compute pial surface: mris_place_surface
#################
if ((TAG<=9))
then
for (( i=0; i<2; i++ ));
do
	if ((HEMI>=0 && HEMI!=i)); #Cases when the current hemi ($i) is not to be processed
	then
		continue;
	else	
	cmd "${H[$i]} Computes pial surface" \
	"mris_place_surface --i ${ORIG[$i]} --o ${RIBBON_EDIT_PIAL[$i]} --nsmooth 0 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${ORIG[$i]} --invol ${BRAIN_FINALSURFS_NO_CEREB_EDITED[$i]} --threads 6 --white-surf ${ORIG[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg ${ASEG_PRESURF[$i]} --no-rip" #--rip-label $LH_CORTEX_HIPAMYG_LABEL"
	#Second pass
	cmd "${H[$i]} Computes pial surface - second pass" \
	"mris_place_surface --i ${RIBBON_EDIT_PIAL[$i]} --o ${RIBBON_EDIT_PIAL_SECOND_PASS_i[$i]} --nsmooth 0 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${ORIG[$i]} --invol ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --threads 6 --white-surf ${ORIG[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg $ASEG_PRESURF --no-rip"
		#Second pass
	cmd "${H[$i]} Computes pial surface - second pass" \
	"mris_place_surface --i ${ORIG[$i]} --o ${RIBBON_EDIT_PIAL_SECOND_PASS_r[$i]} --nsmooth 0 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${RIBBON_EDIT_PIAL[$i]} --invol ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --threads 6 --white-surf ${ORIG[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg $ASEG_PRESURF --no-rip"
		#Second pass
	cmd "${H[$i]} Computes pial surface - second pass" \
	"mris_place_surface --i ${RIBBON_EDIT_PIAL[$i]} --o ${RIBBON_EDIT_PIAL_SECOND_PASS_i_r[$i]} --nsmooth 0 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${RIBBON_EDIT_PIAL[$i]} --invol ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --threads 6 --white-surf ${ORIG[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg $ASEG_PRESURF --no-rip"
		#Second pass BEST RESULT (i_w)
	cmd "${H[$i]} Computes pial surface - second pass" \
	"mris_place_surface --i ${RIBBON_EDIT_PIAL[$i]} --o ${RIBBON_EDIT_PIAL_SECOND_PASS_i_w[$i]} --nsmooth 0 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${ORIG[$i]} --invol ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --threads 6 --white-surf ${RIBBON_EDIT_PIAL[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg $ASEG_PRESURF --no-rip"
		#Second pass
	cmd "${H[$i]} Computes pial surface - second pass" \
	"mris_place_surface --i ${RIBBON_EDIT_PIAL[$i]} --o ${RIBBON_EDIT_PIAL_SECOND_PASS_i_w_r[$i]} --nsmooth 0 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${RIBBON_EDIT_PIAL[$i]} --invol ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --threads 6 --white-surf ${RIBBON_EDIT_PIAL[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg $ASEG_PRESURF --no-rip"
	fi
done
fi

#################
## Add smoothing to pial surface (mris_place_surface)
# #################
cat << EOF
EOF

if ((TAG<=10))
then
for (( i=0; i<2; i++ ));
do
	if ((HEMI>=0 && HEMI!=i)); #Cases when the current hemi ($i) is not to be processed
	then
		continue;
	else		
	cmd "${H[$i]} Smooths pial surface" \
	"mris_place_surface --i ${RIBBON_EDIT_PIAL_SECOND_PASS_i_w[$i]} --o ${RIBBON_EDIT_PIAL_THIRD_PASS_SMOOTH[$i]} --nsmooth 1 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${RIBBON_EDIT_PIAL_SECOND_PASS_i_w[$i]} --invol $BRAIN_FINALSURFS_NO_CEREB_EDITED --threads 6 --white-surf ${ORIG[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg $ASEG_PRESURF --no-rip" #--rip-label --rip-bg
	
	cmd "${H[$i]} Copy white surface to ${WHITE[$i]}" \
	"cp ${ORIG[$i]} ${WHITE[$i]}" 
	
	cmd "${H[$i]} Copy pial surface to ${PIAL[$i]}" \
	"cp ${RIBBON_EDIT_PIAL_THIRD_PASS_SMOOTH[$i]} ${PIAL[$i]}" 
	fi
done
fi

#################
## Add last steps of autorecon3: stats, aseg, labels
# #################
if ((TAG<=11))
then
for (( i=0; i<2; i++ ));
do
	if ((HEMI>=0 && HEMI!=i)); #Cases when the current hemi ($i) is not to be processed
	then
		continue;
	else
 	cmd "${H[$i]} pial curv" \
 	"mris_place_surface --curv-map ${PIAL[$i]} 2 10 ${PIAL_CURV[$i]}"
 	cmd "${H[$i]} pial area" \
 	"mris_place_surface --area-map ${PIAL[$i]} ${PIAL_AREA[$i]}"
 	cmd "${H[$i]} thickness" \
 	"mris_place_surface --thickness ${ORIG[$i]} ${PIAL[$i]} 20 5 ${THICKNESS[$i]}"
 	# same command again for "area and vertex vol rh" ??? 	
 	
 	cmd "${H[$i]} Curvature Stats" \
 	"mris_curvature_stats -m --writeCurvatureFiles -G -o ${CURV_STATS[$i]} -F smoothwm $SUBJID/$OUTPUT_FOLDER ${H[$i]} curv sulc"
	
	#Copies only for next command line (mris_volmask searching for surf/rh.white and surf/rh.pial)
	cmd "${H[$i]} Copy ${ORIG[$i]} to ${WHITE[$i]}" \
	"cp ${ORIG[$i]} ${WHITE[$i]}"
 	cmd "${H[$i]} Cortical ribbon mask" \
 	"mris_volmask --aseg_name aseg.presurf --label_left_white ${LABEL_RIBBON_WM[$i]} --label_left_ribbon ${LABEL_RIBBON_GM[$i]} --label_right_white ${LABEL_RIBBON_WM[$i]} --label_right_ribbon ${LABEL_RIBBON_GM[$i]} --save_ribbon --${H[$i]}-only $SUBJID/$OUTPUT_FOLDER" #Searching for surf/rh.white and surf/rh.pial, and --surf_white and --surf_pial don't help
	#cmd "Copy $RH_SMOOTHW_NOFIX to $RH_SMOOTHW" \
	#"cp $RH_SMOOTHW_NOFIX $RH_SMOOTHW"
  	#cmd "Inflation2 rh" \
 	#"mris_inflate -n 30 $RH_SMOOTHW $RH_INFLATED"
  	cmd "${H[$i]} Sphere" \
 	"mris_sphere -seed 1234 ${INFLATED[$i]} ${SPHERE[$i]}"
 	cmd "${H[$i]} Surf Reg" \
 	"mris_register -curv ${SPHERE[$i]} ${FOLDING_ATLAS_ACFB40[$i]} ${SPHERE_REG[$i]}"
 	
 	
 	cmd "${H[$i]} Cortical Parc" \
 	"mris_ca_label -l ${CORTEX_LABEL[$i]} -aseg $ASEG_PRESURF -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]} ${SPHERE_REG[$i]} ${DKAPARC_ATLAS_ACFB40[$i]} ${APARC_ANNOT[$i]}"
 	cmd "${H[$i]} Cortical Parc 2" \
 	"mris_ca_label -l ${CORTEX_LABEL[$i]} -aseg $ASEG_PRESURF -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]} ${SPHERE_REG[$i]} ${CD_APARC_ATLAS[$i]} ${CD_APARC_ANNOT[$i]}"
 	#Use surf/rh.smoothwm and surf/rh.sphere.reg
 	cmd "${H[$i]} Cortical Parc 3" \
 	"mris_ca_label -l ${CORTEX_LABEL[$i]} -aseg $ASEG_PRESURF -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]} ${SPHERE_REG[$i]} ${DKT_APARC_ATLAS[$i]} ${DKT_APARC_ANNOT[$i]}"
 	cmd "${H[$i]} Copy $RAWAVG to $RAWAVG_MASKED" \
 	"cp $RAWAVG $RAWAVG_MASKED"
 	cmd "${H[$i]} Mask $RAWAVG_MASKED with $BRAIN_MASK into $RAWAVG_MASKED" \
"mri_mask $RAWAVG_MASKED $BRAIN_MASK $RAWAVG_MASKED"
 	cmd "${H[$i]} Copy $ORIG to $ORIG_MASKED" \
 	"cp $ORIG $ORIG_MASKED"
 	cmd "${H[$i]} Mask $ORIG_MASKED with $BRAIN_MASK into $ORIG_MASKED" \
"mri_mask $ORIG_MASKED $BRAIN_MASK $ORIG_MASKED"
 	cmd "${H[$i]} Change SUBJECTS_DIR to SUBJECTS_DIR/SUBJID" \
 	"export SUBJECTS_DIR=$SUBJECTS_DIR/$SUBJID"
 	cmd "${H[$i]} WM/GM Contrast" \
 	"pctsurfcon --s $OUTPUT_FOLDER --${H[$i]}-only"
 	#Need to change $SUBJECTS_DIR because doesn't accept $SUBJID/$OUTPUT_FOLDER 
 	#Need $RAWAVG_COPY, $ORIG_COPY and $RH_APARC_ANNOT
 	#PROBLEM in $RH_APARC_ANNOT: # elements (127231) in rh.aparc.annot does not match # vertices (138041)
 	cmd "${H[$i]} Change back SUBJECTS_DIR/SUBJID to SUBJECTS_DIR" \
 	"export SUBJECTS_DIR=$(dirname $SUBJECTS_DIR)"
 	
 	cmd "${H[$i]} Relabel Hypointensities" \
 	"mri_relabel_hypointensities -${H[$i]} $ASEG_PRESURF $O/surf $ASEG_PRESURF_HYPOS"
 	cmd "${H[$i]} APas-to-ASeg" \
 	"mri_surf2volseg --o $ASEG --i $ASEG_PRESURF_HYPOS --fix-presurf-with-ribbon $RIBBON_EDIT --threads 1 --${H[$i]}-cortex-mask ${CORTEX_LABEL[$i]} --${H[$i]}-white ${ORIG[$i]} --${H[$i]}-pial ${RIBBON_EDIT_PIAL[$i]} --${H[$i]}"
 	cmd "${H[$i]} AParc-to-ASeg aparc" \
 	"mri_surf2volseg --o $APARC_PLUS_ASEG --label-cortex --i $ASEG --threads 1 --${H[$i]}-annot ${APARC_ANNOT[$i]} 2000 --${H[$i]}-cortex-mask ${CORTEX_LABEL[$i]} --${H[$i]}-white ${ORIG[$i]} --${H[$i]}-pial ${RIBBON_EDIT_PIAL[$i]} --${H[$i]}"
 	cmd "${H[$i]} AParc-to-ASeg aparc.a2009s" \
 	"mri_surf2volseg --o $APARC_A2009S_ASEG --label-cortex --i $ASEG --threads 1 --${H[$i]}-annot ${APARC_A2009S_ANNOT[$i]} 12100 --${H[$i]}-cortex-mask ${CORTEX_LABEL[$i]} --${H[$i]}-white ${ORIG[$i]} --${H[$i]}-pial ${RIBBON_EDIT_PIAL[$i]} --${H[$i]}"
 	cmd "${H[$i]} AParc-to-ASeg aparc.DKTatlas" \
 	"mri_surf2volseg --o $APARC_DKT_ATLAS_ASEG --label-cortex --i $ASEG --threads 1 --${H[$i]}-annot ${DKT_APARC_ANNOT[$i]} 2000 --${H[$i]}-cortex-mask ${CORTEX_LABEL[$i]} --${H[$i]}-white ${ORIG[$i]} --${H[$i]}-pial ${RIBBON_EDIT_PIAL[$i]} --${H[$i]}"

 	cmd "${H[$i]} WMParc" \
 	"mri_surf2volseg --o $WMPARC --label-wm --i $APARC_PLUS_ASEG --threads 1 --${H[$i]}-annot ${APARC_ANNOT[$i]} 4000 --${H[$i]}-cortex-mask ${CORTEX_LABEL[$i]} --${H[$i]}-white ${ORIG[$i]} --${H[$i]}-pial ${RIBBON_EDIT_PIAL[$i]} --${H[$i]}"	
 	
	cmd "${H[$i]} Copy $TALAIRACH_XFM to $TALAIRACH_XFM_COPY" \
	"cp $TALAIRACH_XFM $TALAIRACH_XFM_COPY"
 	cmd "${H[$i]} WMParc stats" \
 	"mri_segstats --seed 1234 --seg $WMPARC --sum $WMPARC_STATS --pv $NORM --excludeid 0 --brainmask $BRAIN_FINALSURFS --in $NORM --in-intensity-name norm --in-intensity-units MR --subject $SUBJID/$OUTPUT_FOLDER --surf-wm-vol --ctab $WMPARC_STATS_LUT --etiv --no-global-stats"
 	#Needs mri/transforms/talairach.xfm
 	
 	cmd "${H[$i]} Change SUBJID" \
 	"export SUBJID=$SUBJID/$OUTPUT_FOLDER"
 	cmd "${H[$i]} Parcellation Stats White" \
 	"mris_anatomical_stats -th3 -mgz -noglobal -cortex ${CORTEX_LABEL[$i]} -f ${APARC_STATS[$i]} -b -a ${APARC_ANNOT[$i]} -c $APARC_ANNOT_CTAB $SUBJID ${H[$i]} white"
 	cmd "${H[$i]} Parcellation Stats Pial" \
 	"mris_anatomical_stats -th3 -mgz -noglobal -cortex ${CORTEX_LABEL[$i]} -f ${APARC_PIAL_STATS[$i]} -b -a ${APARC_ANNOT[$i]} -c $APARC_ANNOT_CTAB $SUBJID ${H[$i]} pial"
 	cmd "${H[$i]} Parcellation Stats 2 " \
 	"mris_anatomical_stats -th3 -mgz -noglobal -cortex ${CORTEX_LABEL[$i]} -f ${APARC_A2009S_STATS[$i]} -b -a ${APARC_A2009S_ANNOT[$i]} -c $APARC_A2009S_CTAB $SUBJID ${H[$i]} white"
 	cmd "${H[$i]} Parcellation Stats 3" \
 	"mris_anatomical_stats -th3 -mgz -noglobal -cortex ${CORTEX_LABEL[$i]} -f ${DKT_APARC_STATS[$i]} -b -a ${DKT_APARC_ANNOT[$i]} -c $DKT_APARC_CTAB $SUBJID ${H[$i]} white"
 	cmd "${H[$i]} Change back SUBJID" \
 	"export SUBJID=`echo "$SUBJID" | cut -d/ -f1-1`"
 	
 	cmd "ASeg Stats" \
 	"mri_segstats --seed 1234 --seg $ASEG --sum $ASEG_STATS --pv $NORM --empty --brainmask $BRAIN_FINALSURFS --brain-vol-from-seg --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in $NORM --in-intensity-name norm --in-intensity-units MR --etiv --euler --ctab $ASEG_STATS_LUT --subject $SUBJID/$OUTPUT_FOLDER --no-global-stats"
 	
 	cmd "Create Symlink of $FSAVERAGE folder in SUBJECTS_DIR" \
 	"ln -s $FSAVERAGE $SUBJECTS_DIR"
 	
 	cmd "${H[$i]} BA_exvivo Labels" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA1_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA1_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label BA2_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA2_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA2_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	
 	cmd "${H[$i]} mri_label2label BA3a_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA3a_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA3A_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label BA3b_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA3b_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA3B_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label BA4a_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA4a_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA4A_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label BA4p_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA4p_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA4P_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label BA6_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA6_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA6_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label BA44_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA44_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA44_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]}  mri_label2label BA45_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA45_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA45_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label V1_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.V1_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${V1_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label V2_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.V2_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${V2_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label MT_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.MT_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${MT_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label entorhinal_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.entorhinal_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${ENTORHINAL_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label perirhinal_exvivo" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.perirhinal_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel vPERIRHINAL_EXVIVO_LABEL --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label fg1_mpm_vpnl" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.FG1.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${FG1_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label fg2_mpm_vpnl" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.FG2.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${FG2_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label fg3_mpm_vpnl" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.FG3.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${FG3_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label fg4_mpm_vpnl" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.FG4.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${FG4_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label h0c1_mpm_vpnl" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.hOc1.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${H0C1_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label h0c2_mpm_vpnl" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.hOc2.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${H0C2_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label h0c3v_mpm_vpnl" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.hOc3v.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${H0C3V_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label h0c4v_mpm_vpnl" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.hOc4v.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${H0C4V_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	
 	cmd "${H[$i]} mri_label2label ctab" \
 	"mris_label2annot --s $SUBJID/$OUTPUT_FOLDER --ctab $COLORTABLE_VPNL_TXT --hemi ${H[$i]} --a mpm.vpnl --maxstatwinner --noverbose --l ${FG1_MPM_VPNL_LABEL[$i]} --l ${FG2_MPM_VPNL_LABEL[$i]} --l ${FG3_MPM_VPNL_LABEL[$i]} --l ${FG4_MPM_VPNL_LABEL[$i]} --l ${H0C1_MPM_VPNL_LABEL[$i]} --l ${H0C2_MPM_VPNL_LABEL[$i]} --l ${H0C3V_MPM_VPNL_LABEL[$i]} --l ${H0C4V_MPM_VPNL_LABEL[$i]}"
 	
 	cmd "${H[$i]} mri_label2label ba1_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA1_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA1_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label ba2_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA2_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA2_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	
 	cmd "${H[$i]} mri_label2label ba3a_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA3a_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA3A_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label ba3b_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA3b_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA3B_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label ba4a_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA4a_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA4A_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label ba4p_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA4p_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA4P_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label ba6_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA6_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA6_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label ba44_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA44_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA44_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label ba45_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA45_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA45_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label v1_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.V1_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${V1_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface "
 	cmd "${H[$i]} mri_label2label v2_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.V2_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${V2_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label mt_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.MT_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${MT_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label entorhinal_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.entorhinal_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${ENTORHINAL_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	cmd "${H[$i]} mri_label2label perirhinal_exvivo_thresh" \
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.perirhinal_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${PERIRHINAL_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
 	
 	cmd "${H[$i]} mri_label2label ctab" \
 	"mris_label2annot --s $SUBJID/$OUTPUT_FOLDER --hemi ${H[$i]} --ctab $COLORTABLE_BA_TXT --l ${BA1_EXVIVO_LABEL[$i]} --l ${BA2_EXVIVO_LABEL[$i]} --l ${BA3A_EXVIVO_LABEL[$i]} --l ${BA3B_EXVIVO_LABEL[$i]} --l ${BA4A_EXVIVO_LABEL[$i]} --l ${BA4P_EXVIVO_LABEL[$i]} --l ${BA6_EXVIVO_LABEL[$i]} --l ${BA44_EXVIVO_LABEL[$i]} --l ${BA45_EXVIVO_LABEL[$i]} --l ${V1_EXVIVO_LABEL[$i]} --l ${V2_EXVIVO_LABEL[$i]} --l ${MT_EXVIVO_LABEL[$i]} --l ${PERIRHINAL_EXVIVO_LABEL[$i]} --l ${ENTORHINAL_EXVIVO_LABEL[$i]} --a BA_exvivo --maxstatwinner --noverbose"
 	
 	cmd "${H[$i]} mri_label2label ctab thresh" \
 	"mris_label2annot --s $SUBJID/$OUTPUT_FOLDER --hemi ${H[$i]} --ctab $COLORTABLE_BA_TXT --l ${BA1_EXVIVO_THRESH_LABEL[$i]} --l ${BA2_EXVIVO_THRESH_LABEL[$i]} --l ${BA3A_EXVIVO_THRESH_LABEL[$i]} --l ${BA3B_EXVIVO_THRESH_LABEL[$i]} --l ${BA4A_EXVIVO_THRESH_LABEL[$i]} --l ${BA4P_EXVIVO_THRESH_LABEL[$i]} --l ${BA6_EXVIVO_THRESH_LABEL[$i]} --l ${BA44_EXVIVO_THRESH_LABEL[$i]} --l ${BA45_EXVIVO_THRESH_LABEL[$i]} --l ${V1_EXVIVO_THRESH_LABEL[$i]} --l ${V2_EXVIVO_THRESH_LABEL[$i]} --l ${MT_EXVIVO_THRESH_LABEL[$i]} --l ${PERIRHINAL_EXVIVO_THRESH_LABEL[$i]} --l ${ENTORHINAL_EXVIVO_THRESH_LABEL[$i]} --a BA_exvivo.thresh --maxstatwinner --noverbose"

 	cmd "${H[$i]} mris_anatomical_stats ctab" \
 	"mris_anatomical_stats -th3 -mgz -noglobal -f ${BA_EXVIVO_STATS[$i]} -b -a ${BA_EXVIVO_ANNOT[$i]} -c $BA_EXVIVO_CTAB $SUBJID/$OUTPUT_FOLDER ${H[$i]} white"
 	cmd "${H[$i]}mris_anatomical_stats ctab thresh" \
 	"mris_anatomical_stats -th3 -mgz -f ${BA_EXVIVO_THRESH_STATS[$i]} -noglobal -b -a ${BA_EXVIVO_THRESH_ANNOT[$i]} -c $BA_EXVIVO_THRESH_CTAB $SUBJID/$OUTPUT_FOLDER ${H[$i]} white"
	fi
done
fi
