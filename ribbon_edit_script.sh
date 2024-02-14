#!/bin/bash

#################
## Help
#################
Help ()
{
builtin echo "
AUTHOR: Benoît Verreman

LAST UPDATE: 2024-02-14

DESCRIPTION: 
Use ribbon and subcortical NIFTI files to recompute pial surface,
based on previously created <subjid> folder using Freesurfer 7.4.1
Create a log 'report.sh'.
Create a folder 'outputs' with all the output files.

PREREQUISITE:
Export SUBJECTS_DIR and FREESURFER_HOME correctly

Test if you have access to python: 'which python'
Install two python libraries: 'pip install nibabel scipy'

If you want to launch Freesurfer 7.4.1 recon-all pipeline using the script:
	Add argument -i
Else:
	Launch Freesurfer 7.4.1 command before using script: 
$ recon-all -s <subjid> -i <subject_image> -autorecon1 -autorecon2 -hires -parallel -openmp 4 -expert expert_file.txt
Prepare ribbon and subcortical NIFTI files (for step 0)
Put the script inside <subjid> folder

EXAMPLES:
$ bash ribbon_edit_script.sh -i 133019_T1w_acpc_dc_restore.nii.gz -s 133019 -r 133019_ribbon.nii.gz -c 133019_subcortical.nii.gz

$ bash ribbon_edit_script.sh -s <subjid> -t 8 -r #start from pial computation step (8), right hemisphere only (r)

PARAMETERS:
-i: Relative or absolute path to T1w image file

-s: Relative or absolute path to subjid folder (Necessary)
-r: Relative or absolute path to ribbon file
-c: Relative or absolute path to subcortical file

-h: Print this string, and exit

TAG
-t 0: (ribbons) Start with resizing RIBBON_EDIT and SUBCORTICAL
-t 1: (bmask) Start with BRAIN_MASK
-t 2: (maskT1) Start with T1_MASKED
-t 3: (brain.finalsurfs) Start with skull-stripping up to BRAIN_FINALSURFS
-t 4: (wm-bmask) Start the creation of WM_BMASK based on RIBBON_EDIT
-t 5: (wm) Start from computing WM based on WM_BMASK
-t 6: (orig) Start from computing orig surface based on wm from RIBBON_EDIT
-t 7: (stats) Start from computing stats
-t 8: (pial) Start from computing pial surface
-t 9: (smooth) Start from smoothing pial surface

HEMI
-r: Compute only right hemisphere surface
-l: Compute only left hemisphere surface

-d: Reset outputs folder and report.sh script

"
}

#################
## Default global variables
#################
current_date_time=$(date)
TAG=-1 # Start from beginning
HEMI=0 # Both hemispheres
FS=0 # Default: No recon-all
OUTPUT_FOLDER="outputs"
LABELS_SUBCORTICAL="5 15 29 30 32 31"
LABEL_RIBBON_WM_LH=2
LABEL_RIBBON_WM_RH=41
PIAL_BORDER_LOW=5

#################
## Manage flags
#################
unset -v IMAGE
unset -v SUBJID
unset -v RIBBON
unset -v SUBCORTICAL

#If a character is followed by :, then it needs an argument
VALID_ARGS="i:s:r:c:t:hlrd"

echo "$VALID_ARGS"

while getopts ${VALID_ARGS} opt; do
  case ${opt} in
    i)
        IMAGE=${OPTARG}
        FS=1
        ;;
    s)
        SUBJID=${OPTARG}
        ;;
    r)
        RIBBON=${OPTARG}
        ;;
    c)
        SUBCORTICAL=${OPTARG}
        ;;     
    t)
	TAG=${OPTARG}
	;;
    h)
	Help
	exit 1
	;;
    l)
	HEMI=1 #left hemisphere only
	;;
    r)
	HEMI=-1 #right hemisphere only
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
: ${SUBJID:?Missing argument --subjid or -s}

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
	mkdir $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/label;"
fi
}

Delete()
{
cmd "Reset $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER" \
"rm -r $SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER;"

rm $SUBJECTS_DIR/report.sh
}

#################
## Create a python script "mri_convert_correction_by_translation.py"
#################
script_nifti_padding()
{
if [ ! -f "$SUBJECTS_DIR/nifti_padding.py" ]
then
cat > $SUBJECTS_DIR/nifti_padding.py <<EOF
import os
import nibabel as nib
import nibabel.processing #Used in conform
import scipy.ndimage #Used in conform
import sys #To add parameters

# SUBJID directory
img_in = sys.argv[1]
img_out = sys.argv[2]

#Load image to be treated
if not os.path.isfile(img_in):
    raise FileNotFoundError("The following path doesn't exist: " + img_in)
else:
    img = nib.load(img_in)

#New MRI image
def padding(img, new_name):
    new_img = nib.processing.conform(img, out_shape=(311, 311, 311), \
    voxel_size=(0.7, 0.7, 0.7), order=0, cval=0, orientation='LAS', out_class=None)
    nib.save(new_img, new_name)
    
padding(img, img_out)
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
## Use python script "nifti_padding.py"
#################
nifti_padding()
{
python $SUBJECTS_DIR/nifti_padding.py $1 $2
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
ASEG_PRESURF="$SUBJECTS_DIR/$SUBJID/mri/aseg.presurf.mgz"

WM="$SUBJECTS_DIR/$SUBJID/mri/wm.mgz"

SUBCORTICALMASSLUT="$FREESURFER_HOME/SubCorticalMassLUT.txt"

#################
## Output files
#################
IMAGE_PADDED="$SUBJECTS_DIR/image-padded.mgz"

RIBBON_PADDED="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/ribbon-precorrection.mgz"
SUBCORTICAL_PADDED="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/subcortical-precorrection.mgz"
RIBBON_EDIT="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/ribbon-edit.mgz"
SUBCORTICAL_EDIT="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/subcortical-edit.mgz"

SUBCORTICAL_MASK="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/subcortical-mask.mgz"
BRAIN_MASK="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/brain-mask.mgz"

T1_CORRECTED="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/T1-corrected.mgz"
T1_MASKED="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/T1-masked.mgz"

TALAIRACH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/transforms/talairach.lta"
NORM="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/norm.mgz"
TALAIRACH_M3Z="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/transforms/talairach.m3z"
ASEG_AUTO_NOCCSEG="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/aseg.auto_noCCseg.mgz"
ASEG_AUTO="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/aseg.auto.mgz"
ASEG_PRESURF="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/aseg.presurf.mgz"

BRAIN="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/brain.mgz"
BRAIN_FINALSURFS="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/brain.finalsurfs.mgz"

WM_BMASK="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/wm-bmask.mgz"
WM_MASK="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/wm-mask.mgz"
WM_CONCAT="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/wm-concat.mgz"
WM_BMASK_250="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/wm-bmask-250.mgz"
WM_ASEGEDIT="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/wm-asegedit.mgz"
WM_EDITED="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/wm.mgz" # Use this name for mri_fix_topology

FILLED_PRETRESS_LH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/filled_pretress_lh.mgz"
FILLED_PRETRESS_RH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/mri/filled_pretress_rh.mgz"
LH_ORIG_NOFIX_PREDEC="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.orig.nofix.predec"
RH_ORIG_NOFIX_PREDEC="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.orig.nofix.predec"
LH_ORIG_NOFIX="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.orig.nofix"
RH_ORIG_NOFIX="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.orig.nofix"

LH_SMOOTHW_NOFIX="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.smoothwm.nofix"
RH_SMOOTHW_NOFIX="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.smoothwm.nofix"
LH_INFLATED_NOFIX="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.inflated.nofix"
RH_INFLATED_NOFIX="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.inflated.nofix"
LH_QSPHERE_NOFIX="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.qsphere.nofix"
RH_QSPHERE_NOFIX="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.qsphere.nofix"
LH_ORIG_PREMESH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.orig.premesh"
RH_ORIG_PREMESH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.orig.premesh"
LH_ORIG="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.orig"
RH_ORIG="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.orig"

LH_CORTEX_LABEL="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/label/lh.cortex.label"
RH_CORTEX_LABEL="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/label/rh.cortex.label"
LH_CORTEX_HIPAMYG_LABEL="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/label/lh.cortex+hipamyg.label"
RH_CORTEX_HIPAMYG_LABEL="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/label/rh.cortex+hipamyg.label"

AUTODET_NEW_GW_STATS_LH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/autodet-new.gw.stats.lh.dat"
AUTODET_NEW_GW_STATS_RH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/autodet-new.gw.stats.rh.dat"
LH_RIBBON_EDIT_PIAL="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.ribbon_edit.pial"
RH_RIBBON_EDIT_PIAL="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.ribbon_edit.pial"

LH_RIBBON_EDIT_PIAL_SMOOTH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.ribbon_edit.smooth.pial"
RH_RIBBON_EDIT_PIAL_SMOOTH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.ribbon_edit.smooth.pial"

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
: ${IMAGE:?Missing argument --image or -i}
Echo "# Given image: $IMAGE"

if [ -d "$SUBJECTS_DIR/$SUBJID" ]
then
	Echo "Do not re-run FreeSurfer on same SUBJID: $SUBJECTS_DIR/$SUBJID"
	exit 1
fi

cmd "Compensate for future translation in FreeSurfer" \
"nifti_padding $IMAGE $IMAGE_PADDED"

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
: ${RIBBON:?Missing argument --ribbon or -r} ${SUBCORTICAL:?Missing argument --subcortical or -c}
Echo "# Given ribbon: $RIBBON"
Echo "# Given subcortical: $SUBCORTICAL"

CreateFolders

cmd "Use script nifti_padding.py on $RIBBON" \
"nifti_padding $RIBBON $RIBBON_PADDED"

cmd "Use script nifti_padding.py on $SUBCORTICAL" \
"nifti_padding $SUBCORTICAL $SUBCORTICAL_PADDED"

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
# Extract white matter from ribbon-edit to create wm-bmask.mgz
if ((TAG<=4))
then
cmd "Extract WM from $RIBBON_EDIT" \
"mri_extract_label $RIBBON_EDIT $LABEL_RIBBON_WM_LH $LABEL_RIBBON_WM_RH $WM_BMASK" #0/128 binary mask
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
if ((HEMI>=0))
then
	# Compute directly ORIG_NOFIX
	cmd "Pretress lh WM from $RIBBON_EDIT" \
	"mri_pretess $RIBBON_EDIT $LABEL_RIBBON_WM_LH $NORM $FILLED_PRETRESS_LH"

	cmd "Tessellate lh WM surf" \
	"mri_tessellate $FILLED_PRETRESS_LH $LABEL_RIBBON_WM_LH $LH_ORIG_NOFIX_PREDEC"

	cmd "Extract main component lh WM surf" \
	"mris_extract_main_component $LH_ORIG_NOFIX_PREDEC $LH_ORIG_NOFIX_PREDEC"

	cmd "Remesh lh WM surf" \
	"mris_remesh --desired-face-area 0.5 --input $LH_ORIG_NOFIX_PREDEC --output $LH_ORIG_NOFIX"

	# Smooth 1
	cmd "Smooth lh WM surf" \
	"mris_smooth -n 3 -nw -seed 1234 $LH_ORIG_NOFIX $LH_SMOOTHW_NOFIX"

	# Inflate 1
	cmd "Inflate lh WM surf" \
	"mris_inflate -no-save-sulc -n 30 $LH_SMOOTHW_NOFIX $LH_INFLATED_NOFIX"

	# Sphere 1
	cmd "Make spherical lh WM surf" \
	"mris_sphere -q -p 6 -a 128 -seed 1234 $LH_INFLATED_NOFIX $LH_QSPHERE_NOFIX"

	# Fix topology
	cmd "Fix tolpology lh WM surf" \
	"mris_fix_topology -mgz -sphere qsphere.nofix -inflated inflated.nofix -orig orig.nofix -out orig.premesh -ga -seed 1234 $SUBJID/$OUTPUT_FOLDER lh"

	# Remesh
	cmd "Remesh lh WM surf" \
	"mris_remesh --remesh --iters 3 --input $LH_ORIG_PREMESH --output $LH_ORIG"

	# Remove intersection
	cmd "Remove intersection lh" \
	"mris_remove_intersection $LH_ORIG $LH_ORIG"
fi

if ((HEMI<=0))
then
	# Compute directly ORIG_NOFIX
	cmd "Pretress rh WM from $RIBBON_EDIT" \
	"mri_pretess $RIBBON_EDIT $LABEL_RIBBON_WM_RH $NORM $FILLED_PRETRESS_RH"

	cmd "Tessellate rh WM surf" \
	"mri_tessellate $FILLED_PRETRESS_RH $LABEL_RIBBON_WM_RH $RH_ORIG_NOFIX_PREDEC"

	cmd "Extract main component rh WM surf" \
	"mris_extract_main_component $RH_ORIG_NOFIX_PREDEC $RH_ORIG_NOFIX_PREDEC"

	cmd "Remesh rh WM surf" \
	"mris_remesh --desired-face-area 0.5 --input $RH_ORIG_NOFIX_PREDEC --output $RH_ORIG_NOFIX"

	# Smooth 1
	cmd "Smooth rh WM surf" \
	"mris_smooth -n 3 -nw -seed 1234 $RH_ORIG_NOFIX $RH_SMOOTHW_NOFIX"

	# Inflate 1
	cmd "Inflate rh WM surf" \
	"mris_inflate -no-save-sulc -n 30 $RH_SMOOTHW_NOFIX $RH_INFLATED_NOFIX"

	# Sphere 1
	cmd "Make spherical rh WM surf" \
	"mris_sphere -q -p 6 -a 128 -seed 1234 $RH_INFLATED_NOFIX $RH_QSPHERE_NOFIX"

	# Fix topology
	cmd "Fix tolpology rh WM surf" \
	"mris_fix_topology -mgz -sphere qsphere.nofix -inflated inflated.nofix -orig orig.nofix -out orig.premesh -ga -seed 1234 $SUBJID/$OUTPUT_FOLDER rh"

	# Remesh
	cmd "Remesh rh WM surf" \
	"mris_remesh --remesh --iters 3 --input $RH_ORIG_PREMESH --output $RH_ORIG"

	# Remove intersection
	cmd "Remove intersection rh" \
	"mris_remove_intersection $RH_ORIG $RH_ORIG"
fi
fi

#################
## Compute stats for surface
#################
if ((TAG<=7))
then
if ((HEMI>=0))
then
	# Compute stats
	cmd "Computes stats for lh pial surface" \
	"mris_autodet_gwstats --o $AUTODET_NEW_GW_STATS_LH --i $BRAIN_FINALSURFS --wm $WM_EDITED --surf $LH_ORIG"
	
	# In order to improve pial surface, you can lower 'pial_border_low' to 20 
	# Change stats
	cmd "Change stats" \
	"sed -i'' -e 's/^pial_border_low[^/n]*/pial_border_low $PIAL_BORDER_LOW/' $AUTODET_NEW_GW_STATS_LH"
	#ex -s -c '%s/^pial_border_low.*/pial_border_low   $PIAL_BORDER_LOW/g|x' $AUTODET_NEW_GW_STATS_LH
	
	# Compute labels for pin-medial-wall
	cmd "Label2label for lh cortex" \
	"mri_label2label --label-cortex $LH_ORIG $ASEG_PRESURF 0 $LH_CORTEX_LABEL"
	
	# Compute labels to remove HIPOCAMPUS AND AMYGDALA from pial surface in mris_place_surface
	cmd "Label2label for lh cortex" \
	"mri_label2label --label-cortex $LH_ORIG $ASEG_PRESURF 1 $LH_CORTEX_HIPAMYG_LABEL"

fi

if ((HEMI<=0))
then
	# Compute stats
	cmd "Computes stats for rh pial surface" \
	"mris_autodet_gwstats --o $AUTODET_NEW_GW_STATS_RH --i $BRAIN_FINALSURFS --wm $WM_EDITED --surf $RH_ORIG"
	
	# In order to improve pial surface, you can lower 'pial_border_low' to 20 	
	# Change stats
	cmd "Change stats" \
	"sed -i'' -e 's/^pial_border_low[^/n]*/pial_border_low $PIAL_BORDER_LOW/' $AUTODET_NEW_GW_STATS_RH"	
	#ex -s -c '%s/^pial_border_low.*/pial_border_low   $PIAL_BORDER_LOW/g|x' $AUTODET_NEW_GW_STATS_RH
	
	# Compute labels for pin-medial-wall
	cmd "Label2label for rh cortex" \
	"mri_label2label --label-cortex $RH_ORIG $ASEG_PRESURF 0 $RH_CORTEX_LABEL"

	# Compute labels to remove HIPOCAMPUS AND AMYGDALA from pial surface in mris_place_surface
	cmd "Label2label for rh cortex" \
	"mri_label2label --label-cortex $RH_ORIG $ASEG_PRESURF 1 $RH_CORTEX_HIPAMYG_LABEL"
fi
fi

#################
## Compute pial surface: mris_place_surface
#################
if ((TAG<=8))
then
if ((HEMI>=0))
then	
	cmd "Computes lh pial surface" \
	"mris_place_surface --i $LH_ORIG --o $LH_RIBBON_EDIT_PIAL --nsmooth 0 --adgws-in $AUTODET_NEW_GW_STATS_LH --pial --lh --repulse-surf $LH_ORIG --invol $BRAIN_FINALSURFS --threads 6 --white-surf $LH_ORIG --pin-medial-wall $LH_CORTEX_LABEL --seg $ASEG_PRESURF --no-rip" #--rip-label $LH_CORTEX_HIPAMYG_LABEL"
	#"mris_place_surface --i $LH_ORIG --o $LH_RIBBON_EDIT_PIAL --nsmooth 0 --adgws-in $AUTODET_NEW_GW_STATS_LH --pial --lh --repulse-surf $LH_ORIG --invol $BRAIN_FINALSURFS --threads 6 --white-surf $LH_ORIG --pin-medial-wall $LH_CORTEX_LABEL --seg $ASEG_PRESURF --no-rip"
fi

if ((HEMI<=0))
then
	cmd "Computes rh pial surface" \
	"mris_place_surface --i $RH_ORIG --o $RH_RIBBON_EDIT_PIAL --nsmooth 0 --adgws-in $AUTODET_NEW_GW_STATS_RH --pial --rh --repulse-surf $RH_ORIG --invol $BRAIN_FINALSURFS --threads 6 --white-surf $RH_ORIG --pin-medial-wall $RH_CORTEX_LABEL --seg $ASEG_PRESURF --no-rip" #--rip-label $RH_CORTEX_HIPAMYG_LABEL"
	#"mris_place_surface --i $RH_ORIG --o $RH_RIBBON_EDIT_PIAL --nsmooth 0 --adgws-in $AUTODET_NEW_GW_STATS_RH --pial --rh --repulse-surf $RH_ORIG --invol $BRAIN_FINALSURFS --threads 6 --white-surf $RH_ORIG --pin-medial-wall $RH_CORTEX_LABEL --seg $ASEG_PRESURF --no-rip"
fi
fi

#################
## Add smoothing to pial surface (mris_place_surface)
# #################
if ((TAG<=9))
then
if ((HEMI>=0))
then	
	cmd "Smooths lh pial surface" \
	"mris_place_surface --i $LH_RIBBON_EDIT_PIAL --o $LH_RIBBON_EDIT_PIAL_SMOOTH --nsmooth 1 --adgws-in $AUTODET_NEW_GW_STATS_LH --pial --lh --repulse-surf $LH_RIBBON_EDIT_PIAL --invol $BRAIN_FINALSURFS --threads 6 --white-surf $LH_RIBBON_EDIT_PIAL --pin-medial-wall $LH_CORTEX_LABEL --seg $ASEG_PRESURF --no-rip"
fi

if ((HEMI<=0))
then
	cmd "Smooths rh pial surface" \
 	"mris_place_surface --i $RH_RIBBON_EDIT_PIAL --o $RH_RIBBON_EDIT_PIAL_SMOOTH --nsmooth 1 --adgws-in $AUTODET_NEW_GW_STATS_RH --pial --rh --repulse-surf $RH_RIBBON_EDIT_PIAL --invol $BRAIN_FINALSURFS --threads 6 --white-surf $RH_ORIG --pin-medial-wall $RH_CORTEX_LABEL --seg $ASEG_PRESURF --no-rip"
fi
fi
