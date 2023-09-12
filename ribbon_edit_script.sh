#!/bin/bash

# Global variables
current_date_time=$(date)
TAG=1
HEMI=0 # Both hemispheres
DEL=0 # Don't reset outputs folder
OUTPUT_FOLDER="outputs"
SUBJID="198451_hires"
SUBJECTS_DIR="/home/bverreman/Desktop/FS741_198451"
FREESURFER_HOME="/usr/local/freesurfer/7.4.1"

## Function to print both on terminal and on script report.sh
echo ()
{
    builtin echo "$@" | tee -a report.sh
}

## Function to both print and launch commands ($2) with a description ($1)
cmd () {
if [ -z "$1" ] #First variable is empty
then
	echo "
$2"
else
	echo "
#---------------------------------
#@# $1: $current_date_time

$2"
fi
eval $2
}

## Associate inputs given to script with a command or a starting TAG
if [ $# -ge 1 ]
then
    options="$@"
else
    options="empty"
fi

for opt in $options
do
case $opt in
	-ribbons)
	TAG=2
	;;

	-bmask)
	TAG=3
	;;

	-maskT1)
	TAG=4
	;;	

	-brain.finalsurfs)
	TAG=5
	;;

	-wm-bmask)
	TAG=6
	;;

	-wm)
	TAG=7
	;;

	-orig)
	TAG=8
	;;
	
	-pial)
	TAG=9
	;;
	
	-lh)
	HEMI=1
	;;
	
	-rh)
	HEMI=-1
	;;
	
	-del)
	DEL=1
	;;
	
	-help | --h)
	builtin echo "
Author: Benoît Verreman

Last update: 2023-09-12

Descrpition: 
Process your RIBBON_EDIT and RIBBON_SUBCORTICAL in subjid folder, previously created by 'recon-all -all -hires' from Freesurfer 7.4.1
Create a bash script report.sh in subjid folder with commands used by this script.
Create a folder OUTPUT_FOLDER with all the output files.

Parameters:
-ribbons: Start with resizing RIBBON_EDIT and RIBBON_SUBCORTICAL
-bmask: Start with BRAIN_BMASK
-maskT1: Start with T1_MASKED
-brain.finalsurfs: Start with skull-stripping up to BRAIN_FINALSURFS
-wm-bmask: Start the creation of WM_BMASK based on RIBBON_EDIT
-wm: Start from computing WM based on WM_BMASK
-orig: Start from computing orig surface based on wm from RIBBON_EDIT
-pial: Start from computing pial surface
-rh: Compute only right hemisphere surface
-lh: Compute only left hemisphere surface
-del: Reset folder output
-help or --h: Print a description of the script, and exit
*: Wrong parameter given

Example:
bash ribbon_edit_script.sh -pial -rh
"
	exit 1
	;;

	empty)
	;;
	
	*)
	builtin echo "
Wrong parameter given: $opt
See 'bash ribbon_edit_script.sh -help' or 'bash ribbon_edit_script.sh --h'
"
	exit 1
	;;
esac 
done

## Input files

RIBBON_EDIT_ORIGINAL="../198451_ribbon-edit.nii.gz"
RIBBON_SUBCORTICAL_ORIGINAL="../198451_cerebellum_and_subcortical_label.nii.gz"

T1="mri/T1.mgz"

TALAIRACH_WITH_SKULL="mri/transforms/talairach_with_skull.lta"
NU="mri/nu.mgz"
RB_ALL="$FREESURFER_HOME/average/RB_all_2020-01-02.gca"
CTRL_PTS="mri/ctrl_pts.mgz"
CC_UP="$SUBJECTS_DIR/$SUBJID/mri/transforms/cc_up.lta"

WM="mri/wm.mgz"

SUBCORTICALMASSLUT="$FREESURFER_HOME/SubCorticalMassLUT.txt"


## Output files

RIBBON_EDIT="$OUTPUT_FOLDER/mri/ribbon-edit.mgz"
RIBBON_SUBCORTICAL="$OUTPUT_FOLDER/mri/ribbon-subcortical.mgz"

SUBCORTICAL_MASK="$OUTPUT_FOLDER/mri/subcortical-mask.mgz"
BRAIN_MASK="$OUTPUT_FOLDER/mri/brain-mask.mgz"
BRAIN_BMASK="$OUTPUT_FOLDER/mri/brain-bmask.mgz"

T1_MASKED="$OUTPUT_FOLDER/mri/T1-masked.mgz"

BRAINMASK_AUTO="$OUTPUT_FOLDER/mri/brainmask.auto.mgz"
BRAINMASK="$OUTPUT_FOLDER/mri/brainmask.mgz"

TALAIRACH="$OUTPUT_FOLDER/transforms/talairach.lta"
NORM="$OUTPUT_FOLDER/mri/norm.mgz"
TALAIRACH_M3Z="$OUTPUT_FOLDER/transforms/talairach.m3z"
ASEG_AUTO_NOCCSEG="$OUTPUT_FOLDER/mri/aseg.auto_noCCseg.mgz"
ASEG_AUTO="$OUTPUT_FOLDER/mri/aseg.auto.mgz"
ASEG_PRESURF="$OUTPUT_FOLDER/mri/aseg.presurf.mgz"
BRAIN="$OUTPUT_FOLDER/mri/brain.mgz"
BRAIN_FINALSURFS="$OUTPUT_FOLDER/mri/brain.finalsurfs.mgz"

WM_BMASK="$OUTPUT_FOLDER/mri/wm-bmask.mgz"
WM_MASK="$OUTPUT_FOLDER/mri/wm-mask.mgz"
WM_CONCAT="$OUTPUT_FOLDER/mri/wm-concat.mgz"
WM_BMASK_250="$OUTPUT_FOLDER/mri/wm-bmask-250.mgz"
WM_ASEGEDIT="$OUTPUT_FOLDER/mri/wm-asegedit.mgz"
WM_EDITED="$OUTPUT_FOLDER/mri/wm.mgz" # Use this name for mri_fix_topology

FILLED_PRETRESS_LH="$OUTPUT_FOLDER/mri/filled_pretress_lh.mgz"
FILLED_PRETRESS_RH="$OUTPUT_FOLDER/mri/filled_pretress_rh.mgz"
LH_ORIG_NOFIX_PREDEC="$OUTPUT_FOLDER/surf/lh.orig.nofix.predec"
RH_ORIG_NOFIX_PREDEC="$OUTPUT_FOLDER/surf/rh.orig.nofix.predec"
LH_ORIG_NOFIX="$OUTPUT_FOLDER/surf/lh.orig.nofix"
RH_ORIG_NOFIX="$OUTPUT_FOLDER/surf/rh.orig.nofix"

AUTODET_NEW_GW_STATS_LH="$OUTPUT_FOLDER/surf/autodet-new.gw.stats.lh.dat"
AUTODET_NEW_GW_STATS_RH="$OUTPUT_FOLDER/surf/autodet-new.gw.stats.rh.dat"
LH_RIBBON_EDIT_PIAL="$OUTPUT_FOLDER/surf/lh.ribbon_edit.pial"
RH_RIBBON_EDIT_PIAL="$OUTPUT_FOLDER/surf/rh.ribbon_edit.pial"
LH_SMOOTHW_NOFIX="$OUTPUT_FOLDER/surf/lh.smoothwm.nofix"
RH_SMOOTHW_NOFIX="$OUTPUT_FOLDER/surf/rh.smoothwm.nofix"
LH_INFLATED_NOFIX="$OUTPUT_FOLDER/surf/lh.inflated.nofix"
RH_INFLATED_NOFIX="$OUTPUT_FOLDER/surf/rh.inflated.nofix"
LH_QSPHERE_NOFIX="$OUTPUT_FOLDER/surf/lh.qsphere.nofix"
RH_QSPHERE_NOFIX="$OUTPUT_FOLDER/surf/rh.qsphere.nofix"
LH_ORIG_PREMESH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.orig.premesh"
RH_ORIG_PREMESH="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.orig.premesh"
LH_ORIG="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/lh.orig"
RH_ORIG="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER/surf/rh.orig"


## Refresh report.sh and outputs if exist
if ((DEL))
then
rm -f report.sh

echo "#!/bin/bash"

cmd "Reset $OUTPUT_FOLDER" \
"rm -r $OUTPUT_FOLDER; \
mkdir $OUTPUT_FOLDER;
mkdir $OUTPUT_FOLDER/transforms;
mkdir $OUTPUT_FOLDER/scripts;
mkdir $OUTPUT_FOLDER/surf;
mkdir $OUTPUT_FOLDER/mri;
export SUBJECTS_DIR=$SUBJECTS_DIR"
fi

## Indicate a new invocation of this script in report.sh
cmd "
#*******************
#*******************
#*******************
# New invocation of ribbon_edit_script.sh" \
""

## Use ribbons
# Convert ribbons
if ((TAG<=2))
then
cmd "Convert $RIBBON_EDIT and $RIBBON_SUBCORTICAL for same dimensions" \
"mri_convert $RIBBON_EDIT_ORIGINAL $RIBBON_EDIT -rt nearest -ns 1 --conform_min"

cmd "" \
"mri_convert $RIBBON_SUBCORTICAL_ORIGINAL $RIBBON_SUBCORTICAL -rt nearest -ns 1 --conform_min"
fi

# Get corrected bmask of the whole brain
if ((TAG<=3))
then
cmd "Extract labels from $RIBBON_SUBCORTICAL (Cerebellum, Medulla oblongata, Pons and Midbrain) into $SUBCORTICAL_MASK" \
"mri_extract_label $RIBBON_SUBCORTICAL 5 15 29 30 33 34 $SUBCORTICAL_MASK"

cmd "Concatenate $RIBBON_EDIT with $SUBCORTICAL_MASK into $BRAIN_MASK" \
"mri_concat --i $RIBBON_EDIT --i $SUBCORTICAL_MASK --o $BRAIN_MASK --combine"

cmd "Binarize $BRAIN_MASK into $BRAIN_BMASK" \
"mri_binarize --i $BRAIN_MASK --o $BRAIN_BMASK --min 1"
fi

## Recompute brain.finalsurfs.mgz
# Mask T1
if ((TAG<=4))
then
cmd "Mask $T1 with $BRAIN_BMASK into $T1_MASKED" \
"mri_mask $T1 $BRAIN_BMASK $T1_MASKED"
fi

# Recompute recon-all steps starting at skull-stripping/watershed up to brain.finalsurfs.mgz
if ((TAG<=5))
then
cmd "Skull Strip" \
"mri_watershed -T1 -brain_atlas $FREESURFER_HOME/average/RB_all_withskull_2016-05-10.vc700.gca $TALAIRACH_WITH_SKULL $T1_MASKED $BRAINMASK_AUTO"

cmd "" \
"cp $BRAINMASK_AUTO $BRAINMASK"

cmd "EM Register: mri_em_register" \
"mri_em_register -uns 3 -mask $BRAINMASK $NU $RB_ALL $TALAIRACH"

cmd "CA Normalize: mri_ca_normalize" \
"mri_ca_normalize -c $CTRL_PTS -mask $BRAINMASK $NU $RB_ALL $TALAIRACH $NORM"

cmd "CA Register: mri_ca_register" \
"mri_ca_register -nobigventricles -T $TALAIRACH -align-after -mask $BRAINMASK $NORM $RB_ALL $TALAIRACH_M3Z"

cmd "Subcortical Segment: mri_ca_label" \
"mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align $NORM $TALAIRACH_M3Z $RB_ALL $ASEG_AUTO_NOCCSEG"

cmd "CC Segment: mri_cc" \
"mri_cc -aseg ../$ASEG_AUTO_NOCCSEG -o ../$ASEG_AUTO -lta $CC_UP $SUBJID" # Function "mri_cc" add "mri/" to $ASEG_AUTO_NOCCSEG and $ASEG_AUTO

cmd "Merge ASeg" \
"cp $ASEG_AUTO $ASEG_PRESURF"

cmd "Intensity Normalize" \
"mri_normalize -seed 1234 -mprage -noconform -aseg $ASEG_PRESURF -mask $BRAINMASK $NORM $BRAIN"

cmd "Mask BFS" \
"mri_mask -T 5 $BRAIN $BRAINMASK $BRAIN_FINALSURFS"
fi

## Compute wm.mgz : wm-bmask AND if(wm == 250 & wm-bmask != 0), then wm-mask = 250
# Extract white matter from ribbon-edit to create wm-bmask.mgz
if ((TAG<=6))
then
cmd "Extract WM from $RIBBON_EDIT" \
"mri_extract_label $RIBBON_EDIT 1 3 $WM_BMASK" #0/128 binary mask
fi

# Compute WM_EDIT based on BRAIN_FINALSURFS masked by WM_BMASK
if ((TAG<=7))
then
cmd "Concatenate $WM_BMASK with $WM into $WM_CONCAT" \
"mri_concat --i $WM_BMASK --i $WM --o $WM_CONCAT --sum" #ROI at 378 (128+250)

cmd "Binarize $WM_CONCAT at 251 into $WM_BMASK_250" \
"mri_binarize --i $WM_CONCAT --o $WM_BMASK_250 --match 378"

cmd "Replace 1 by 250 into $WM_BMASK_250" \
"mri_binarize --i $WM_BMASK_250 --o $WM_BMASK_250 --replace 1 250"

cmd "Mask $BRAIN with $WM_BMASK into $WM_MASK" \ # May also use $BRAIN_FINALSURFS
"mri_mask -T 5 $BRAIN $WM_BMASK $WM_MASK"

cmd "Concatenate $WM_MASK with $WM_BMASK_250 into $WM_ASEGEDIT" \
"mri_concat --i $WM_MASK --i $WM_BMASK_250 --o $WM_ASEGEDIT --max"

cmd "Pretess $WM_ASEGEDIT: Solve connectivity issue" \
"mri_pretess $WM_ASEGEDIT wm $NORM $WM_EDITED"
fi

## Compute ORIG: Don't need mri_fill, use ribbon-edit wm
if ((TAG<=8))
then
if ((HEMI>=0))
then
	# Compute directly ORIG_NOFIX
	cmd "Pretress lh WM from $RIBBON_EDIT" \
	"mri_pretess $RIBBON_EDIT 1 $NORM $FILLED_PRETRESS_LH"

	cmd "Tessellate lh WM surf" \
	"mri_tessellate $FILLED_PRETRESS_LH 1 $LH_ORIG_NOFIX_PREDEC"

	cmd "Extract main component lh WM surf" \
	"mris_extract_main_component $LH_ORIG_NOFIX_PREDEC $LH_ORIG_NOFIX_PREDEC"

	cmd "Remesh lh WM surf" \
	"mris_remesh --desired-face-area 0.5 --input $LH_ORIG_NOFIX_PREDEC --output $LH_ORIG_NOFIX"

	# Smooth 1
	cmd "Smooth lh WM surf" \
	"mris_smooth -nw -seed 1234 $LH_ORIG_NOFIX $LH_SMOOTHW_NOFIX"

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
	"mri_pretess $RIBBON_EDIT 3 $NORM $FILLED_PRETRESS_RH"

	cmd "Tessellate rh WM surf" \
	"mri_tessellate $FILLED_PRETRESS_RH 3 $RH_ORIG_NOFIX_PREDEC"

	cmd "Extract main component rh WM surf" \
	"mris_extract_main_component $RH_ORIG_NOFIX_PREDEC $RH_ORIG_NOFIX_PREDEC"

	cmd "Remesh rh WM surf" \
	"mris_remesh --desired-face-area 0.5 --input $RH_ORIG_NOFIX_PREDEC --output $RH_ORIG_NOFIX"

	# Smooth 1
	cmd "Smooth rh WM surf" \
	"mris_smooth -nw -seed 1234 $RH_ORIG_NOFIX $RH_SMOOTHW_NOFIX"

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

## Compute pial surface
if ((TAG<=9))
then
if ((HEMI>=0))
then
	# Compute stats
	cmd "Computes stats for lh pial surface" \
	"mris_autodet_gwstats --o $AUTODET_NEW_GW_STATS_LH --i $BRAIN_FINALSURFS --wm $WM_EDITED --surf $LH_ORIG"
	
	cmd "Computes lh pial surface" \
	"mris_place_surface --i $LH_ORIG --o $LH_RIBBON_EDIT_PIAL --nsmooth 1 --adgws-in $AUTODET_NEW_GW_STATS_LH --pial --lh --repulse-surf $LH_ORIG --invol $BRAIN_FINALSURFS --seg $ASEG_PRESURF --threads 6 --wm $WM_EDITED --white-surf $LH_ORIG"
fi

if ((HEMI<=0))
then
	# Compute stats
	cmd "Computes stats for rh pial surface" \
	"mris_autodet_gwstats --o $AUTODET_NEW_GW_STATS_RH --i $BRAIN_FINALSURFS --wm $WM_EDITED --surf $RH_ORIG"
	
	cmd "Computes rh pial surface" \
	"mris_place_surface --i $RH_ORIG --o $RH_RIBBON_EDIT_PIAL --nsmooth 1 --adgws-in $AUTODET_NEW_GW_STATS_RH --pial --rh --repulse-surf $RH_ORIG --invol $BRAIN_FINALSURFS --seg $ASEG_PRESURF --threads 6 --wm $WM_EDITED --white-surf $RH_ORIG"
fi
fi
