#!/bin/bash

#################
## Help
#################
Help ()
{
builtin echo "
AUTHOR: Benoît Verreman

LAST UPDATE: 2023-10-30

DESCRIPTION: 
Use ribbon and subcortical NIFTI files to compute pial surface.
Also based on previously created <subjid> folder using Freesurfer 7.4.1
Create a log 'report.sh'.
Create a folder 'outputs' with all the output files.

PREREQUISITE:
Export SUBJECTS_DIR and FREESURFER_HOME correctly
Launch Freesurfer 7.4.1 command: 
$ recon-all -s <subjid> -i <subject_image> -autorecon1 -autorecon2 -hires -parallel -openmp 4 -expert expert_file.txt
Prepare ribbon and subcortical NIFTI files (first step in the pipeline: convert)

EXAMPLES:
$ bash ribbon_edit_script.sh --subjid <subjid> --ribbon <ribbon-edit.nii.gz> --subcortical <subcortical.nii.gz>

$ bash ribbon_edit_script.sh -s <subjid> --pial --rh

PARAMETERS:
--subjid or -s: Relative or absolute path to subjid folder 
--ribbon or -r: Relative or absolute path to ribbon file
--subcortical or -c: Relative or absolute path to subcortical file

--help or -h: Print this string, and exit

TAG
--ribbons: Start with resizing RIBBON_EDIT and SUBCORTICAL
--bmask: Start with BRAIN_MASK
--maskT1: Start with T1_MASKED
--brain.finalsurfs: Start with skull-stripping up to BRAIN_FINALSURFS
--wm-bmask: Start the creation of WM_BMASK based on RIBBON_EDIT
--wm: Start from computing WM based on WM_BMASK
--orig: Start from computing orig surface based on wm from RIBBON_EDIT
--stats: Start from computing stats
--pial: Start from computing pial surface
--smooth: Start from smoothing pial surface

HEMI
--rh: Compute only right hemisphere surface
--lh: Compute only left hemisphere surface

--del: Reset outputs folder and report.sh script

"
}

#################
## Default global variables
#################
current_date_time=$(date)
TAG=1 # Start from beginning
HEMI=0 # Both hemispheres
OUTPUT_FOLDER="outputs"

#################
## Function to print both on terminal and on script report.sh
#################
Echo ()
{
    builtin echo "$@" | tee -a report.sh
}

#################
## Function to reset report.sh and outputs
#################
Create()
{
if [ ! -d "$OUTPUT_FOLDER" ]
then
Echo "#!/bin/bash"

cmd "Create $OUTPUT_FOLDER" \
"mkdir $OUTPUT_FOLDER;
mkdir $OUTPUT_FOLDER/transforms;
mkdir $OUTPUT_FOLDER/scripts;
mkdir $OUTPUT_FOLDER/surf;
mkdir $OUTPUT_FOLDER/mri;
mkdir $OUTPUT_FOLDER/label;"
fi
}
Delete()
{
rm -f report.sh

cmd "Reset $OUTPUT_FOLDER" \
"rm -r $OUTPUT_FOLDER;"
Create
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
## Manage flags
#################
unset -v SUBJID
unset -v RIBBON
unset -v SUBCORTICAL

VALID_ARGS=$(getopt -o s:r:c:h --long subjid:,ribbon:,subcortical:,help,convert,bmask,maskT1,brain.finalsurfs,wm-bmask,wm,orig,stats,pial,smooth,rh,lh,del -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -s | --subjid)
        SUBJID=$2
        shift 2
        ;;
    -r | --ribbon)
        RIBBON=$2
        shift 2
        ;;
    -c | --subcortical)
        SUBCORTICAL=$2
        shift 2
        ;;     
    -h | --help)
	Help
	exit 1
	;;
    --convert)
	TAG=2
	shift
	;;
    --bmask)
	TAG=3
	shift
	;;
    --maskT1)
	TAG=4
	shift
	;;	
    --brain.finalsurfs)
	TAG=5
	shift
	;;
    --wm-bmask)
	TAG=6
	shift
	;;
    --wm)
	TAG=7
	shift
	;;
    --orig)
	TAG=8
	shift
	;;
    --stats)
	TAG=9
	shift
	;;
    --pial)
	TAG=10
	shift
	;;
    --smooth)
	TAG=11
	shift
	;;
    --lh)
	HEMI=1
	shift
	;;
    --rh)
	HEMI=-1
	shift
	;;
    --del)
	Delete
	shift
	;;
    --) shift; 
        break 
        ;;
  esac
done

shift "$(( OPTIND - 1 ))"

# Test if user provided SUBJID
: ${SUBJID:?Missing argument --subjid or -s}

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
RIBBON_EDIT="$OUTPUT_FOLDER/mri/ribbon-edit.mgz"
SUBCORTICAL_EDIT="$OUTPUT_FOLDER/mri/subcortical-edit.mgz"

SUBCORTICAL_MASK="$OUTPUT_FOLDER/mri/subcortical-mask.mgz"
BRAIN_MASK="$OUTPUT_FOLDER/mri/brain-mask.mgz"

T1_MASKED="$OUTPUT_FOLDER/mri/T1-masked.mgz"

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

LH_CORTEX_LABEL="$OUTPUT_FOLDER/label/lh.cortex.label"
RH_CORTEX_LABEL="$OUTPUT_FOLDER/label/rh.cortex.label"
#LH_CORTEX_HIPAMYG_LABEL="$OUTPUT_FOLDER/label/lh.cortex+hipamyg.label"
#RH_CORTEX_HIPAMYG_LABEL="$OUTPUT_FOLDER/label/rh.cortex+hipamyg.label"

AUTODET_NEW_GW_STATS_LH="$OUTPUT_FOLDER/surf/autodet-new.gw.stats.lh.dat"
AUTODET_NEW_GW_STATS_RH="$OUTPUT_FOLDER/surf/autodet-new.gw.stats.rh.dat"
LH_RIBBON_EDIT_PIAL="$OUTPUT_FOLDER/surf/lh.ribbon_edit.pial"
RH_RIBBON_EDIT_PIAL="$OUTPUT_FOLDER/surf/rh.ribbon_edit.pial"

LH_RIBBON_EDIT_PIAL_SMOOTH="$OUTPUT_FOLDER/surf/lh.ribbon_edit.smooth.pial"
RH_RIBBON_EDIT_PIAL_SMOOTH="$OUTPUT_FOLDER/surf/rh.ribbon_edit.smooth.pial"

#################
## Indicate a new invocation of this script in report.sh
#################
Create #Test if outputs folder already exist, and if it does not, create one
Echo "
#*******************
#*******************
#*******************
# New invocation of ribbon_edit_script.sh

# Given subjid: $SUBJID"

#################
## Convert ribbon and subcortical
#################
if ((TAG<=2))
then
# Test if user provided RIBBON and SUBCORTICAL
: ${RIBBON:?Missing argument --ribbon or -r} ${SUBCORTICAL:?Missing argument --subcortical or -c}
Echo "# Given ribbon: $RIBBON"
Echo "# Given subcortical: $SUBCORTICAL"

cmd "Convert $RIBBON and $SUBCORTICAL for same dimensions" \
"mri_convert $RIBBON $RIBBON_EDIT -rt nearest -ns 1 --conform_min"

cmd "" \
"mri_convert $SUBCORTICAL $SUBCORTICAL_EDIT -rt nearest -ns 1 --conform_min"
fi

# Get corrected bmask of the whole brain
if ((TAG<=3))
then
cmd "Extract labels from $SUBCORTICAL_EDIT (Cerebellum, Medulla oblongata, Pons and Midbrain) into $SUBCORTICAL_MASK" \
"mri_extract_label $SUBCORTICAL_EDIT 5 15 29 30 33 34 $SUBCORTICAL_MASK"

cmd "Concatenate $RIBBON_EDIT with $SUBCORTICAL_MASK into $BRAIN_MASK" \
"mri_concat --i $RIBBON_EDIT --i $SUBCORTICAL_MASK --o $BRAIN_MASK --combine"
fi

#################
## Recompute brain.finalsurfs.mgz
#################
# Mask T1
if ((TAG<=4))
then
cmd "Mask $T1 with $BRAIN_MASK into $T1_MASKED" \
"mri_mask $T1 $BRAIN_MASK $T1_MASKED"
fi

# Recompute recon-all steps starting at EM Register up to brain.finalsurfs.mgz
if ((TAG<=5))
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
"mri_cc -aseg ../$ASEG_AUTO_NOCCSEG -o ../$ASEG_AUTO -lta $CC_UP $SUBJID" # Function "mri_cc" add "mri/" to $ASEG_AUTO_NOCCSEG and $ASEG_AUTO

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
	"mri_pretess $RIBBON_EDIT 3 $NORM $FILLED_PRETRESS_RH"

	cmd "Tessellate rh WM surf" \
	"mri_tessellate $FILLED_PRETRESS_RH 3 $RH_ORIG_NOFIX_PREDEC"

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
if ((TAG<=9))
then
if ((HEMI>=0))
then
	# Compute stats
	cmd "Computes stats for lh pial surface" \
	"mris_autodet_gwstats --o $AUTODET_NEW_GW_STATS_LH --i $BRAIN_FINALSURFS --wm $WM_EDITED --surf $LH_ORIG"
	
	# In order to improve pial surface, you can lower 'pial_border_low' to 20 
	# Change stats
	cmd "Change stats" \
	"ex -s -c '%s/^pial_border_low.*/pial_border_low   5/g|x' $AUTODET_NEW_GW_STATS_LH"	
	
	# Compute labels for pin-medial-wall
	cmd "Label2label for lh cortex" \
	"mri_label2label --label-cortex $LH_ORIG $ASEG_PRESURF 0 $LH_CORTEX_LABEL"
	
	# Compute labels to remove HIPOCAMPUS AND AMYGDALA from pial surface in mris_place_surface
	#cmd "Label2label for lh cortex" \
	#"mri_label2label --label-cortex $LH_ORIG $ASEG_PRESURF 1 $LH_CORTEX_HIPAMYG_LABEL"

fi

if ((HEMI<=0))
then
	# Compute stats
	cmd "Computes stats for rh pial surface" \
	"mris_autodet_gwstats --o $AUTODET_NEW_GW_STATS_RH --i $BRAIN_FINALSURFS --wm $WM_EDITED --surf $RH_ORIG"
	
	# In order to improve pial surface, you can lower 'pial_border_low' to 20 	
	# Change stats
	cmd "Change stats" \
	"ex -s -c '%s/^pial_border_low.*/pial_border_low   5/g|x' $AUTODET_NEW_GW_STATS_RH"	
	
	# Compute labels for pin-medial-wall
	cmd "Label2label for rh cortex" \
	"mri_label2label --label-cortex $RH_ORIG $ASEG_PRESURF 0 $RH_CORTEX_LABEL"

	# Compute labels to remove HIPOCAMPUS AND AMYGDALA from pial surface in mris_place_surface
	#cmd "Label2label for rh cortex" \
	#"mri_label2label --label-cortex $RH_ORIG $ASEG_PRESURF 1 $RH_CORTEX_HIPAMYG_LABEL"
fi
fi

#################
## Compute pial surface: mris_place_surface
#################
if ((TAG<=10))
then
if ((HEMI>=0))
then	
	cmd "Computes lh pial surface" \
	"mris_place_surface --i $LH_ORIG --o $LH_RIBBON_EDIT_PIAL --nsmooth 0 --adgws-in $AUTODET_NEW_GW_STATS_LH --pial --lh --repulse-surf $LH_ORIG --invol $BRAIN_FINALSURFS --threads 6 --white-surf $LH_ORIG --pin-medial-wall $LH_CORTEX_LABEL --seg $ASEG_PRESURF --no-rip"
fi

if ((HEMI<=0))
then
	cmd "Computes rh pial surface" \
	"mris_place_surface --i $RH_ORIG --o $RH_RIBBON_EDIT_PIAL --nsmooth 0 --adgws-in $AUTODET_NEW_GW_STATS_RH --pial --rh --repulse-surf $RH_ORIG --invol $BRAIN_FINALSURFS --threads 6 --white-surf $RH_ORIG --pin-medial-wall $RH_CORTEX_LABEL --seg $ASEG_PRESURF --no-rip"
fi
fi

#################
## Add smoothing to pial surface (mris_place_surface)
#################
if ((TAG<=11))
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
