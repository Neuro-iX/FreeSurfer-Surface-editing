#!/bin/bash

#################
## Help
#################
Help ()
{
builtin echo "
AUTHOR: Benoît Verreman

LAST UPDATE: 2026-06-24

DESCRIPTION:
Use a ground truth aseg (in FreeSurfer label space) and a T1 image to recompute pial surface,
based on Freesurfer 7.4.1 intermediate files.
The aseg must already contain FreeSurfer label values (WM=2/41, GM=3/42, subcortical, etc.)
such as those produced by keep_main_components.py.
Ribbon (labels 2,3,41,42) is derived automatically from the aseg.
Create a log 'report.sh'.
Create a folder 'outputs' with all the output files.

PREREQUISITES:
*Freesurfer variables:
Export SUBJECTS_DIR and FREESURFER_HOME correctly

*Python:
pip install nibabel scipy
OR:
conda create --name env_ribbon_edit_script python=3.10 nibabel scipy -c conda-forge

*Aseg file (FreeSurfer label space):
  WM:          2 (lh), 41 (rh)
  GM:          3 (lh), 42 (rh)
  Subcortical: 4,5,7,8,10-18,24,26,28,43-47,49-54,58,60,77,85,251-255...

*Freesurfer output folder:
If you want to run recon-all from inside the script, use -i.
Otherwise run beforehand:
  recon-all -s <subjid>_freesurfer -i <T1> -autorecon1 -autorecon2 -autorecon3 -hires -parallel -openmp 4

EXAMPLES:
$ export SUBJECTS_DIR=\`pwd\`
$ bash /home/bverreman/Documents/FreeSurfer-Surface-editing/aseg_surface_script.sh -i /project/hippocampus/common/bverreman/nnUNet/nnUNet/nnUNet_raw/Dataset110_Brain9/imagesTs_not_masked_yet/100307_0000.nii.gz -s 100307 -a /project/hippocampus/common/bverreman/nnUNet/nnUNet/nnUNet_raw/Dataset110_Brain9/imagesTs_not_masked_yet/100307.nii.gz

$ bash aseg_surface_script.sh -s 100307 -t 4 -r
#restart from orig surface (TAG 4), right hemisphere only

$ bash aseg_surface_script.sh -f /data/subjects -i * -a *
#whole folder mode: T1 detected by '*T1*', aseg detected by '*aseg*'

PARAMETERS:

HELP
-h: Print this string, and exit

INPUT FILES
-i: Relative or absolute path to T1w image file (also triggers recon-all)
-s: Subject ID / subfolder name (REQUIRED)
-a: Relative or absolute path to ground truth aseg (FreeSurfer label space)

TAG
-t 0: (aseg)          Pad ASEG, extract RIBBON_EDIT (labels 2,3,41,42)
-t 1: (bmask)         Brain mask from ASEG, mask T1
-t 2: (aseg.presurf)  Create ASEG_PRESURF and ASEG_PRESURF_WO_SUBC
-t 3: (wm)            WM from RIBBON_EDIT, pretess
-t 4: (orig)          Orig surface tessellation and topology fix
-t 5: (bf-edit)       Edit brain.finalsurfs with GM from RIBBON_EDIT
-t 6: (stats)         Stats, cortex label, sphere, aparc
-t 7: (pial)          Pial surface (mris_place_surface)
-t 8: (non-hemi)      Non-hemisphere files (copy, denoise, segment)
-t 9: (autorecon3)    Curvature and thickness stats
-t 10: (aparc+aseg)   Cortical ribbon, aparc+aseg, wmparc, segstats
-t 11: (bonus)        Extended parcellations (BA_exvivo, vpnl)

VALUES
-p: Value of PIAL_BORDER_LOW (default: 5)

HEMI
-r: Right hemisphere only
-l: Left hemisphere only

RESET
-k: Reset outputs folder and report.sh

WHOLE FOLDER
-f: Run on every subfolder in the given folder.
    Each subfolder must contain:
      - an aseg file (*aseg*)
      - a T1 image (*T1*) if -i * is used
"
}

#################
## Default global variables
#################
TAG=-1
HEMI=-1
FS=0
MULTICASE=0
OUTPUT_FOLDER="outputs"
# Subcortical labels replaced by WM in ASEG_PRESURF_WO_SUBC
# LH: lat.vent(4) thalamus(10) caudate(11) putamen(12) pallidum(13) accumbens(26) ventralDC(28)
# RH: lat.vent(43) thalamus(49) caudate(50) putamen(51) pallidum(52) accumbens(58) ventralDC(60)
LABELS_SUBCORTICAL_WM_SRC="4 10 11 12 13 26 28 43 49 50 51 52 58 60"
LABELS_SUBCORTICAL_WM_DST="2  2  2  2  2  2  2  41 41 41 41 41 41 41"
declare -a H=("lh" "rh")
declare -a LABEL_RIBBON_WM=("2" "41")
declare -a LABEL_RIBBON_GM=("3" "42")
declare -i N_PARALLEL_COMPUTING=3

CHANGE_AUTODET=1
PIAL_BORDER_LOW=5

#################
## Manage flags
#################
string_arguments=""

unset -v IMAGE
unset -v SUBJID
unset -v ASEG

VALID_ARGS="i:s:a:t:p:f:hlrk"

while getopts ${VALID_ARGS} opt; do
  case ${opt} in
    i)
        IMAGE=${OPTARG}
        FS=1
        string_arguments+="-i ${OPTARG} "
        ;;
    s)
        SUBJID=${OPTARG}
        string_arguments+="-s ${OPTARG} "
        ;;
    a)
        ASEG=${OPTARG}
        string_arguments+="-a ${OPTARG} "
        ;;
    t)
        TAG=${OPTARG}
        string_arguments+="-t ${OPTARG} "
        ;;
    p)
        PIAL_BORDER_LOW=${OPTARG}
        string_arguments+="-p ${OPTARG} "
        ;;
    f)
        SUBJECTS_DIR=${OPTARG}
        MULTICASE=1
        string_arguments+="-f ${OPTARG} "
        ;;
    h)
        Help
        exit 1
        ;;
    l)
        HEMI=0
        string_arguments+="-l "
        ;;
    r)
        HEMI=1
        string_arguments+="-r "
        ;;
    k)
        Delete
        string_arguments+="-k "
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

#################
## Remove trailing slash from SUBJECTS_DIR
#################
export var="${SUBJECTS_DIR: -1}"
if [[ "$var" == "/" ]]; then
export SUBJECTS_DIR="${SUBJECTS_DIR:0:-1}"
fi


#################
## Main function
#################
main()
{

: ${SUBJID:?Missing argument -s}

export var="${SUBJID: -1}"
if [[ "$var" == "/" ]]; then export SUBJID="${SUBJID:0:-1}"; fi
export var="${SUBJID:0:1}"
if [[ "$var" == "/" ]]; then export SUBJID="${SUBJID:1}"; fi

O="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER"

#################
## Logging helper
#################
Echo ()
{
    builtin echo "$@" | tee -a $O/report.sh
}

#################
## Folder creation
#################
CreateFolders()
{
echo "Create or confirm existence of $SUBJECTS_DIR/$SUBJID"
mkdir -p $SUBJECTS_DIR/$SUBJID

echo "Create or confirm existence of $O and subfolders"
mkdir -p $O $O/scripts $O/surf $O/mri $O/mri/orig $O/mri/transforms $O/label $O/stats $O/tmp $O/touch $O/trash
}

CreateScripts()
{
if [ ! -f "$O/report.sh" ]
then
touch $O/report.sh
Echo "#!/bin/bash"
fi

script_nifti_padding
script_extract_labels
script_brain_finalsurfs_edit
script_remap_labels
script_expert_file
}

Delete()
{
cmd "Reset $O" \
"rm -r $O;"
}

#################
## nifti_padding.py
#################
script_nifti_padding()
{
if [ ! -f "$O/nifti_padding.py" ]
then
cat > $O/nifti_padding.py <<'EOF'
#!/usr/bin/env python3
import os
import numpy as np
import nibabel as nib
import nibabel.processing
import scipy.ndimage
import sys

def padding(img, new_name):
    d = max(img.header.get_data_shape())
    new_img = nib.processing.conform(img, out_shape=(d, d, d), voxel_size=img.header.get_zooms(), order=0, cval=0, orientation='RAS', out_class=None)
    nib.save(new_img, new_name)

img_in = sys.argv[1]
img_padded = sys.argv[2]

if not os.path.isfile(img_in):
    raise FileNotFoundError("The following path doesn't exist: " + img_in)
else:
    img = nib.load(img_in)

padding(img, img_padded)
EOF
fi
}

#################
## extract_labels.py
## Extracts specified labels from a segmentation, preserving label values.
## Usage: python extract_labels.py input.mgz output.mgz 2 3 41 42
#################
script_extract_labels()
{
if [ ! -f "$O/extract_labels.py" ]
then
cat > $O/extract_labels.py <<'EOF'
#!/usr/bin/env python3
import sys
import numpy as np
import nibabel as nib

path_in  = sys.argv[1]
path_out = sys.argv[2]
labels   = [int(l) for l in sys.argv[3:]]

img  = nib.load(path_in)
data = img.get_fdata().astype(np.int32)
out  = np.where(np.isin(data, labels), data, 0)
nib.save(nib.Nifti1Image(out, img.affine, img.header), path_out)
print(f"Extracted labels {labels} from {path_in} -> {path_out}")
EOF
fi
}

#################
## brain-finalsurfs-edit.py  (HA labels from aseg: 17,18,53,54)
#################
script_brain_finalsurfs_edit()
{
if [ ! -f "$O/brain-finalsurfs-edit.py" ]
then
cat > $O/brain-finalsurfs-edit.py <<'EOF'
import os
import numpy as np
import nibabel as nib
import sys
import copy

path_bf   = sys.argv[1]  # brain.finalsurfs (uniform WM 110)
path_gmbm = sys.argv[2]  # GM binary mask at 128
path_out  = sys.argv[3]  # output brain.finalsurfs (edited)
path_out2 = sys.argv[4]  # output brain.finalsurfs (edited, HA zeroed)
path_ha   = sys.argv[5]  # aseg (hippocampus=17/53, amygdala=18/54)

HA_LABELS = [17, 18, 53, 54]

if not os.path.isfile(path_bf):
    raise FileNotFoundError(path_bf)
img_bf = nib.load(path_bf)
data_bf = img_bf.get_fdata()
(a,b,c) = img_bf.header.get_data_shape()

if not os.path.isfile(path_gmbm):
    raise FileNotFoundError(path_gmbm)
img_gmbm = nib.load(path_gmbm)
data_gmbm = img_gmbm.get_fdata()
(e,f,g) = img_gmbm.header.get_data_shape()

if (a!=e or b!=f or c!=g):
    sys.exit("Shape mismatch between brain.finalsurfs and GM mask")

if not os.path.isfile(path_ha):
    raise FileNotFoundError(path_ha)
img_ha = nib.load(path_ha)
data_ha = img_ha.get_fdata()

motion = np.transpose(np.indices((3,3,3)) - 1).reshape(-1, 3)

data_bf_new  = copy.deepcopy(data_bf)
data_bf_new2 = copy.deepcopy(data_bf)

for x in range(a):
    for y in range(b):
        for z in range(c):
            if int(data_gmbm[x,y,z]) == 128:
                val = data_bf[x,y,z]
                n_coordinates = motion + [[x, y, z]]
                mean = 0
                nn = 0
                for (k,n,m) in n_coordinates:
                    if int(data_gmbm[k,n,m]) == 128:
                        mean += data_bf[k,n,m]
                        nn += 1
                res = 80.0/(mean/nn)*val
                data_bf_new[x,y,z]  = res
                data_bf_new2[x,y,z] = res
            if int(data_ha[x,y,z]) in HA_LABELS:
                data_bf_new2[x,y,z] = 0

nib.save(nib.Nifti1Image(data_bf_new,  img_bf.affine.copy()), path_out)
nib.save(nib.Nifti1Image(data_bf_new2, img_bf.affine.copy()), path_out2)
EOF
fi
}

#################
## remap_labels.py
#################
script_remap_labels()
{
if [ ! -f "$O/remap_labels.py" ]
then
cat > $O/remap_labels.py <<'EOF'
#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import nibabel as nib

def parse_label_string(label_str):
    try:
        values = [float(v) for v in label_str[0].split()]
        return [int(v) if v == int(v) else v for v in values]
    except ValueError:
        print(f"Error: Could not parse label string '{label_str[0]}'.")
        sys.exit(1)

def remap_labels(data, src_labels, dst_labels, remove_unlisted=False):
    output = np.zeros_like(data) if remove_unlisted else data.copy()
    for src, dst in zip(src_labels, dst_labels):
        mask = (data == src)
        n = np.sum(mask)
        if n == 0:
            print(f"  Warning: Source label {src} not found.")
        else:
            output[mask] = dst
            print(f"  {src:>10.0f}  -->  {dst:<10.0f}  ({n} voxels)")
    return output

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--src", type=str, nargs=1, required=True)
parser.add_argument("--dst", type=str, nargs=1, required=True)
parser.add_argument("--remove-unlisted", action="store_true", default=False)
args = parser.parse_args()

src_labels = parse_label_string(args.src)
dst_labels = parse_label_string(args.dst)

if len(src_labels) != len(dst_labels):
    print("Error: --src and --dst must have same length.")
    sys.exit(1)

print(f"Loading: {args.input}")
img = nib.load(args.input)
data = img.get_fdata().astype(np.int32)
out = remap_labels(data, src_labels, dst_labels, args.remove_unlisted)
nib.save(nib.Nifti1Image(out, img.affine, img.header), args.output)
print(f"Saved: {args.output}")
EOF
fi
}

#################
## expert_file.txt
#################
script_expert_file()
{
if [ ! -f "$SUBJECTS_DIR/expert_file.txt" ]
then
cat > $SUBJECTS_DIR/expert_file.txt <<'EOF'
mris_inflate -n 30
EOF
fi
}

#################
## Command runner: prints + logs + executes
#################
cmd () {
if [ -z "$1" ]
then
    Echo "
$2"
else
    Echo "
#---------------------------------
#@# $1: $(date)

$2"
fi
eval $2
if [ $? -ne 0 ]; then
  Echo "#---------------------------------
ERROR DETECTED
#---------------------------------
"
  exit 1
fi
}

#################
## FreeSurfer input files (from recon-all run)
#################
IMAGE_ORIG_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/orig/001.mgz"
RAWAVG_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/rawavg.mgz"
T1_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/T1.mgz"
NORM_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/norm.mgz"
BRAIN_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/brain.mgz"
BRAIN_FINALSURFS_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/brain.finalsurfs.mgz"
BRAIN_FINALSURFS_MANEDIT_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/brain.finalsurfs.manedit.mgz"
WM_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/wm.mgz"
ORIG_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/orig.mgz"
TALAIRACH_XFM_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/transforms/talairach.xfm"
BRAINMASK_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/brainmask.mgz"

declare -a CD_APARC_ATLAS=("$FREESURFER_HOME/average/lh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs" "$FREESURFER_HOME/average/rh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs")
declare -a DKT_APARC_ATLAS=("$FREESURFER_HOME/average/lh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs" "$FREESURFER_HOME/average/rh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs")
WMPARC_STATS_LUT="$FREESURFER_HOME/WMParcStatsLUT.txt"
ASEG_STATS_LUT="$FREESURFER_HOME/ASegStatsLUT.txt"
FSAVERAGE="$FREESURFER_HOME/subjects/fsaverage"
COLORTABLE_VPNL_TXT="$FREESURFER_HOME/average/colortable_vpnl.txt"
COLORTABLE_BA_TXT="$FREESURFER_HOME/average/colortable_BA.txt"

declare -a FOLDING_ATLAS_ACFB40=("$FREESURFER_HOME/average/lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif" "$FREESURFER_HOME/average/rh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif")
declare -a DKAPARC_ATLAS_ACFB40=("$FREESURFER_HOME/average/lh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs" "$FREESURFER_HOME/average/rh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs")

#################
## Output files
#################
IMAGE_PADDED="$SUBJECTS_DIR/$SUBJID/RS_image-padded.mgz"
IMAGE_ORIG="$O/mri/orig/001.mgz"

ASEG_PADDED="$O/mri/RS_aseg-padded.mgz"

# Ribbon derived from aseg (labels 2, 3, 41, 42 only)
RIBBON_EDIT="$O/mri/ribbon-edit.mgz"
RIBBON_WO_EDIT="$O/mri/ribbon-wo-edit.mgz"  # same file: no HA editing needed

BRAIN="$O/mri/brain.mgz"
BRAIN_FINALSURFS="$O/mri/brain.finalsurfs.mgz"
BRAIN_FINALSURFS_MANEDIT="$O/mri/brain.finalsurfs.manedit.mgz"

BRAIN_MASK="$O/mri/RS_brain-mask.mgz"

T1="$O/mri/T1.mgz"
T1_MASKED="$O/mri/RS_T1-masked.mgz"

NORM="$O/mri/norm.mgz"

ASEG_PRESURF="$O/mri/aseg.presurf.mgz"
ASEG_PRESURF_WO_SUBC="$O/mri/RS_aseg.presurf_wo_subc.mgz"

WM_BMASK_ALL="$O/mri/RS_wm-bmask.mgz"
WM_MASK="$O/mri/RS_wm-mask.mgz"
WM_CONCAT="$O/mri/RS_wm-concat.mgz"
WM_BMASK_250="$O/mri/RS_wm-bmask-250.mgz"
WM_ASEGEDIT="$O/mri/RS_wm-asegedit.mgz"
WM="$O/mri/wm.mgz"

BRAINMASK="$O/mri/brainmask.mgz"

declare -a FILLED_PRETRESS=("$O/mri/RS_filled_pretress_lh.mgz" "$O/mri/RS_filled_pretress_rh.mgz")
declare -a ORIG_NOFIX_PREDEC=("$O/surf/lh.orig.nofix.predec" "$O/surf/rh.orig.nofix.predec")
declare -a ORIG_NOFIX=("$O/surf/lh.orig.nofix" "$O/surf/rh.orig.nofix")
declare -a SMOOTHW_NOFIX=("$O/surf/lh.smoothwm.nofix" "$O/surf/rh.smoothwm.nofix")
declare -a INFLATED_NOFIX=("$O/surf/lh.inflated.nofix" "$O/surf/rh.inflated.nofix")
declare -a QSPHERE_NOFIX=("$O/surf/lh.qsphere.nofix" "$O/surf/rh.qsphere.nofix")
declare -a ORIG_PREMESH=("$O/surf/lh.orig.premesh" "$O/surf/rh.orig.premesh")
declare -a ORIG=("$O/surf/lh.orig" "$O/surf/rh.orig")
declare -a INFLATED=("$O/surf/lh.inflated" "$O/surf/rh.inflated")
declare -a SMOOTHW=("$O/surf/lh.smoothwm" "$O/surf/rh.smoothwm")
declare -a CURV=("$O/surf/lh.curv" "$O/surf/rh.curv")
declare -a CORTEX_LABEL=("$O/label/lh.cortex.label" "$O/label/rh.cortex.label")
declare -a CORTEX_HIPAMYG_LABEL=("$O/label/lh.cortex+hipamyg.label" "$O/label/rh.cortex+hipamyg.label")
declare -a GM_BMASK=("$O/mri/RS_gm-bmask_lh.mgz" "$O/mri/RS_gm-bmask_rh.mgz")
declare -a WM_BMASK=("$O/mri/RS_wm-bmask_lh.mgz" "$O/mri/RS_wm-bmask_rh.mgz")
declare -a BMASK=("$O/mri/RS_bmask_lh.mgz" "$O/mri/RS_bmask_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB=("$O/mri/RS_brain.finalsurfs_no_cereb_lh.mgz" "$O/mri/RS_brain.finalsurfs_no_cereb_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110=("$O/mri/RS_brain.finalsurfs_no_cereb_uniform_wm_110_lh.mgz" "$O/mri/RS_brain.finalsurfs_no_cereb_uniform_wm_110_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80=("$O/mri/RS_brain.finalsurfs_no_cereb_uniform_gm_80_lh.mgz" "$O/mri/RS_brain.finalsurfs_no_cereb_uniform_gm_80_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB_EDITED=("$O/mri/RS_brain.finalsurfs_no_cereb_edited_lh.mgz" "$O/mri/RS_brain.finalsurfs_no_cereb_edited_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB_EDITED2=("$O/mri/RS_brain.finalsurfs_no_cereb_edited2_lh.mgz" "$O/mri/RS_brain.finalsurfs_no_cereb_edited2_rh.mgz")
declare -a AUTODET_NEW_GW_STATS=("$O/surf/autodet-new.gw.stats.lh.dat" "$O/surf/autodet-new.gw.stats.rh.dat")
declare -a RIBBON_EDIT_PIAL=("$O/surf/RS_lh.ribbon_edit.pial" "$O/surf/RS_rh.ribbon_edit.pial")
declare -a CURV=("$O/surf/lh.curv" "$O/surf/rh.curv")
declare -a AREA=("$O/surf/lh.area" "$O/surf/rh.area")
declare -a PIAL_CURV=("$O/surf/lh.curv.pial" "$O/surf/rh.curv.pial")
declare -a PIAL_AREA=("$O/surf/lh.area.pial" "$O/surf/rh.area.pial")
declare -a THICKNESS=("$O/surf/lh.thickness" "$O/surf/rh.thickness")
declare -a CURV_STATS=("$O/stats/lh.curv.stats" "$O/stats/rh.curv.stats")
declare -a WHITE=("$O/surf/lh.white" "$O/surf/rh.white")
declare -a PIAL=("$O/surf/lh.pial" "$O/surf/rh.pial")
declare -a SPHERE=("$O/surf/lh.sphere" "$O/surf/rh.sphere")
declare -a SPHERE_REG=("$O/surf/lh.sphere.reg" "$O/surf/rh.sphere.reg")

RAWAVG="$O/mri/rawavg.mgz"
RAWAVG_MASKED="$O/mri/RS_rawavg-masked.mgz"
ORIG_VOLUME="$O/mri/orig.mgz"
ORIG_MASKED="$O/mri/RS_orig-masked.mgz"

declare -a CD_APARC_ANNOT=("$O/label/lh.aparc.a2009s.annot" "$O/label/rh.aparc.a2009s.annot")
declare -a DKT_APARC_ANNOT=("$O/label/lh.aparc.DKTatlas.annot" "$O/label/rh.aparc.DKTatlas.annot")

RIBBON_NEW="$O/mri/ribbon.mgz"
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
declare -a HOC1_MPM_VPNL_LABEL=("$O/label/lh.h0c1.mpm.vpnl.label" "$O/label/rh.h0c1.mpm.vpnl.label")
declare -a HOC2_MPM_VPNL_LABEL=("$O/label/lh.h0c2.mpm.vpnl.label" "$O/label/rh.h0c2.mpm.vpnl.label")
declare -a HOC3V_MPM_VPNL_LABEL=("$O/label/lh.h0c3v.mpm.vpnl.label" "$O/label/rh.h0c3v.mpm.vpnl.label")
declare -a HOC4V_MPM_VPNL_LABEL=("$O/label/lh.h0c4v.mpm.vpnl.label" "$O/label/rh.h0c4v.mpm.vpnl.label")
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
TALAIRACH_XFM="$O/mri/transforms/talairach.xfm"
ANTSDN_BRAIN="$O/mri/antsdn.brain.mgz"
WM_SEG="$O/mri/wm.seg.mgz"

#################
## Initialize folders and scripts
#################
CreateFolders
CreateScripts

Echo "
#*******************
#*******************
#*******************
# New invocation: $(date)

bash $(basename "$0") $string_arguments

# Given subjid: $SUBJID"

#################
## TAG FS: recon-all -autorecon1 -autorecon2 -autorecon3
#################
if ((FS == 1 && TAG < 0))
then
: ${IMAGE:?Missing argument -i}
Echo "# Given image: $IMAGE"

if [ -d "$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer" ]
then
    Echo "Do not re-run FreeSurfer on same SUBJID: $SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer"
    exit 1
fi

cmd "Pad T1 to cubic RAS" \
"python $O/nifti_padding.py $IMAGE $IMAGE_PADDED"

cmd "Temporarily change SUBJECTS_DIR to SUBJECTS_DIR/SUBJID" \
"export SUBJECTS_DIR=$SUBJECTS_DIR/$SUBJID"

cmd "recon-all autorecon1+2+3 on $IMAGE_PADDED" \
"recon-all -autorecon1 -autorecon2 -autorecon3 -s ${SUBJID}_freesurfer -i $IMAGE_PADDED -hires -parallel -openmp 4 -expert expert_file.txt -xopts-overwrite -cw256"

cmd "Restore SUBJECTS_DIR" \
"export SUBJECTS_DIR=$(dirname $SUBJECTS_DIR)"
fi

#################
## TAG 0: Pad aseg, extract ribbon (labels 2,3,41,42), set RIBBON_WO_EDIT = RIBBON_EDIT
#################
if ((TAG<=0))
then
: ${ASEG:?Missing argument -a}
Echo "# Given aseg: $ASEG"

cmd "Pad ASEG to cubic RAS" \
"python $O/nifti_padding.py $ASEG $ASEG_PADDED"

cmd "Conform-min ASEG (nearest, no smoothing)" \
"mri_convert $ASEG_PADDED $ASEG_PADDED -rt nearest -ns 1 --conform_min"

cmd "Extract ribbon labels (2,3,41,42) from ASEG into $RIBBON_EDIT" \
"python $O/extract_labels.py $ASEG_PADDED $RIBBON_EDIT 2 3 41 42"

# No HA editing: RIBBON_WO_EDIT == RIBBON_EDIT (same file)
cmd "Copy $RIBBON_EDIT to $RIBBON_WO_EDIT (no HA editing needed)" \
"cp $RIBBON_EDIT $RIBBON_WO_EDIT"
fi

#################
## TAG 1: Brain mask (all nonzero in ASEG), mask T1
#################
if ((TAG<=1))
then
cmd "Binarize ASEG (all nonzero) into $BRAIN_MASK" \
"mri_binarize --i $ASEG_PADDED --o $BRAIN_MASK --min 1 --replace 1 128"

cmd "Copy $T1_FS into $T1" \
"cp $T1_FS $T1"

cmd "Mask $T1 with $BRAIN_MASK into $T1_MASKED" \
"mri_mask $T1 $BRAIN_MASK $T1_MASKED"
fi

#################
## TAG 2: ASEG_PRESURF and ASEG_PRESURF_WO_SUBC
#################
if ((TAG<=2))
then
cmd "Convert padded ASEG to ASEG_PRESURF (MGZ)" \
"mri_convert $ASEG_PADDED $ASEG_PRESURF"

# ASEG_PRESURF_WO_SUBC: replace subcortical structures that lie within WM territory with WM
# LH subcortical (lat.vent, thalamus, caudate, putamen, pallidum, accumbens, ventralDC) -> 2
# RH subcortical -> 41
cmd "Create $ASEG_PRESURF_WO_SUBC: replace subcortical WM-region labels with WM" \
"python $O/remap_labels.py $ASEG_PRESURF $ASEG_PRESURF_WO_SUBC \
  --src '${LABELS_SUBCORTICAL_WM_SRC}' \
  --dst '${LABELS_SUBCORTICAL_WM_DST}'"
fi

#################
## TAG 3: WM from RIBBON_EDIT, pretess
#################
if ((TAG<=3))
then
cmd "Extract WM from $RIBBON_WO_EDIT (labels 2 and 41) -> $WM_BMASK_ALL" \
"mri_extract_label $RIBBON_WO_EDIT ${LABEL_RIBBON_WM[0]} ${LABEL_RIBBON_WM[1]} $WM_BMASK_ALL"

cmd "Copy $WM_FS into $WM" \
"cp $WM_FS $WM"

cmd "Concatenate $WM_BMASK_ALL and $WM -> $WM_CONCAT (sum: ribbon=128, FS-wm=250)" \
"mri_concat --i $WM_BMASK_ALL --i $WM --o $WM_CONCAT --sum"

cmd "Extract intersection (128+250=378) -> $WM_BMASK_250" \
"mri_binarize --i $WM_CONCAT --o $WM_BMASK_250 --match 378"

cmd "Replace 1 by 250 in $WM_BMASK_250" \
"mri_binarize --i $WM_BMASK_250 --o $WM_BMASK_250 --replace 1 250"

cmd "Copy $BRAIN_FS into $BRAIN" \
"cp $BRAIN_FS $BRAIN"

cmd "Mask $BRAIN with $WM_BMASK_ALL -> $WM_MASK" \
"mri_mask -T 5 $BRAIN $WM_BMASK_ALL $WM_MASK"

cmd "Concatenate $WM_MASK and $WM_BMASK_250 -> $WM_ASEGEDIT (max)" \
"mri_concat --i $WM_MASK --i $WM_BMASK_250 --o $WM_ASEGEDIT --max"

cmd "Copy $NORM_FS into $NORM" \
"cp $NORM_FS $NORM"

cmd "Pretess $WM_ASEGEDIT: fix WM connectivity" \
"mri_pretess $WM_ASEGEDIT wm $NORM $WM"
fi

#################
## TAG 4: Orig surface (tessellate, topology fix, remesh)
#################
if ((TAG<=4))
then
for (( i=0; i<2; i++ ));
do
    if ((HEMI>=0 && HEMI!=i)); then continue; fi

    cmd "${H[$i]} Pretess WM from $RIBBON_WO_EDIT" \
    "mri_pretess $RIBBON_WO_EDIT ${LABEL_RIBBON_WM[$i]} $NORM_FS ${FILLED_PRETRESS[$i]}"

    cmd "${H[$i]} Tessellate WM surface" \
    "mri_tessellate ${FILLED_PRETRESS[$i]} ${LABEL_RIBBON_WM[$i]} ${ORIG_NOFIX_PREDEC[$i]}"

    cmd "${H[$i]} Extract main component" \
    "mris_extract_main_component ${ORIG_NOFIX_PREDEC[$i]} ${ORIG_NOFIX_PREDEC[$i]}"

    cmd "${H[$i]} Remesh (desired-face-area 0.5)" \
    "mris_remesh --desired-face-area 0.5 --input ${ORIG_NOFIX_PREDEC[$i]} --output ${ORIG_NOFIX[$i]}"

    cmd "${H[$i]} Smooth 1" \
    "mris_smooth -n 1 -nw -seed 1234 ${ORIG_NOFIX[$i]} ${SMOOTHW_NOFIX[$i]}"

    cmd "${H[$i]} Inflate 1" \
    "mris_inflate -no-save-sulc -n 30 ${SMOOTHW_NOFIX[$i]} ${INFLATED_NOFIX[$i]}"

    cmd "${H[$i]} Sphere 1 (qsphere)" \
    "mris_sphere -q -p 6 -a 128 -seed 1234 ${INFLATED_NOFIX[$i]} ${QSPHERE_NOFIX[$i]}"

    cmd "${H[$i]} Copy $BRAIN_FS for mris_fix_topology" \
    "if [ ! -f $BRAIN ]; then cp $BRAIN_FS $BRAIN; fi"

    cmd "${H[$i]} Fix topology" \
    "mris_fix_topology -mgz -sphere qsphere.nofix -inflated inflated.nofix -orig orig.nofix -out orig.premesh -ga -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]}"

    cmd "${H[$i]} Fix orig.premesh header" \
    "mris_copy_header ${ORIG_PREMESH[$i]} ${ORIG_NOFIX[$i]} ${ORIG_PREMESH[$i]}"

    cmd "${H[$i]} Remesh (remesh 3 iters)" \
    "mris_remesh --remesh --iters 3 --input ${ORIG_PREMESH[$i]} --output ${ORIG[$i]}"

    cmd "${H[$i]} Remove self-intersections" \
    "mris_remove_intersection ${ORIG[$i]} ${ORIG[$i]}"

    cmd "${H[$i]} Smooth 2" \
    "mris_smooth -n 1 -nw -seed 1234 ${ORIG[$i]} ${SMOOTHW[$i]}"

    cmd "${H[$i]} Inflate 2 (with sulc)" \
    "mris_inflate ${ORIG[$i]} ${INFLATED[$i]}"

    cmd "${H[$i]} WM curv" \
    "mris_place_surface --curv-map ${ORIG[$i]} 2 10 ${CURV[$i]}"

    cmd "${H[$i]} WM area" \
    "mris_place_surface --area-map ${ORIG[$i]} ${AREA[$i]}"
done
fi

#################
## TAG 5: Edit brain.finalsurfs with GM from RIBBON_EDIT
#################
if ((TAG<=5))
then
for (( i=0; i<2; i++ ));
do
    if ((HEMI>=0 && HEMI!=i)); then continue; fi

    cmd "${H[$i]} Extract GM+WM from $RIBBON_WO_EDIT -> ${BMASK[$i]}" \
    "mri_extract_label $RIBBON_WO_EDIT ${LABEL_RIBBON_GM[$i]} ${LABEL_RIBBON_WM[$i]} ${BMASK[$i]}"

    cmd "${H[$i]} Replace 128->1 in ${BMASK[$i]}" \
    "mri_binarize --i ${BMASK[$i]} --o ${BMASK[$i]} --replace 128 1"

    cmd "${H[$i]} Mask $BRAIN_FINALSURFS_FS with ${BMASK[$i]} -> ${BRAIN_FINALSURFS_NO_CEREB[$i]}" \
    "mri_mask $BRAIN_FINALSURFS_FS ${BMASK[$i]} ${BRAIN_FINALSURFS_NO_CEREB[$i]}"

    cmd "${H[$i]} Extract WM from $RIBBON_WO_EDIT -> ${WM_BMASK[$i]}" \
    "mri_extract_label $RIBBON_WO_EDIT ${LABEL_RIBBON_WM[$i]} ${WM_BMASK[$i]}"

    cmd "${H[$i]} Concat WM_BMASK + brain.finalsurfs -> uniform WM 110" \
    "mri_concat --i ${WM_BMASK[$i]} --i ${BRAIN_FINALSURFS_NO_CEREB[$i]} --o ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]} --max"

    cmd "${H[$i]} Replace 128->110 in ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]}" \
    "mri_binarize --i ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]} --o ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]} --replace 128 110"

    cmd "${H[$i]} Extract GM from $RIBBON_WO_EDIT -> ${GM_BMASK[$i]}" \
    "mri_extract_label $RIBBON_WO_EDIT ${LABEL_RIBBON_GM[$i]} ${GM_BMASK[$i]}"

    cmd "${H[$i]} Concat GM_BMASK + uniform_wm_110 -> uniform GM 80" \
    "mri_concat --i ${GM_BMASK[$i]} --i ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]} --o ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --max"

    cmd "${H[$i]} Replace 128->80 in ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]}" \
    "mri_binarize --i ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --o ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --replace 128 80"
done
fi

#################
## TAG 6: Stats, cortex labels, sphere, aparc
#################
if ((TAG<=6))
then
for (( i=0; i<2; i++ ));
do
    if ((HEMI>=0 && HEMI!=i)); then continue; fi

    # Edit brain.finalsurfs: normalize GM intensity, zero HA region (labels 17,18,53,54 in aseg)
    cmd "${H[$i]} brain-finalsurfs-edit" \
    "python $O/brain-finalsurfs-edit.py ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]} ${GM_BMASK[$i]} ${BRAIN_FINALSURFS_NO_CEREB_EDITED[$i]} ${BRAIN_FINALSURFS_NO_CEREB_EDITED2[$i]} $ASEG_PADDED"

    cmd "${H[$i]} mris_autodet_gwstats" \
    "mris_autodet_gwstats --o ${AUTODET_NEW_GW_STATS[$i]} --i ${BRAIN_FINALSURFS_NO_CEREB_EDITED[$i]} --wm $WM --surf ${ORIG[$i]}"

    if ((CHANGE_AUTODET==1))
    then
        cmd "${H[$i]} Change pial_border_low to $PIAL_BORDER_LOW" \
        "sed -i'' -e 's/^pial_border_low[^/n]*/pial_border_low $PIAL_BORDER_LOW/' ${AUTODET_NEW_GW_STATS[$i]}"
    fi

    cmd "${H[$i]} cortex label (excl. subcortical)" \
    "mri_label2label --label-cortex ${ORIG[$i]} $ASEG_PRESURF_WO_SUBC 0 ${CORTEX_LABEL[$i]}"

    cmd "${H[$i]} cortex+hipamyg label" \
    "mri_label2label --label-cortex ${ORIG[$i]} $ASEG_PRESURF_WO_SUBC 1 ${CORTEX_HIPAMYG_LABEL[$i]}"

    cmd "${H[$i]} Sphere" \
    "mris_sphere -seed 1234 ${INFLATED[$i]} ${SPHERE[$i]}"

    cmd "${H[$i]} Surface registration" \
    "mris_register -curv ${SPHERE[$i]} ${FOLDING_ATLAS_ACFB40[$i]} ${SPHERE_REG[$i]}"

    cmd "${H[$i]} Cortical parcellation (DK atlas)" \
    "mris_ca_label -l ${CORTEX_LABEL[$i]} -aseg $ASEG_PRESURF -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]} ${SPHERE_REG[$i]} ${DKAPARC_ATLAS_ACFB40[$i]} ${APARC_ANNOT[$i]}"
done
fi

#################
## TAG 7: Pial surface
#################
if ((TAG<=7))
then
for (( i=0; i<2; i++ ));
do
    if ((HEMI>=0 && HEMI!=i)); then continue; fi

    cmd "${H[$i]} mris_place_surface --pial" \
    "mris_place_surface --i ${ORIG[$i]} --o ${RIBBON_EDIT_PIAL[$i]} --nsmooth 2 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${ORIG[$i]} --invol ${BRAIN_FINALSURFS_NO_CEREB_EDITED2[$i]} --threads 6 --white-surf ${ORIG[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg $ASEG_PRESURF_WO_SUBC --no-rip"
done
fi

#################
## TAG 8: Non-hemisphere files
#################
if ((TAG<=8))
then
cmd "Copy $RAWAVG_FS -> $RAWAVG" \
"cp $RAWAVG_FS $RAWAVG"
cmd "Copy $RAWAVG -> $RAWAVG_MASKED, then mask" \
"cp $RAWAVG $RAWAVG_MASKED"
cmd "Mask $RAWAVG_MASKED" \
"mri_mask $RAWAVG_MASKED $BRAIN_MASK $RAWAVG_MASKED"

cmd "Copy $ORIG_FS -> $ORIG_VOLUME" \
"cp $ORIG_FS $ORIG_VOLUME"
cmd "Copy $ORIG_FS -> $ORIG_MASKED, then mask" \
"cp $ORIG_FS $ORIG_MASKED"
cmd "Mask $ORIG_MASKED" \
"mri_mask $ORIG_MASKED $BRAIN_MASK $ORIG_MASKED"

cmd "Copy $IMAGE_ORIG_FS -> $IMAGE_ORIG" \
"cp $IMAGE_ORIG_FS $IMAGE_ORIG"

cmd "Copy $BRAIN_FINALSURFS_FS -> $BRAIN_FINALSURFS" \
"cp $BRAIN_FINALSURFS_FS $BRAIN_FINALSURFS"

cmd "Copy $BRAIN_FINALSURFS_MANEDIT_FS -> $BRAIN_FINALSURFS_MANEDIT (if exists)" \
"if [ -f $BRAIN_FINALSURFS_MANEDIT_FS ]; then cp $BRAIN_FINALSURFS_MANEDIT_FS $BRAIN_FINALSURFS_MANEDIT; fi"

cmd "AntsDenoiseImageFs $BRAIN -> $ANTSDN_BRAIN" \
"AntsDenoiseImageFs -i $BRAIN -o $ANTSDN_BRAIN"
cmd "mri_segment $ANTSDN_BRAIN -> $WM_SEG" \
"mri_segment -wsizemm 13 -mprage $ANTSDN_BRAIN $WM_SEG"

cmd "Copy $BRAINMASK_FS -> $BRAINMASK" \
"cp $BRAINMASK_FS $BRAINMASK"
fi

#################
## TAG 9: autorecon3 curvature and thickness stats
#################
if ((TAG<=9))
then
for (( i=0; i<2; i++ ));
do
    if ((HEMI>=0 && HEMI!=i)); then continue; fi

    cmd "${H[$i]} Copy orig -> white, pial -> pial" \
    "cp ${ORIG[$i]} ${WHITE[$i]}"
    cmd "${H[$i]} Copy pial" \
    "cp ${RIBBON_EDIT_PIAL[$i]} ${PIAL[$i]}"

    cmd "${H[$i]} Pial curv" \
    "mris_place_surface --curv-map ${PIAL[$i]} 2 10 ${PIAL_CURV[$i]}"
    cmd "${H[$i]} Pial area" \
    "mris_place_surface --area-map ${PIAL[$i]} ${PIAL_AREA[$i]}"
    cmd "${H[$i]} Thickness" \
    "mris_place_surface --thickness ${ORIG[$i]} ${PIAL[$i]} 20 5 ${THICKNESS[$i]}"

    cmd "${H[$i]} Curvature stats" \
    "mris_curvature_stats -m --writeCurvatureFiles -G -o ${CURV_STATS[$i]} -F smoothwm $SUBJID/$OUTPUT_FOLDER ${H[$i]} curv sulc"
done
fi

#################
## TAG 10: Cortical ribbon, aparc+aseg, wmparc, segstats
#################
if ((TAG<=10))
then
cmd "Cortical ribbon mask" \
"mris_volmask --aseg_name aseg.presurf --label_left_white ${LABEL_RIBBON_WM[0]} --label_left_ribbon ${LABEL_RIBBON_GM[0]} --label_right_white ${LABEL_RIBBON_WM[1]} --label_right_ribbon ${LABEL_RIBBON_GM[1]} --save_ribbon --out_root ribbon $SUBJID/$OUTPUT_FOLDER"

for (( i=0; i<2; i++ ));
do
    if ((HEMI>=0 && HEMI!=i)); then continue; fi

    cmd "${H[$i]} Cortical Parc 2 (CD atlas)" \
    "mris_ca_label -l ${CORTEX_LABEL[$i]} -aseg $ASEG_PRESURF -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]} ${SPHERE_REG[$i]} ${CD_APARC_ATLAS[$i]} ${CD_APARC_ANNOT[$i]}"

    cmd "${H[$i]} Cortical Parc 3 (DKT atlas)" \
    "mris_ca_label -l ${CORTEX_LABEL[$i]} -aseg $ASEG_PRESURF -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]} ${SPHERE_REG[$i]} ${DKT_APARC_ATLAS[$i]} ${DKT_APARC_ANNOT[$i]}"

    cmd "${H[$i]} Change SUBJECTS_DIR to SUBJECTS_DIR/SUBJID" \
    "export SUBJECTS_DIR=$SUBJECTS_DIR/$SUBJID"
    cmd "${H[$i]} WM/GM Contrast" \
    "pctsurfcon --s $OUTPUT_FOLDER --${H[$i]}-only"
    cmd "${H[$i]} Restore SUBJECTS_DIR" \
    "export SUBJECTS_DIR=$(dirname $SUBJECTS_DIR)"
done

cmd "Relabel hypointensities" \
"mri_relabel_hypointensities $ASEG_PRESURF $O/surf $ASEG_PRESURF_HYPOS"

cmd "APas-to-ASeg" \
"mri_surf2volseg --o $ASEG --i $ASEG_PRESURF_HYPOS --fix-presurf-with-ribbon $RIBBON_NEW --threads 1 --${H[0]}-cortex-mask ${CORTEX_LABEL[0]} --${H[0]}-white ${ORIG[0]} --${H[0]}-pial ${RIBBON_EDIT_PIAL[0]} --${H[1]}-cortex-mask ${CORTEX_LABEL[1]} --${H[1]}-white ${ORIG[1]} --${H[1]}-pial ${RIBBON_EDIT_PIAL[1]}"

cmd "AParc+ASeg" \
"mri_surf2volseg --o $APARC_PLUS_ASEG --label-cortex --i $ASEG --threads 1 --${H[0]}-annot ${APARC_ANNOT[0]} 1000 --${H[0]}-cortex-mask ${CORTEX_LABEL[0]} --${H[0]}-white ${ORIG[0]} --${H[0]}-pial ${RIBBON_EDIT_PIAL[0]} --${H[1]}-annot ${APARC_ANNOT[1]} 2000 --${H[1]}-cortex-mask ${CORTEX_LABEL[1]} --${H[1]}-white ${ORIG[1]} --${H[1]}-pial ${RIBBON_EDIT_PIAL[1]}"
fi

#################
## TAG 11: Extended parcellations (BA_exvivo, vpnl) + anatomical stats
#################
if ((TAG<=11))
then
for (( i=0; i<2; i++ ));
do
    if ((HEMI>=0 && HEMI!=i)); then continue; fi

    cmd "${H[$i]} aparc.a2009s" \
    "mri_surf2volseg --o $APARC_A2009S_ASEG --label-cortex --i $ASEG --threads 4 --${H[$i]}-annot ${APARC_A2009S_ANNOT[$i]} 11100 --${H[$i]}-cortex-mask ${CORTEX_LABEL[$i]} --${H[$i]}-white ${ORIG[$i]} --${H[$i]}-pial ${RIBBON_EDIT_PIAL[$i]} --${H[$((1-i))]}-annot ${APARC_A2009S_ANNOT[$((1-i))]} 12100 --${H[$((1-i))]}-cortex-mask ${CORTEX_LABEL[$((1-i))]} --${H[$((1-i))]}-white ${ORIG[$((1-i))]} --${H[$((1-i))]}-pial ${RIBBON_EDIT_PIAL[$((1-i))]}"

    cmd "${H[$i]} aparc.DKTatlas" \
    "mri_surf2volseg --o $APARC_DKT_ATLAS_ASEG --label-cortex --i $ASEG --threads 4 --${H[$i]}-annot ${DKT_APARC_ANNOT[$i]} 1000 --${H[$i]}-cortex-mask ${CORTEX_LABEL[$i]} --${H[$i]}-white ${ORIG[$i]} --${H[$i]}-pial ${RIBBON_EDIT_PIAL[$i]} --${H[$((1-i))]}-annot ${APARC_A2009S_ANNOT[$((1-i))]} 2000 --${H[$((1-i))]}-cortex-mask ${CORTEX_LABEL[$((1-i))]} --${H[$((1-i))]}-white ${ORIG[$((1-i))]} --${H[$((1-i))]}-pial ${RIBBON_EDIT_PIAL[$((1-i))]}"

    cmd "${H[$i]} WMParc" \
    "mri_surf2volseg --o $WMPARC --label-wm --i $APARC_PLUS_ASEG --threads 4 --${H[$i]}-annot ${APARC_ANNOT[$i]} 3000 --${H[$i]}-cortex-mask ${CORTEX_LABEL[$i]} --${H[$i]}-white ${ORIG[$i]} --${H[$i]}-pial ${RIBBON_EDIT_PIAL[$i]} --${H[$((1-i))]}-annot ${APARC_A2009S_ANNOT[$((1-i))]} 4000 --${H[$((1-i))]}-cortex-mask ${CORTEX_LABEL[$((1-i))]} --${H[$((1-i))]}-white ${ORIG[$((1-i))]} --${H[$((1-i))]}-pial ${RIBBON_EDIT_PIAL[$((1-i))]}"

    cmd "${H[$i]} Copy talairach.xfm" \
    "cp $TALAIRACH_XFM_FS $TALAIRACH_XFM"

    cmd "${H[$i]} WMParc stats" \
    "mri_segstats --seed 1234 --seg $WMPARC --sum $WMPARC_STATS --pv $NORM_FS --excludeid 0 --brainmask $BRAIN_FINALSURFS_FS --in $NORM_FS --in-intensity-name norm --in-intensity-units MR --subject $SUBJID/$OUTPUT_FOLDER --surf-wm-vol --ctab $WMPARC_STATS_LUT --etiv"

    cmd "${H[$i]} Change SUBJID" \
    "export SUBJID=$SUBJID/$OUTPUT_FOLDER"

    cmd "${H[$i]} Parcellation stats white" \
    "mris_anatomical_stats -th3 -mgz -noglobal -cortex ${CORTEX_LABEL[$i]} -f ${APARC_STATS[$i]} -b -a ${APARC_ANNOT[$i]} -c $APARC_ANNOT_CTAB $SUBJID ${H[$i]} white"
    cmd "${H[$i]} Parcellation stats pial" \
    "mris_anatomical_stats -th3 -mgz -noglobal -cortex ${CORTEX_LABEL[$i]} -f ${APARC_PIAL_STATS[$i]} -b -a ${APARC_ANNOT[$i]} -c $APARC_ANNOT_CTAB $SUBJID ${H[$i]} pial"
    cmd "${H[$i]} Parcellation stats a2009s" \
    "mris_anatomical_stats -th3 -mgz -noglobal -cortex ${CORTEX_LABEL[$i]} -f ${APARC_A2009S_STATS[$i]} -b -a ${APARC_A2009S_ANNOT[$i]} -c $APARC_A2009S_CTAB $SUBJID ${H[$i]} white"
    cmd "${H[$i]} Parcellation stats DKT" \
    "mris_anatomical_stats -th3 -mgz -noglobal -cortex ${CORTEX_LABEL[$i]} -f ${DKT_APARC_STATS[$i]} -b -a ${DKT_APARC_ANNOT[$i]} -c $DKT_APARC_CTAB $SUBJID ${H[$i]} white"

    cmd "${H[$i]} Restore SUBJID" \
    "export SUBJID=\`echo \"\$SUBJID\" | cut -d/ -f1-1\`"

    cmd "ASeg stats" \
    "mri_segstats --seed 1234 --seg $ASEG --sum $ASEG_STATS --pv $NORM_FS --empty --brainmask $BRAIN_FINALSURFS_FS --brain-vol-from-seg --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in $NORM_FS --in-intensity-name norm --in-intensity-units MR --etiv --euler --ctab $ASEG_STATS_LUT --subject $SUBJID/$OUTPUT_FOLDER --no-global-stats"

    cmd "Create fsaverage symlink if needed" \
    "if [ ! -d \"$SUBJECTS_DIR/fsaverage\" ]; then ln -s $FSAVERAGE $SUBJECTS_DIR; fi"

    cmd "${H[$i]} BA1_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA1_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA1_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA2_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA2_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA2_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA3a_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA3a_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA3A_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA3b_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA3b_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA3B_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA4a_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA4a_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA4A_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA4p_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA4p_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA4P_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA6_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA6_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA6_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA44_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA44_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA44_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA45_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA45_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA45_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} V1_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.V1_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${V1_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} V2_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.V2_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${V2_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} MT_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.MT_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${MT_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} entorhinal_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.entorhinal_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${ENTORHINAL_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} perirhinal_exvivo" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.perirhinal_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${PERIRHINAL_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} FG1 vpnl" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.FG1.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${FG1_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} FG2 vpnl" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.FG2.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${FG2_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} FG3 vpnl" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.FG3.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${FG3_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} FG4 vpnl" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.FG4.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${FG4_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} hOc1 vpnl" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.hOc1.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${HOC1_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} hOc2 vpnl" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.hOc2.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${HOC2_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} hOc3v vpnl" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.hOc3v.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${HOC3V_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} hOc4v vpnl" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.hOc4v.mpm.vpnl.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${HOC4V_MPM_VPNL_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"

    if [ -f "$O/label/${H[$i]}.mpm.vpnl.annot" ]; then
    var=$(date +%F_%H-%M-%S)
    cmd "Rename existing ${H[$i]}.mpm.vpnl.annot" \
    "mv $O/label/${H[$i]}.mpm.vpnl.annot $O/label/${H[$i]}.mpm.vpnl.annot_$var"
    fi
    cmd "${H[$i]} vpnl annot" \
    "mris_label2annot --s $SUBJID/$OUTPUT_FOLDER --ctab $COLORTABLE_VPNL_TXT --hemi ${H[$i]} --a mpm.vpnl --maxstatwinner --noverbose --l ${FG1_MPM_VPNL_LABEL[$i]} --l ${FG2_MPM_VPNL_LABEL[$i]} --l ${FG3_MPM_VPNL_LABEL[$i]} --l ${FG4_MPM_VPNL_LABEL[$i]} --l ${HOC1_MPM_VPNL_LABEL[$i]} --l ${HOC2_MPM_VPNL_LABEL[$i]} --l ${HOC3V_MPM_VPNL_LABEL[$i]} --l ${HOC4V_MPM_VPNL_LABEL[$i]}"

    cmd "${H[$i]} BA1_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA1_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA1_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA2_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA2_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA2_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA3a_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA3a_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA3A_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA3b_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA3b_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA3B_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA4a_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA4a_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA4A_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA4p_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA4p_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA4P_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA6_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA6_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA6_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA44_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA44_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA44_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} BA45_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.BA45_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${BA45_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} V1_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.V1_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${V1_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} V2_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.V2_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${V2_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} MT_exvivo thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.MT_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${MT_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} entorhinal thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.entorhinal_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${ENTORHINAL_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
    cmd "${H[$i]} perirhinal thresh" \
    "mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.perirhinal_exvivo.thresh.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${PERIRHINAL_EXVIVO_THRESH_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"

    cmd "Rename ${BA_EXVIVO_ANNOT[$i]} if exists" \
    "if [ -f ${BA_EXVIVO_ANNOT[$i]} ]; then mv ${BA_EXVIVO_ANNOT[$i]} ${BA_EXVIVO_ANNOT[$i]}_old_\$(date +%F_%H-%M-%S); fi"
    cmd "${H[$i]} BA_exvivo annot" \
    "mris_label2annot --s $SUBJID/$OUTPUT_FOLDER --hemi ${H[$i]} --ctab $COLORTABLE_BA_TXT --l ${BA1_EXVIVO_LABEL[$i]} --l ${BA2_EXVIVO_LABEL[$i]} --l ${BA3A_EXVIVO_LABEL[$i]} --l ${BA3B_EXVIVO_LABEL[$i]} --l ${BA4A_EXVIVO_LABEL[$i]} --l ${BA4P_EXVIVO_LABEL[$i]} --l ${BA6_EXVIVO_LABEL[$i]} --l ${BA44_EXVIVO_LABEL[$i]} --l ${BA45_EXVIVO_LABEL[$i]} --l ${V1_EXVIVO_LABEL[$i]} --l ${V2_EXVIVO_LABEL[$i]} --l ${MT_EXVIVO_LABEL[$i]} --l ${PERIRHINAL_EXVIVO_LABEL[$i]} --l ${ENTORHINAL_EXVIVO_LABEL[$i]} --a BA_exvivo --maxstatwinner --noverbose"

    cmd "Rename ${BA_EXVIVO_THRESH_ANNOT[$i]} if exists" \
    "if [ -f ${BA_EXVIVO_THRESH_ANNOT[$i]} ]; then mv ${BA_EXVIVO_THRESH_ANNOT[$i]} ${BA_EXVIVO_THRESH_ANNOT[$i]}_old_\$(date +%F_%H-%M-%S); fi"
    cmd "${H[$i]} BA_exvivo.thresh annot" \
    "mris_label2annot --s $SUBJID/$OUTPUT_FOLDER --hemi ${H[$i]} --ctab $COLORTABLE_BA_TXT --l ${BA1_EXVIVO_THRESH_LABEL[$i]} --l ${BA2_EXVIVO_THRESH_LABEL[$i]} --l ${BA3A_EXVIVO_THRESH_LABEL[$i]} --l ${BA3B_EXVIVO_THRESH_LABEL[$i]} --l ${BA4A_EXVIVO_THRESH_LABEL[$i]} --l ${BA4P_EXVIVO_THRESH_LABEL[$i]} --l ${BA6_EXVIVO_THRESH_LABEL[$i]} --l ${BA44_EXVIVO_THRESH_LABEL[$i]} --l ${BA45_EXVIVO_THRESH_LABEL[$i]} --l ${V1_EXVIVO_THRESH_LABEL[$i]} --l ${V2_EXVIVO_THRESH_LABEL[$i]} --l ${MT_EXVIVO_THRESH_LABEL[$i]} --l ${PERIRHINAL_EXVIVO_THRESH_LABEL[$i]} --l ${ENTORHINAL_EXVIVO_THRESH_LABEL[$i]} --a BA_exvivo.thresh --maxstatwinner --noverbose"

    cmd "${H[$i]} BA_exvivo stats" \
    "mris_anatomical_stats -th3 -mgz -noglobal -f ${BA_EXVIVO_STATS[$i]} -b -a ${BA_EXVIVO_ANNOT[$i]} -c $BA_EXVIVO_CTAB $SUBJID/$OUTPUT_FOLDER ${H[$i]} white"
    cmd "${H[$i]} BA_exvivo.thresh stats" \
    "mris_anatomical_stats -th3 -mgz -f ${BA_EXVIVO_THRESH_STATS[$i]} -noglobal -b -a ${BA_EXVIVO_THRESH_ANNOT[$i]} -c $BA_EXVIVO_THRESH_CTAB $SUBJID/$OUTPUT_FOLDER ${H[$i]} white"
done
fi

} ### END OF MAIN FUNCTION


if ((MULTICASE==0));
then
    main

elif ((MULTICASE==1));
then
declare -i counter=-1
for SUB in $SUBJECTS_DIR/*/;
do
    if ((FS==1));
    then
        IMAGE="$(find $SUB -maxdepth 1 -name "*T1*")"
    fi
    ASEG="$(find $SUB -maxdepth 1 -name "*aseg*" | head -1)"

    export var="${SUB: -1}"
    if [[ "$var" == "/" ]]; then export SUB="${SUB:0:-1}"; fi

    SUBJID="$(echo ${SUB##*/})"
    if ((SUBJID=="fsaverage")); then continue; fi

    counter+=1
    main &
    if ((counter%N_PARALLEL_COMPUTING==(N_PARALLEL_COMPUTING-1)));
    then
        wait
    fi
done
fi
