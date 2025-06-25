#!/bin/bash

#################
## Help
#################
Help ()
{
builtin echo "
AUTHOR: Benoît Verreman

LAST UPDATE: 2025-06-26

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
Modify the value of some constants in the script if needed: LABEL_RIBBON_WM, LABEL_RIBBON_GM, PIAL_BORDER_LOW

EXAMPLES:
$ bash ribbon_edit_script.sh -i 133019_T1w_acpc_dc_restore.nii.gz -s 133019 -b 133019_ribbon.nii.gz -c 133019_subcortical.nii.gz -n '5 15 29 30 32 31' -o 133019_ribbon.nii.gz -x 21 -y 22

$ bash ribbon_edit_script.sh -s <subjid> -t 8 -r 
#start from pial computation step (8), right hemisphere only (r)

$ bash ribbon_edit_script.sh -f <FOLDER_DATASET> -i * -r -n '5 15 29 30 32 31' -x 21 -y 22
#whole folder (f), apply recon-all (-i *), right hemisphere only (r)

PARAMETERS:

HELP
-h: Print this string, and exit

INPUT FILES
-i: Relative or absolute path to T1w image file
-s: Relative or absolute path to subjid folder (NECESSARY)
-b: Relative or absolute path to ribbon file (you may use labels 21 and 22 for HA in ribbon)

-c: Relative or absolute path to subcortical file
-n: List of labels (LABELS_SUBCORTICAL) to be used in subcortical file (-c)
#for aparc+aseg.nii.gz : LABELS_SUBCORTICAL='7 8 15 16 46 47' (default, no need to use -n option)
#for previous subcortical files : LABELS_SUBCORTICAL='5 15 29 30 32 31'

-o: Origin image of HA labels
-x: List of labels for HA left #'5 17 18' in aparc+aseg / '18 20 22' in HOA_Subcortical_Labels
-y: List of labels for HA right #'44 53 54' in aparc+aseg / '19 21 23' in HOA_Subcortical_Labels

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
-t 10: (aseg+aparc) Compute stats and other files
-t 12: BONUS: test some command lines

HA IMPROVEMENT
-d: Number of dilation in ha_ribbon_edit.py
-e: Number of erosion in ha_ribbon_edit.py

VALUES
-p: Give value of PIAL_BORDER_LOW

HEMI
-r: Compute only right hemisphere surface
-l: Compute only left hemisphere surface

RESET
-d: Reset outputs folder and report.sh script

WHOLE FOLDER
-f: Execute the script for each subfolder in given folder. 
Each subfolder name will be used as SUBJID and the given folder will be SUBJECTS_DIR.
Each subfolder should contain:
	- a ribbon file with a name containing 'ribbon' (the only one)
	- a subcortical file with a name containing 'subcortical' (the only one)
	- a T1 image (containing 'T1') if '-i *' is used as an option (execute recon-all)
	IF NOT: each subfolder should contain recon-all output (mri folder at least)

TROUBLESHOOTS:
-Missing argument: check if you put the necessary option flags, and for each flag, if it needs an argument or not.
"
}

#################
## Default global variables
#################
TAG=-1 # Start from beginning (option -t not used)
HEMI=-1 # Both hemispheres (option -r or -l not used)
FS=0 # Default: No recon-all (option -i not used)
MULTICASE=0 # Default: Only one case (option -f not used)
OUTPUT_FOLDER="outputs"
N_DILATION=3
N_EROSION=2
#LABELS_SUBCORTICAL="7 8 15 16 46 47" #for previous subcortical files (given to -c) : "5 15 29 30 32 31"
#declare -a LABELS_HA=("5 17 18" "44 53 54") #HA in aparc+aseg / ('18 20 22' '19 21 23') in HOA_Subcortical_Labels
declare -a LABELS_HA_RIB=("21" "22") #HA in ribbon (-n flag)
declare -a H=("lh" "rh") #Left then Right hemispheres
declare -a LABEL_RIBBON_WM=("2" "41")
declare -a LABEL_RIBBON_GM=("3" "42")
declare -i N_PARALLEL_COMPUTING=3 #Number of images computed in parallel at the same time

CHANGE_AUTODET=1 # Default: Change autodet with parameters bellow or given as option
PIAL_BORDER_LOW=5

#################
## Manage flags
#################
string_arguments="" #String of the different arguments in the command line

unset -v IMAGE
unset -v SUBJID
unset -v RIBBON
unset -v SUBCORTICAL
unset -v LABELS_SUBCORTICAL
unset -v HA
unset -v LABELS_HA_LEFT
unset -v LABELS_HA_RIGHT

#If a character is followed by :, then it needs an argument just after
VALID_ARGS="i:s:b:c:n:o:x:y:t:p:f:d:e:hlrk"

while getopts ${VALID_ARGS} opt; do
  case ${opt} in
    i) #T1 image path
        IMAGE=${OPTARG}
        FS=1
        string_arguments+="-i ${OPTARG} "
        ;;
    s) #Output folder name
        SUBJID=${OPTARG}
        string_arguments+="-s ${OPTARG} "
        ;;
    b) #ribbon image path
        RIBBON=${OPTARG}
        string_arguments+="-b ${OPTARG} "
        ;;
    c) #subcortical image path
        SUBCORTICAL=${OPTARG}
        string_arguments+="-c ${OPTARG} "
        ;;
    n) #subcortical labels
        LABELS_SUBCORTICAL=${OPTARG} #for previous subcortical files : "5 15 29 30 32 31"
        string_arguments+="-n ${OPTARG} "
        ;;
    o) #HA image path
        HA=${OPTARG}
        string_arguments+="-c ${OPTARG} "
        ;;
    x) #HA labels left
        LABELS_HA_LEFT=${OPTARG}
        string_arguments+="-c ${OPTARG} "
        ;;
    y) #HA labels right
        LABELS_HA_RIGHT=${OPTARG}
        string_arguments+="-c ${OPTARG} "
        ;;
    t) #Tag to (re)start script from
	TAG=${OPTARG}
	string_arguments+="-t ${OPTARG} "
	;;
    p) #Value to give to stat PIAL_BORDER_LOW
	PIAL_BORDER_LOW=${OPTARG}
	string_arguments+="-p ${OPTARG} "
	;;
    f) #Do script on all folder
	SUBJECTS_DIR=${OPTARG}
	MULTICASE=1
	string_arguments+="-f ${OPTARG} "
	;;
    d) #Number of dilation in ha_ribbon_edit.py
	N_DILATION=${OPTARG}
	string_arguments+="-f ${OPTARG} "
	;;
    e) #Number of erosion in ha_ribbon_edit.py
	N_EROSION=${OPTARG}
	string_arguments+="-f ${OPTARG} "
	;;
    h) #Help
	Help
	exit 1
	;;
    l) #Left hemisphere only
	HEMI=0 #left hemisphere only
	string_arguments+="-l "
	;;
    r) #Right hemisphere only
	HEMI=1 #right hemisphere only
	string_arguments+="-r "
	;;
    k) #Kill last OUTPUT_FODLER
	Delete #Kill report.sh and $OUTPUT_FOLDER
	string_arguments+="-d "
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
## Remove last character of $SUBJECTS_DIR if /
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

#################
## Remove slaches from $SUBJID
#################

# Test if user provided SUBJID (with -s, or used -f)
: ${SUBJID:?Missing argument -s}

#remove last / if any
export var="${SUBJID: -1}"
if [[ "$var" == "/" ]]; then
export SUBJID="${SUBJID:0:-1}"
fi

#remove first character if /
export var="${SUBJID:0:1}"
if [[ "$var" == "/" ]]; then
export SUBJID="${SUBJID:1}"
fi
	
O="$SUBJECTS_DIR/$SUBJID/$OUTPUT_FOLDER"

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
    builtin echo "$@" | tee -a $O/report.sh
}

#################
## Function to reset report.sh, $OUTPUT_FOLDER and mri_convert_correction_by_translation.py
#################
CreateFolders()
{
if [ ! -d "$SUBJECTS_DIR/$SUBJID" ]
then
	echo "Create $SUBJECTS_DIR/$SUBJID"
	mkdir $SUBJECTS_DIR/$SUBJID;
fi

if [ ! -d "$O" ]
then
	echo "Create $O and subfolders"
	mkdir $O;
	mkdir $O/scripts;
	mkdir $O/surf;
	mkdir $O/mri;
	mkdir $O/mri/orig;
	mkdir $O/mri/transforms;
	mkdir $O/label;
	mkdir $O/stats;
	mkdir $O/tmp;
	mkdir $O/touch;
	mkdir $O/trash;
fi
}

CreateScripts()
{
if [ ! -f "$O/report.sh" ]
then
touch $O/report.sh
Echo "#!/bin/bash"
fi

script_ha_ribbon_edit
script_nifti_padding
script_brain-finalsurfs-edit
script_edit_aseg_presurf_based_on_ribbon
script_expert_file
}

Delete()
{
cmd "Reset $O" \
"rm -r $O;"
}

#################
## Create a python script "ha_ribbon_edit.py"
#################
script_ha_ribbon_edit()
{
if [ ! -f "$O/ha_ribbon_edit.py" ]
then

cat > $O/ha_ribbon_edit.py <<EOF
import os
import nibabel as nib
import nibabel.processing #Used in nib.processing.conform
import scipy.ndimage #Used in nib.processing.conform
import sys #To add arguments
import copy #For deepcopy
import numpy #For motion1

#Arguments
path_ribbon = sys.argv[1]
path_aparc = sys.argv[2]
path_out = sys.argv[3]
labels_ha_lh = sys.argv[4]
labels_ha_rh = sys.argv[5]
labels_ha_rib_lh = sys.argv[6]
labels_ha_rib_rh = sys.argv[7]
n_dilation = sys.argv[8]
n_erosion = sys.argv[9]


#Load ribbon
if not os.path.isfile(path_ribbon):
    raise FileNotFoundError("The following path doesn't exist: " + path_ribbon)
else:
    img_ribbon = nib.load(path_ribbon)
data_ribbon = img_ribbon.get_fdata()
(a,b,c)=img_ribbon.header.get_data_shape()

#Load -o file (ie aparc+aseg)
if not os.path.isfile(path_aparc):
    raise FileNotFoundError("The following path doesn't exist: " + path_aparc)
else:
    img_aparc = nib.load(path_aparc)
data_aparc = img_aparc.get_fdata()

#Deepcopy of data_ribbon
data_out= copy.deepcopy(data_ribbon)

#labels into list
labels_ha_lh = list(map(int,labels_ha_lh.split()))
labels_ha_rh = list(map(int,labels_ha_rh.split()))
labels_ha_rib_lh = int(labels_ha_rib_lh)
labels_ha_rib_rh = int(labels_ha_rib_rh)

#Edit ribbon with HA from -o file
list_HA=[]
for x in range(a):
    for y in range(b):
        for z in range(c):
            aparc=data_aparc[x,y,z]
            if aparc in labels_ha_lh:
                list_HA.append((labels_ha_rib_lh,[x,y,z]))
                data_out[x,y,z]=labels_ha_rib_lh #21 #lh HA
            elif aparc in labels_ha_rh:
                list_HA.append((labels_ha_rib_rh,[x,y,z]))
                data_out[x,y,z]=labels_ha_rib_rh #22 #rh HA

###Add three diffusion and three erosion steps of HA in GM of ribbon
### Create labels class
class L:    
    l_wm = 2 # lh white matter
    r_wm = 41 # rh white matter
    
    l_gm = 3 # lh gray matter
    r_gm = 42 # rh gray matter
    
### Motions
d = numpy.eye(3, dtype=int).reshape(-1, 3)
motion1 = numpy.concatenate((d, -d), axis=0)

def dilation(data, list):
    data_out = copy.deepcopy(data)
    list_out = list[:]
    for (label,[x,y,z]) in list:
        n_coordinates = motion1 + [[x, y, z]]
        for (k,n,m) in n_coordinates:
            out_voxel = int(data[k,n,m])
            if out_voxel in [L.l_gm,L.r_gm]:
                data_out[k,n,m]=label
                list_out.append((label,[k,n,m]))
    return data_out,list_out

def erosion(data, list):
    data_out = copy.deepcopy(data)
    list_out = list[:]
    for (label,[x,y,z]) in list:
        n_coordinates = motion1 + [[x, y, z]]
        next_to_gm = False
        for (k,n,m) in n_coordinates:
            out_voxel = int(data[k,n,m])
            if out_voxel in [L.l_gm,L.r_gm]:
                next_to_gm = True
                data_out[x,y,z] = out_voxel
                break
        if not(next_to_gm):
            list_out.append((label,[x,y,z]))
    return data_out,list_out

data = copy.deepcopy(data_out)
list = list_HA[:]
for i in range(int(n_dilation)):
    data, list = dilation(data, list)
for i in range(int(n_erosion)):
    data, list = erosion(data, list)

###Create and save new image
img_out_rib = nib.Nifti1Image(data, img_ribbon.affine.copy())
nib.save(img_out_rib, path_out)
EOF
fi
}

#################
## Create a python script "nifti_padding.py"
#################
script_nifti_padding()
{
if [ ! -f "$O/nifti_padding.py" ]
then

cat > $O/nifti_padding.py <<EOF
import os
import nibabel as nib
import nibabel.processing #Used in nib.processing.conform
import scipy.ndimage #Used in nib.processing.conform
import sys #To add arguments
import subprocess

# SUBJID directory
img_in = sys.argv[1]
img_padded = sys.argv[2]
is_t1 = int(sys.argv[3])

#Load image to be treated
if not os.path.isfile(img_in):
    raise FileNotFoundError("The following path doesn't exist: " + img_in)
else:
    img = nib.load(img_in)

#Padding function: reshape the image to (max_dim, max_dim, max_dim) with same resolution and an orientation of 'LAS'
def padding(img, new_name):
    d = max(img.header.get_data_shape())
    new_img = nib.processing.conform(img, out_shape=(d, d, d),     voxel_size = img.header.get_zooms(), order=0, cval=0, orientation='LAS', out_class=None) #d
    nib.save(new_img, new_name)

if is_t1:
    padding(img, img_padded)
    print('apply_recon_all')
else:
    a = img.affine.copy()
    if int(a[0,3]) == 90:
        padding(img, img_padded)
        print('convert')
    else:
        print('no_convert')

EOF
fi
}

#################
## Create a python script "brain-finalsurfs-edit.py"
#################
script_brain-finalsurfs-edit()
{
if [ ! -f "$O/brain-finalsurfs-edit.py" ]
then

cat > $O/brain-finalsurfs-edit.py <<EOF
import os
import numpy as np #To compute motion
import nibabel as nib #To edit MRI images
import sys #To add arguments
import copy #For deepcopy

#Outside parameters
path_bf = sys.argv[1] #Brain.finalsurfs without the cerebellum
path_gmbm = sys.argv[2] #Gray Matter binary mask at 128 (by default) (based on ribbon-edit.mgz)
path_out = sys.argv[3] #Absolute path of the output brain.finalsurfs
path_out2 = sys.argv[4] #Absolute path of the output brain.finalsurfs
path_ha = sys.argv[5] #Put HA as background in brain.finalsurfs (in ribbon-edit)

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

#Load ribbon
if not os.path.isfile(path_ha):
    raise FileNotFoundError("Make sure the following path is correct: " + path_ha)
else:
    img_ha = nib.load(path_ha)
data_ha = img_ha.get_fdata()

#List of coordinates of the 26-nearest-neighbors around (0,0,0) plus itself
motion = np.transpose(np.indices((3,3,3)) - 1).reshape(-1, 3)

#Deepcopy of data_bf
data_bf_new = copy.deepcopy(data_bf)
data_bf_new2 = copy.deepcopy(data_bf)

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
                res = 80.0/(mean/nn)*val
                data_bf_new[x,y,z] = res
                data_bf_new2[x,y,z]= res
            if int(data_ha[x,y,z]) in [21,22]: #l_edit and r_edit
                data_bf_new2[x,y,z]=0

#Create and save new images
img_bf_new = nib.Nifti1Image(data_bf_new, img_bf.affine.copy())
nib.save(img_bf_new, path_out)

img_bf_new2 = nib.Nifti1Image(data_bf_new2, img_bf.affine.copy())
nib.save(img_bf_new2, path_out2)
EOF
fi
}

#################
## Create a python script "edit_aseg_presurf_based_on_ribbon.py"
#################
script_edit_aseg_presurf_based_on_ribbon()
{
if [ ! -f "$O/edit_aseg_presurf_based_on_ribbon.py" ]
then

cat > $O/edit_aseg_presurf_based_on_ribbon.py <<EOF
####
#Create aseg.presurf.mgz and aseg.presurf_wo_subc.mgz based on ribbon-edit.mgz
####

import os
import nibabel as nib
import nibabel.processing #Used in nib.processing.conform
import scipy.ndimage #Used in nib.processing.conform
import sys #To add arguments
import copy
import numpy

from itertools import product  #for motion2

#Outside parameters
path_aseg = sys.argv[1] #Path to aseg.presurf.old.mgz
path_ribbon = sys.argv[2] #Path to ribbon-edit.mgz
path_out = sys.argv[3] #Path to output aseg.presurf.mgz
path_out2 = sys.argv[4] #Path to output aseg.presurf_wo_subc.mgz

### Get images
if not os.path.isfile(path_aseg):
    raise FileNotFoundError("Make sure the following path is correct: " + path_aseg)
else:
    img_aseg = nib.load(path_aseg)

if not os.path.isfile(path_ribbon):
    raise FileNotFoundError("Make sure the following path is correct: " + path_ribbon)
else:
    img_ribbon = nib.load(path_ribbon)

### Get data
data_aseg = img_aseg.get_fdata() #Not to be edited
data_ribbon = img_ribbon.get_fdata() #Not to be edited

### Copy aseg
data_new_aseg = copy.deepcopy(data_aseg) #Copy to be edited
data_new_mask = copy.deepcopy(data_ribbon) #Copy to be edited

### Get dimensions
(a,b,c)=img_aseg.header.get_data_shape()

### Create labels class
class L:
    r_ce = 45 # right cerebellum exterior (in ribbon)
    cc_ant = 255 # CC_Anterior
    wmh = 77 # WM Hypointensities
    oc = 85 # optic chasm
    csf = 24 # Cerebrospinal fluid
    
    l_wm = 2 # lh white matter
    r_wm = 41 # rh white matter
    
    l_gm = 3 # lh gray matter
    r_gm = 42 # rh gray matter
    
    l_p = 12 # lh putamen
    r_p = 51 # rh putamen
    
    l_pd = 13 # lh pallidum
    r_pd = 52 # rh pallidum

    l_t = 10 # lh thalamus
    r_t = 49 # rh thalamus

    l_c = 11 # lh caudate
    r_c = 50 # rh caudate
    
    l_v = 28 # lh ventralDC
    r_v = 60 # rh ventralDC
    
    l_vs = 30 # lh left vessel
    r_vs = 62 # rh right vessel
    
    l_lv = 4 # lh lateral ventricle
    r_lv = 43 # rh lateral ventricle
    
    l_h = 17 # lh hippocampus
    r_h = 53 # rh hippocampus
    
    l_a = 18 # lh amygdala
    r_a = 54 # rh amygdala
    
    l_aa = 26 # lh Accumbens area
    r_aa = 58 # rh Accumbens area
    
    l_cp = 31 # lh choroid plexus
    r_cp = 63 # rh choroid plexus
    
    l_ilv = 5 # rh Inf lat vent
    r_ilv = 44 # rh Inf lat vent
    
    #CEREBELUM
    
    bs = 16 # brain stem
    
    fv = 15 # 4th ventricle
    
    l_cw = 7 # lh cerebellum white matter
    r_cw = 46 # rh cerebellum white matter
    
    l_cc = 8 # lh cerebellum cortex
    r_cc = 47 # rh cerebellum cortex
    
    #EDIT RIBBON
    
    l_edit = 21 #lh manuel edit in ribbon
    r_edit = 22 #rh manuel edit in ribbon

### Instantiate lists
list_LH=[] #Subcortical structures in lh
list_RH=[] #Subcortical structures in rh

### Modify BG, WM, GM, CC_anterior and subcortical structures based on ribbon
for x in range(a):
    for y in range(b):
        for z in range(c):
            ribbon_voxel = int(data_ribbon[x,y,z]) #ribbon-edit.mgz
            aseg_voxel = int(data_aseg[x,y,z])
            
            if aseg_voxel != ribbon_voxel: #Voxel for which label may have to be changed
                match ribbon_voxel:
            	    case 0 if aseg_voxel not in [L.l_cc,L.r_cc,L.l_cw,L.r_cw,L.bs,L.fv]: #correction of CEREBELLUM
            	        data_new_aseg[x,y,z] = ribbon_voxel
            	    case L.l_gm | L.r_gm: 
            	        data_new_aseg[x,y,z] = ribbon_voxel
            	    case L.l_edit:
            	        data_new_aseg[x,y,z] = L.l_a
            	    case L.r_edit:
            	        data_new_aseg[x,y,z] = L.r_a
            	    case L.l_wm:
            	        if aseg_voxel in [0, L.l_gm, L.bs]:
            	            data_new_aseg[x,y,z] = L.l_wm
            	        elif aseg_voxel in [L.l_p,L.l_pd,L.l_v,L.l_lv,L.l_aa,L.l_t,L.l_c,L.l_cp,L.l_ilv,L.l_h,L.l_a,L.l_vs,L.wmh,L.oc,L.csf]: #SUBCORTICAL structure lh only
            	            list_LH.append((aseg_voxel,[x,y,z]))
            	    case L.r_wm:
            	        if aseg_voxel in [0, L.r_gm, L.bs]:
            	            data_new_aseg[x,y,z] = L.r_wm
            	        elif aseg_voxel in [L.r_p,L.r_pd,L.r_v,L.r_lv,L.r_aa,L.r_t,L.r_c,L.r_cp,L.r_ilv,L.r_h,L.r_a,L.r_vs,L.wmh,L.oc,L.csf]: #SUBCORTICAL structure rh only
            	            list_RH.append((aseg_voxel,[x,y,z]))

### Copy data_new_aseg
data_new_aseg_wo_subc = copy.deepcopy(data_new_aseg) #Copy to be edited

### Replace different subcortical structures by WM (Putamen, VentralDC, CC_Anterior)
for (label,[x,y,z]) in list_LH:
    data_new_aseg_wo_subc[x,y,z] = L.l_wm
for (label,[x,y,z]) in list_RH:
    data_new_aseg_wo_subc[x,y,z] = L.r_wm
    
### Save both new aseg
img_new_aseg = nib.Nifti1Image(data_new_aseg, img_aseg.affine.copy())
nib.save(img_new_aseg, path_out)

img_new_aseg_wo_subc = nib.Nifti1Image(data_new_aseg_wo_subc, img_aseg.affine.copy())
nib.save(img_new_aseg_wo_subc, path_out2)
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
## Input files
#################
IMAGE_ORIG_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/orig/001.mgz"
RAWAVG_FS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/rawavg.mgz"

T1="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/T1.mgz"

NORM="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/norm.mgz"
ASEG_PRESURF_NOFIX="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/aseg.presurf.mgz"

BRAIN="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/brain.mgz"
BRAIN_FINALSURFS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/brain.finalsurfs.mgz"

#RB_ALL_WITHSKULL="$FREESURFER_HOME/average/RB_all_withskull_2016-05-10.vc700.gca"
#TALAIRACH_WITH_SKULL="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/transforms/talairach_with_skull.lta"

#RB_ALL="$FREESURFER_HOME/average/RB_all_2020-01-02.gca"
#CTRL_PTS="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/ctrl_pts.mgz"
#CC_UP="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/transforms/cc_up.lta"

WM="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/wm.mgz"

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

RAWAVG="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/rawavg.mgz"
ORIG_VOLUME="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/orig.mgz"

#ADDED for mri_segstats
TALAIRACH_XFM="$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer/mri/transforms/talairach.xfm"


#################
## Output files
#################
IMAGE_PADDED="$SUBJECTS_DIR/$SUBJID/image-padded.mgz"
$IMAGE_ORIG="$O/mri/orig/001.mgz"
$RAWAVG="$O/mri/rawavg.mgz"

RIBBON_PADDED="$O/mri/ribbon-precorrection.mgz"
SUBCORTICAL_PADDED="$O/mri/subcortical-precorrection.mgz"
RIBBON_CONVERT="$O/mri/ribbon-convert.mgz"
RIBBON_EDIT="$O/mri/ribbon-edit.mgz"
SUBCORTICAL_EDIT="$O/mri/subcortical-edit.mgz"

HA_PADDED="$O/mri/ha-padded.mgz"
HA_CONVERT="$O/mri/ha-convert.mgz"

SUBCORTICAL_MASK="$O/mri/subcortical-mask.mgz"
BRAIN_MASK="$O/mri/brain-mask.mgz"

T1_MASKED="$O/mri/T1-masked.mgz"

ASEG_PRESURF="$O/mri/aseg.presurf.mgz"
ASEG_PRESURF_WO_SUBC="$O/mri/aseg.presurf_wo_subc.mgz"

WM_BMASK_ALL="$O/mri/wm-bmask.mgz"
WM_MASK="$O/mri/wm-mask.mgz"
WM_CONCAT="$O/mri/wm-concat.mgz"
WM_BMASK_250="$O/mri/wm-bmask-250.mgz"
WM_ASEGEDIT="$O/mri/wm-asegedit.mgz"
WM_EDITED="$O/mri/wm.mgz" # Use this name for mri_fix_topology

BRAIN_COPY="$O/mri/brain.mgz" # for mris_fix_topology

declare -a FILLED_PRETRESS=("$O/mri/filled_pretress_lh.mgz" "$O/mri/filled_pretress_rh.mgz")
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

declare -a GM_BMASK=("$O/mri/gm-bmask_lh.mgz" "$O/mri/gm-bmask_rh.mgz")
declare -a WM_BMASK=("$O/mri/wm-bmask_lh.mgz" "$O/mri/wm-bmask_rh.mgz")
declare -a BMASK=("$O/mri/bmask_lh.mgz" "$O/mri/bmask_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB=("$O/mri/brain.finalsurfs_no_cereb_lh.mgz" "$O/mri/brain.finalsurfs_no_cereb_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110=("$O/mri/brain.finalsurfs_no_cereb_uniform_wm_110_lh.mgz" "$O/mri/brain.finalsurfs_no_cereb_uniform_wm_110_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80=("$O/mri/brain.finalsurfs_no_cereb_uniform_gm_80_lh.mgz" "$O/mri/brain.finalsurfs_no_cereb_uniform_gm_80_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB_EDITED=("$O/mri/brain.finalsurfs_no_cereb_edited_lh.mgz" "$O/mri/brain.finalsurfs_no_cereb_edited_rh.mgz")
declare -a BRAIN_FINALSURFS_NO_CEREB_EDITED2=("$O/mri/brain.finalsurfs_no_cereb_edited2_lh.mgz" "$O/mri/brain.finalsurfs_no_cereb_edited2_rh.mgz")

declare -a AUTODET_NEW_GW_STATS=("$O/surf/autodet-new.gw.stats.lh.dat" "$O/surf/autodet-new.gw.stats.rh.dat")

declare -a RIBBON_EDIT_PIAL=("$O/surf/lh.ribbon_edit.pial" "$O/surf/rh.ribbon_edit.pial")
declare -a RIBBON_EDIT_PIAL_SECOND_PASS=("$O/surf/lh.ribbon_edit-second-pass.pial" "$O/surf/rh.ribbon_edit-second-pass.pial")

declare -a RIBBON_EDIT_PIAL_THIRD_PASS=("$O/surf/lh.ribbon_edit.smooth-third-pass.pial" "$O/surf/rh.ribbon_edit.smooth-third-pass.pial")

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

TALAIRACH_XFM_COPY="$O/mri/transforms/talairach.xfm"

#################
## New invocation in report.sh and create
#################
#Test if $OUTPUT_FOLDER folder already exist, and if it does not, create one
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
## FreeSurfer 7.4.1 on $IMAGE creating $SUBJECTS_DIR/$SUBJID folder
#################
if ((FS == 1 && TAG < 0))
then
# Test if user provided $IMAGE
: ${IMAGE:?Missing argument -i}
Echo "# Given image: $IMAGE"

if [ -d "$SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer" ]
then
	Echo "Do not re-run FreeSurfer on same SUBJID: $SUBJECTS_DIR/$SUBJID/${SUBJID}_freesurfer"
	exit 1
fi

cmd "Compensate for future translation in recon-all" \
"result=`python $O/nifti_padding.py $IMAGE $IMAGE_PADDED 1`"

cmd "Add SUBJID to SUBJECTS_DIR" \
"export SUBJECTS_DIR=$SUBJECTS_DIR/$SUBJID"

#-xopts-overwrite is used when expert file already used before
cmd "Apply recon-all -autorecon 1 and 2 on $IMAGE_PADDED" \
"recon-all -autorecon1 -autorecon2 -s ${SUBJID}_freesurfer -i $IMAGE_PADDED -hires -parallel -openmp 4 -expert expert_file.txt -xopts-overwrite -cw256" 

cmd "Change back SUBJECTS_DIR/SUBJID to SUBJECTS_DIR" \
"export SUBJECTS_DIR=$(dirname $SUBJECTS_DIR)"

#Copy image-padded.mgz to orig/001.mgz
cmd "Copy $IMAGE_ORIG_FS" \
"if [ ! -f $IMAGE_ORIG ]; then cp $IMAGE_ORIG_FS $IMAGE_ORIG; fi" 

#Copy image-padded.mgz to orig/001.mgz
cmd "Copy $RAWAVG_FS to get $RAWAVG if not already exists" \
"if [ ! -f $RAWAVG ]; then cp $RAWAVG_FS $RAWAVG; fi" 

fi

#################
## Convert ribbon and subcortical
#################
if ((TAG<=0))
then
# Test if user provided RIBBON and SUBCORTICAL
: ${RIBBON:?Missing argument -b} ${SUBCORTICAL:?Missing argument -c} ${LABELS_SUBCORTICAL:?Missing argument -n}
Echo "# Given ribbon: $RIBBON"
Echo "# Given subcortical: $SUBCORTICAL"

#If needed, correcting the dimensions of the image with padding
cmd "Use script $O/nifti_padding.py on $RIBBON" \
"result=`python $O/nifti_padding.py $RIBBON $RIBBON_PADDED 0`"

if [[ "$result" == "convert" ]]; then
cmd "Convert $RIBBON_PADDED" \
"mri_convert $RIBBON_PADDED $RIBBON_CONVERT -rt nearest -ns 1 --conform_min"
else
cmd "Copy $RIBBON in $RIBBON_CONVERT" \
"cp $RIBBON $RIBBON_CONVERT" 
fi

cmd "Use script $O/nifti_padding.py on $SUBCORTICAL" \
"result=`python $O/nifti_padding.py $SUBCORTICAL $SUBCORTICAL_PADDED 0`"

if [[ "$result" == "convert" ]]; then
cmd "Convert $SUBCORTICAL_PADDED" \
"mri_convert $SUBCORTICAL_PADDED $SUBCORTICAL_EDIT -rt nearest -ns 1 --conform_min"
else
cmd "Copy $SUBCORTICAL in $SUBCORTICAL_EDIT" \
"cp $SUBCORTICAL $SUBCORTICAL_EDIT" 
fi

fi

#################
## Exctract labels from ribbon and subcortical into brain_mask
#################
if ((TAG<=1))
then

# Test if user provided HA information
: ${HA:?Missing argument -o} ${LABELS_HA_LEFT:?Missing argument -x} ${LABELS_HA_RIGHT:?Missing argument -y}

#If needed, correcting the dimensions of the image with padding
cmd "Use script $O/nifti_padding.py on $HA" \
"result=`python $O/nifti_padding.py $HA $HA_PADDED 0`"

if [[ "$result" == "convert" ]]; then
cmd "Convert $HA" \
"mri_convert $HA_PADDED $HA_CONVERT -rt nearest -ns 1 --conform_min"
else
cmd "Copy $HA in $HA_CONVERT" \
"cp $HA $HA_CONVERT" 
fi

cmd "Use script $O/ha_ribbon_edit.py on $RIBBON_CONVERT to add HA from $SUBCORTICAL_EDIT" \
"python $O/ha_ribbon_edit.py $RIBBON_CONVERT $HA_CONVERT $RIBBON_EDIT '${LABELS_HA_LEFT}' '${LABELS_HA_RIGHT}' ${LABELS_HA_RIB[0]} ${LABELS_HA_RIB[1]} $N_DILATION $N_EROSION"

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

#################
## Edit original aseg.presurf based on ribbon-edit
#################
if ((TAG<=4))
then

## Create aseg.presurf.mgz and aseg.presurf_wo_subc.mgz based on ribbon-edit.mgz
cmd "Use script $O/edit_aseg_presurf_based_on_ribbon.py on $ASEG_PRESURF_NOFIX" \
"python $O/edit_aseg_presurf_based_on_ribbon.py $ASEG_PRESURF_NOFIX $RIBBON_EDIT $ASEG_PRESURF $ASEG_PRESURF_WO_SUBC"

fi

#################
## Compute WM_EDIT based on BRAIN_FINALSURFS masked by WM_BMASK_ALL
#################
if ((TAG<=5))
then
# Extract white matter from ribbon-edit to create wm-bmask.mgz
cmd "Extract WM from $RIBBON_EDIT" \
"mri_extract_label $RIBBON_EDIT ${LABEL_RIBBON_WM[0]} ${LABEL_RIBBON_WM[1]} $WM_BMASK_ALL" #0/128 binary mask

cmd "Concatenate $WM_BMASK_ALL with $WM into $WM_CONCAT" \
"mri_concat --i $WM_BMASK_ALL --i $WM --o $WM_CONCAT --sum" #ROI at 378 (128+250)

cmd "Binarize $WM_CONCAT at 251 into $WM_BMASK_250" \
"mri_binarize --i $WM_CONCAT --o $WM_BMASK_250 --match 378"

cmd "Replace 1 by 250 into $WM_BMASK_250" \
"mri_binarize --i $WM_BMASK_250 --o $WM_BMASK_250 --replace 1 250"

# May also use $BRAIN_FINALSURFS
cmd "Mask $BRAIN with $WM_BMASK_ALL into $WM_MASK" \
"mri_mask -T 5 $BRAIN $WM_BMASK_ALL $WM_MASK"

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
        
        # Copy BRAIN for mris_fix_topology
	cmd "${H[$i]} Copy $BRAIN for mris_fix_topology" \
	"if [ ! -f $BRAIN_COPY ]; then cp $BRAIN $BRAIN_COPY; fi" 
	
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
	cmd "${H[$i]} Extract forebrain from $RIBBON_EDIT" \
	"mri_extract_label $RIBBON_EDIT ${LABEL_RIBBON_GM[$i]} ${LABEL_RIBBON_WM[$i]} ${BMASK[$i]}" #0/128 binary mask

	cmd "${H[$i]} Replace 128 by 1 into ${BMASK[$i]}" \
	"mri_binarize --i ${BMASK[$i]} --o ${BMASK[$i]} --replace 128 1"

	cmd "${H[$i]} Mask $BRAIN_FINALSURFS with ${BMASK[$i]} into ${BRAIN_FINALSURFS_NO_CEREB[$i]}" \
	"mri_mask $BRAIN_FINALSURFS ${BMASK[$i]} ${BRAIN_FINALSURFS_NO_CEREB[$i]}"

	# Extract white matter from ribbon-edit to create wm-bmask.mgz
	cmd "${H[$i]} Extract GM from $RIBBON_EDIT" \
	"mri_extract_label $RIBBON_EDIT ${LABEL_RIBBON_WM[$i]} ${WM_BMASK[$i]}" #0/128 binary mask

	cmd "${H[$i]} Concatenate ${WM_BMASK[$i]} with ${BRAIN_FINALSURFS_NO_CEREB[$i]} into ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]}" \
	"mri_concat --i ${WM_BMASK[$i]} --i ${BRAIN_FINALSURFS_NO_CEREB[$i]} --o ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]} --max"

	cmd "${H[$i]} Replace 128 by 110 in ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]}" \
	"mri_binarize --i ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]} --o ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]} --replace 128 110"
	
	# Extract gray matter from ribbon-edit to create gm-bmask.mgz, and create with it bf_80
	cmd "${H[$i]} Extract GM from $RIBBON_EDIT" \
	"mri_extract_label $RIBBON_EDIT ${LABEL_RIBBON_GM[$i]} ${GM_BMASK[$i]}" #0/128 binary mask

	cmd "${H[$i]} Concatenate ${GM_BMASK[$i]} with ${BRAIN_FINALSURFS_NO_CEREB[$i]} into ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]}" \
	"mri_concat --i ${GM_BMASK[$i]} --i ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]} --o ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --max"

	cmd "${H[$i]} Replace 128 by 80 in ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]}" \
	"mri_binarize --i ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --o ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --replace 128 80"
	fi
done
fi

#################
## Compute stats, labels and aparc, to prepare for surface computation 
#################
if ((TAG<=8))
then
for (( i=0; i<2; i++ ));
do
	if ((HEMI>=0 && HEMI!=i)); #Cases when the current hemi ($i) is not to be processed
	then
		continue;
	else	
	
	# Use script brain-finalsurfs-edit.py to edit brain.finalsurfs.mgz
	cmd "${H[$i]} Use script $O/brain-finalsurfs-edit.py on ${BRAIN_FINALSURFS_NO_CEREB[$i]} with ${GM_BMASK[$i]}" \
	"python $O/brain-finalsurfs-edit.py ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_WM_110[$i]} ${GM_BMASK[$i]} ${BRAIN_FINALSURFS_NO_CEREB_EDITED[$i]} ${BRAIN_FINALSURFS_NO_CEREB_EDITED2[$i]} $RIBBON_EDIT"
	
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
	"mri_label2label --label-cortex ${ORIG[$i]} $ASEG_PRESURF_WO_SUBC 0 ${CORTEX_LABEL[$i]}"
	
	# Compute labels to remove HIPOCAMPUS AND AMYGDALA from pial surface in mris_place_surface
	cmd "${H[$i]} Label2label for cortex" \
	"mri_label2label --label-cortex ${ORIG[$i]} $ASEG_PRESURF_WO_SUBC 1 ${CORTEX_HIPAMYG_LABEL[$i]}"
	
	# Compute APARC cortical parcellation for --aparc option in mris_place_surface
	cmd "${H[$i]} Sphere" \
 	"mris_sphere -seed 1234 ${INFLATED[$i]} ${SPHERE[$i]}"
 	cmd "${H[$i]} Surf Reg" \
 	"mris_register -curv ${SPHERE[$i]} ${FOLDING_ATLAS_ACFB40[$i]} ${SPHERE_REG[$i]}"
 	cmd "${H[$i]} Cortical Parc" \
 	"mris_ca_label -l ${CORTEX_LABEL[$i]} -aseg $ASEG_PRESURF -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]} ${SPHERE_REG[$i]} ${DKAPARC_ATLAS_ACFB40[$i]} ${APARC_ANNOT[$i]}"
 	
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
	# PIAL first pass: BRAIN_FINALSURFS_NO_CEREB_EDITED2
	cmd "${H[$i]} Computes pial surface - first pass" \
	"mris_place_surface --i ${ORIG[$i]} --o ${RIBBON_EDIT_PIAL[$i]} --nsmooth 1 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${ORIG[$i]} --invol ${BRAIN_FINALSURFS_NO_CEREB_EDITED2[$i]} --threads 6 --white-surf ${ORIG[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg $ASEG_PRESURF_WO_SUBC --no-rip"
	
	# PIAL second pass: BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80
	cmd "${H[$i]} Computes pial surface - second pass: i+w" \
	"mris_place_surface --i ${RIBBON_EDIT_PIAL[$i]} --o ${RIBBON_EDIT_PIAL_SECOND_PASS[$i]} --nsmooth 0 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${ORIG[$i]} --invol ${BRAIN_FINALSURFS_NO_CEREB_UNIFORM_GM_80[$i]} --threads 6 --white-surf ${RIBBON_EDIT_PIAL[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg $ASEG_PRESURF_WO_SUBC --no-rip"
	
	# PIAL third pass: nsmooth 1
	cmd "${H[$i]} Computes pial surface - third pass: smooth" \
	"mris_place_surface --i ${RIBBON_EDIT_PIAL_SECOND_PASS[$i]} --o ${RIBBON_EDIT_PIAL_THIRD_PASS[$i]} --nsmooth 1 --adgws-in ${AUTODET_NEW_GW_STATS[$i]} --pial --${H[$i]} --repulse-surf ${RIBBON_EDIT_PIAL_SECOND_PASS_i_w[$i]} --invol ${BRAIN_FINALSURFS_NO_CEREB_EDITED2[$i]} --threads 6 --white-surf ${ORIG[$i]} --pin-medial-wall ${CORTEX_LABEL[$i]} --seg $ASEG_PRESURF_WO_SUBC --no-rip"
	fi
done
fi

#################
## Add last steps of autorecon3: stats, aseg, labels
# #################
if ((TAG<=10))
then
for (( i=0; i<2; i++ ));
do
	if ((HEMI>=0 && HEMI!=i)); #Cases when the current hemi ($i) is not to be processed
	then
		continue;
	else
	# Copies to have the right names for the next functions
	cmd "${H[$i]} Copy white surface to ${WHITE[$i]}" \
	"cp ${ORIG[$i]} ${WHITE[$i]}" 
	cmd "${H[$i]} Copy pial surface to ${PIAL[$i]}" \
	"cp ${RIBBON_EDIT_PIAL_THIRD_PASS[$i]} ${PIAL[$i]}" 
	
	# Compute the stats
	cmd "${H[$i]} pial curv" \
 	"mris_place_surface --curv-map ${PIAL[$i]} 2 10 ${PIAL_CURV[$i]}"
 	cmd "${H[$i]} pial area" \
 	"mris_place_surface --area-map ${PIAL[$i]} ${PIAL_AREA[$i]}"
 	cmd "${H[$i]} thickness" \
 	"mris_place_surface --thickness ${ORIG[$i]} ${PIAL[$i]} 20 5 ${THICKNESS[$i]}"
 	# same command again for "area and vertex vol rh" ??? 	
 	
 	cmd "${H[$i]} Curvature Stats" \
 	"mris_curvature_stats -m --writeCurvatureFiles -G -o ${CURV_STATS[$i]} -F smoothwm $SUBJID/$OUTPUT_FOLDER ${H[$i]} curv sulc"

 	cmd "${H[$i]} Cortical ribbon mask" \
 	"mris_volmask --aseg_name aseg.presurf --label_left_white ${LABEL_RIBBON_WM[0]} --label_left_ribbon ${LABEL_RIBBON_GM[0]} --label_right_white ${LABEL_RIBBON_WM[1]} --label_right_ribbon ${LABEL_RIBBON_GM[1]} --save_ribbon --out_root ribbon_script_${H[$i]} --${H[$i]}-only $SUBJID/$OUTPUT_FOLDER" #Searching for surf/rh.white and surf/rh.pial, and --surf_white and --surf_pial don't help, can use arg 
	#cmd "Copy $RH_SMOOTHW_NOFIX to $RH_SMOOTHW" \
	#"cp $RH_SMOOTHW_NOFIX $RH_SMOOTHW"
  	#cmd "Inflation2 rh" \
 	#"mris_inflate -n 30 $RH_SMOOTHW $RH_INFLATED"

 	cmd "${H[$i]} Cortical Parc 2" \
 	"mris_ca_label -l ${CORTEX_LABEL[$i]} -aseg $ASEG_PRESURF -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]} ${SPHERE_REG[$i]} ${CD_APARC_ATLAS[$i]} ${CD_APARC_ANNOT[$i]}"
	#Use surf/rh.smoothwm and surf/rh.sphere.reg
 	cmd "${H[$i]} Cortical Parc 3" \
 	"mris_ca_label -l ${CORTEX_LABEL[$i]} -aseg $ASEG_PRESURF -seed 1234 $SUBJID/$OUTPUT_FOLDER ${H[$i]} ${SPHERE_REG[$i]} ${DKT_APARC_ATLAS[$i]} ${DKT_APARC_ANNOT[$i]}"
 	
 	cmd "${H[$i]} Copy $RAWAVG to $RAWAVG_MASKED" \
 	"cp $RAWAVG $RAWAVG_MASKED"
 	cmd "${H[$i]} Mask $RAWAVG_MASKED with $BRAIN_MASK into $RAWAVG_MASKED" \
"mri_mask $RAWAVG_MASKED $BRAIN_MASK $RAWAVG_MASKED"
 	cmd "${H[$i]} Copy $ORIG_VOLUME to $ORIG_MASKED" \
 	"cp $ORIG_VOLUME $ORIG_MASKED"
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
 	"if [ ! -d "$SUBJECTS_DIR/fsaverage" ]; then ln -s $FSAVERAGE $SUBJECTS_DIR; fi"
 	
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
 	"mri_label2label --srcsubject fsaverage --srclabel $SUBJECTS_DIR/fsaverage/label/${H[$i]}.perirhinal_exvivo.label --trgsubject $SUBJID/$OUTPUT_FOLDER --trglabel ${PERIRHINAL_EXVIVO_LABEL[$i]} --hemi ${H[$i]} --regmethod surface"
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
 	
 	if [ -f "$O/label/${H[$i]}.mpm.vpnl.annot" ]
	then
	var=$(date +%F_%H-%M-%S)
	cmd "Change name of $O/label/${H[$i]}.mpm.vpnl.annot to recompute it: ${H[$i]}.mpm.vpnl.annot_$var" \
	"mv $O/label/${H[$i]}.mpm.vpnl.annot $O/label/${H[$i]}.mpm.vpnl.annot_$var"
	fi
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
 	
	cmd "Change name of ${BA_EXVIVO_ANNOT[$i]} if already exists" \
 	"if [ -f ${BA_EXVIVO_ANNOT[$i]} ]; then mv ${BA_EXVIVO_ANNOT[$i]} ${BA_EXVIVO_ANNOT[$i]}_old_$(date +%F_%H-%M-%S); fi"
 	cmd "${H[$i]} mri_label2label ctab" \
 	"mris_label2annot --s $SUBJID/$OUTPUT_FOLDER --hemi ${H[$i]} --ctab $COLORTABLE_BA_TXT --l ${BA1_EXVIVO_LABEL[$i]} --l ${BA2_EXVIVO_LABEL[$i]} --l ${BA3A_EXVIVO_LABEL[$i]} --l ${BA3B_EXVIVO_LABEL[$i]} --l ${BA4A_EXVIVO_LABEL[$i]} --l ${BA4P_EXVIVO_LABEL[$i]} --l ${BA6_EXVIVO_LABEL[$i]} --l ${BA44_EXVIVO_LABEL[$i]} --l ${BA45_EXVIVO_LABEL[$i]} --l ${V1_EXVIVO_LABEL[$i]} --l ${V2_EXVIVO_LABEL[$i]} --l ${MT_EXVIVO_LABEL[$i]} --l ${PERIRHINAL_EXVIVO_LABEL[$i]} --l ${ENTORHINAL_EXVIVO_LABEL[$i]} --a BA_exvivo --maxstatwinner --noverbose"
 	
 	cmd "Change name of ${BA_EXVIVO_THRESH_ANNOT[$i]} if already exists" \
 	"if [ -f ${BA_EXVIVO_THRESH_ANNOT[$i]} ]; then mv ${BA_EXVIVO_THRESH_ANNOT[$i]} ${BA_EXVIVO_THRESH_ANNOT[$i]}_old_$(date +%F_%H-%M-%S); fi"
 	cmd "${H[$i]} mri_label2label ctab thresh" \
 	"mris_label2annot --s $SUBJID/$OUTPUT_FOLDER --hemi ${H[$i]} --ctab $COLORTABLE_BA_TXT --l ${BA1_EXVIVO_THRESH_LABEL[$i]} --l ${BA2_EXVIVO_THRESH_LABEL[$i]} --l ${BA3A_EXVIVO_THRESH_LABEL[$i]} --l ${BA3B_EXVIVO_THRESH_LABEL[$i]} --l ${BA4A_EXVIVO_THRESH_LABEL[$i]} --l ${BA4P_EXVIVO_THRESH_LABEL[$i]} --l ${BA6_EXVIVO_THRESH_LABEL[$i]} --l ${BA44_EXVIVO_THRESH_LABEL[$i]} --l ${BA45_EXVIVO_THRESH_LABEL[$i]} --l ${V1_EXVIVO_THRESH_LABEL[$i]} --l ${V2_EXVIVO_THRESH_LABEL[$i]} --l ${MT_EXVIVO_THRESH_LABEL[$i]} --l ${PERIRHINAL_EXVIVO_THRESH_LABEL[$i]} --l ${ENTORHINAL_EXVIVO_THRESH_LABEL[$i]} --a BA_exvivo.thresh --maxstatwinner --noverbose"

 	cmd "${H[$i]} mris_anatomical_stats ctab" \
 	"mris_anatomical_stats -th3 -mgz -noglobal -f ${BA_EXVIVO_STATS[$i]} -b -a ${BA_EXVIVO_ANNOT[$i]} -c $BA_EXVIVO_CTAB $SUBJID/$OUTPUT_FOLDER ${H[$i]} white"
 	
 	cmd "${H[$i]}mris_anatomical_stats ctab thresh" \
 	"mris_anatomical_stats -th3 -mgz -f ${BA_EXVIVO_THRESH_STATS[$i]} -noglobal -b -a ${BA_EXVIVO_THRESH_ANNOT[$i]} -c $BA_EXVIVO_THRESH_CTAB $SUBJID/$OUTPUT_FOLDER ${H[$i]} white"
	fi
done
fi

#################
## BONUS FOR TESTING SOMETHING ALONE: use tag -t 12
# #################
if ((TAG<=12))
then
for (( i=0; i<2; i++ ));
do
	if ((HEMI>=0 && HEMI!=i)); #Cases when the current hemi ($i) is not to be processed
	then
		continue;
	else
	echo "..."
	#echo "counter=${counter}" #Put your command lines below
	#sleep 5
	fi
done
fi

} ### END OF MAIN FUNCTION




if ((MULTICASE==0)); # ONLY ONE IMAGE; don't forget to export SUBJECTS_DIR
then
	main
	
elif ((MULTICASE==1)); # SEVERAL IMAGES IN SUBJECTS_DIR GIVEN WITH -f
then
declare -i counter=-1
for SUB in $SUBJECTS_DIR/*/;
do
	if ((FS==1)); # -i was used, need T1 for each subfolder
	then
		IMAGE="$(find $SUB -maxdepth 1 -name "*T1*")"
	fi
        RIBBON="$(find $SUB -maxdepth 1 -name "*ribbon*")"
        SUBCORTICAL="$(find $SUB -maxdepth 1 -name "*subcortical*")"
        
        HA=$RIBBON #TO BE CHANGED IF NECESSARY
        
        #remove last character if /
	export var="${SUB: -1}"
	if [[ "$var" == "/" ]]; then
	export SUB="${SUB:0:-1}"
	fi
	
	#Get string after last /
        SUBJID="$(echo ${SUB##*/})"
        if ((SUBJID=="fsaverage")); then #fsaverage is a symlink of freesurfer own data, should be skipped by this script
        	continue
        fi
        
        counter+=1
	main & #Launch computation in parallel
	#echo $! #Get job PID
	#jobs -l #Get job ID
	if ((counter%N_PARALLEL_COMPUTING==(N_PARALLEL_COMPUTING-1))); 
	then
		wait
	fi
done
fi

