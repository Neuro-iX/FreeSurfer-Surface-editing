# FreeSurfer Surface editing

[![DOI](https://zenodo.org/badge/665245427.svg)](https://zenodo.org/doi/10.5281/zenodo.12784724)

This repository contains:
- **ribbon_edit_script.sh**: \
This shell script reuses Freesurfer 7.4.1 recon-all pipeline to correct the white and pial surfaces based on manually curated ribbon and subcortical files.

Help function:

AUTHOR: Benoît Verreman

LAST UPDATE: 2024-06-25

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
