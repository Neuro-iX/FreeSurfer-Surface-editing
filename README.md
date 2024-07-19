# FreeSurfer Surface editing

[![DOI](https://zenodo.org/badge/665245427.svg)](https://zenodo.org/doi/10.5281/zenodo.12784724)

This repository contains:
- **ribbon_edit_script.sh**: \
This shell script reuses Freesurfer 7.4.1 recon-all pipeline to correct the white and pial surfaces based on manually curated ribbon and subcortical files.

Help function:

AUTHOR: Benoît Verreman

LAST UPDATE: 2024-04-16

DESCRIPTION: \
Use ribbon and subcortical NIFTI files to recompute pial surface, \
based on previously created <subjid> folder using Freesurfer 7.4.1 \
Create a log 'report.sh'. \
Create a folder 'outputs' with all the output files.

PREREQUISITES: \
*Freesurfer variables: \
Export SUBJECTS_DIR and FREESURFER_HOME correctly

*Python script: \
Test if you have access to python: 'which python' \
Install two python libraries: 'pip install nibabel scipy' \
OR use conda environment: \
conda create --name env_ribbon_edit_script \
conda activate env_ribbon_edit_script \
conda install nibabel scipy -c conda-forge \
conda list | grep -E 'nibabel|scipy'

*Freesurfer 'raw' output folder: \
If you want to launch Freesurfer 7.4.1 recon-all pipeline using the script: \
	Add argument -i \
Else: \
	Launch Freesurfer 7.4.1 command before using script:  \
$ recon-all -s <subjid> -i <subject_image> -autorecon1 -autorecon2 -hires -parallel -openmp 4 -expert expert_file.txt \
Prepare ribbon and subcortical NIFTI files (for step 0) \
Put the script inside <subjid> folder

*Specific labels and statical values: \
Modify the value of some constants in the script if needed: LABELS_SUBCORTICAL, LABEL_RIBBON_WM_LH, LABEL_RIBBON_WM_RH, PIAL_BORDER_LOW

EXAMPLES: \
$ bash ribbon_edit_script.sh -i 133019_T1w_acpc_dc_restore.nii.gz -s 133019 -b 133019_ribbon.nii.gz -c 133019_subcortical.nii.gz

$ bash ribbon_edit_script.sh -s <subjid> -t 8 -r #start from pial computation step (8), right hemisphere only (r)

PARAMETERS:

HELP \
-h: Print this string, and exit

INPUT FILES \
-i: Relative or absolute path to T1w image file \
-s: Relative or absolute path to subjid folder (Necessary) \
-b: Relative or absolute path to ribbon file \
-c: Relative or absolute path to subcortical file

TAG \
-t 0: (ribbons) Start with resizing RIBBON_EDIT and SUBCORTICAL \
-t 1: (bmask) Start with BRAIN_MASK \
-t 2: (maskT1) Start with T1_MASKED \
-t 3: (brain.finalsurfs) Start with skull-stripping up to BRAIN_FINALSURFS \
-t 4: (wm-bmask) Start the creation of WM_BMASK based on RIBBON_EDIT \
-t 5: (wm) Start from computing WM based on WM_BMASK \
-t 6: (orig) Start from computing orig surface based on wm from RIBBON_EDIT \
-t 7: (stats) Start from computing stats \
-t 8: (brain.finalsurfs-edit) edit brain.finalsurfs with GM from RIBBON_EDIT \
-t 9: (pial) Start from computing pial surface \
-t 10: (smooth) Start from smoothing pial surface 

VALUES \
-p: Give value of PIAL_BORDER_LOW

HEMI \
-r: Compute only right hemisphere surface  \
-l: Compute only left hemisphere surface

RESET \
-d: Reset outputs folder and report.sh script
