# Legacy tools

Earlier tools, kept unchanged and not maintained: `ribbon_edit_script.sh` and `ribbon_edit_script_only_subcortical.sh` (ribbon +
subcortical files as input; the tool archived on Zenodo), `keep_main_components.py` and `mask_test_images.py` (aparc+aseg / nnUNet
preprocessing). The current tools are in the parent folder, see `../README.md`.

## Legacy: ribbon_edit_script.sh

`legacy/ribbon_edit_script.sh` reuses the FreeSurfer 7.4.1 recon-all pipeline to correct the white and pial surfaces based on manually
curated ribbon and subcortical files. Its help text, unchanged:

```
AUTHOR: Benoît Verreman

LAST UPDATE: 2025-10-08
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
-t 10: Non hemisphere specific files to create
-t 11: (aseg+aparc) Compute stats and other files
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
```

### Comparison: merge_hoa_into_aseg.py vs ribbon_edit_script.sh

- `ribbon_edit_script.sh` is a full FreeSurfer surface reconstruction pipeline: it takes an edited ribbon + subcortical volume and drives
  about 30 FreeSurfer tools to produce white/pial surfaces, parcellations and stats. HOA merging is one step inside a larger pipeline.
- `merge_hoa_into_aseg.py` is a single-purpose volume editing tool: it takes an existing aseg, replaces its subcortical labels with HOA and
  writes a new aseg. No surfaces, no FreeSurfer tool calls.

| | `ribbon_edit_script.sh` | `merge_hoa_into_aseg.py` |
|---|---|---|
| Label mapping | Hardcoded bash arrays (`LABELS_SUBCORTICAL_OLD/NEW`) | External TSV file (`--tsv`) |
| Strategy | Remap, dilate/erode HA into the ribbon, merge into aseg | Remap, carve FS labels, paint HOA labels, refill gaps |
| HA morphology | Dilation/erosion applied | None: the HOA boundary is used as is |
| Carved-but-unpainted voxels | WM fills where the subcortical label was removed (`edit_aseg_presurf_based_on_ribbon.py`) | Configurable: `surround` (nearest non-replaced tissue via EDT), `nearest` (regrows the nucleus), `none` |
| Cerebellum | Always kept from the FS ribbon | Optional via `--include-cerebellum` |
| VDC anterior/posterior | Collapsed to 28/60 | Optionally split into 28/60 + 8028/8060 via `--split-vdc` |
| Reslicing | Handled externally before the script | Built in, `--reslice` (nibabel) |

Only `ribbon_edit_script.sh`: morphologically integrates HA (hippocampus, amygdala, inferior lateral ventricle) into the cortical ribbon with
dilation/erosion, which matters for surface topology and not only for volume labels; recomputes the white and pial surfaces from the edited
anatomy; produces `aseg.presurf_wo_subc.mgz` (subcortical replaced by WM), used as a constraint during surface placement.

Only `merge_hoa_into_aseg.py`: resolution-agnostic (HOA may be on another grid than the aseg); TSV-driven mapping, portable, no hardcoded label
tables; the EDT-based `surround` refill gives a tighter HOA boundary than the WM-fill fallback of the shell script.

The two can be combined: `merge_hoa_into_aseg.py` prepares `aseg.presurf.mgz`, then `aseg_surface_script.sh` (or `ribbon_edit_script.sh`) drives the surface
reconstruction from it.
