# FreeSurfer Surface editing
This repository contains:
- ribbon_edit_script.sh: a shell script that reuses Freesurfer 7.4.1 pipeline to correct the pial surface based on ribbon and subcortical files.  \
  Examples:  \
  bash ribbon_edit_script.sh --help  \
  bash ribbon_edit_script.sh --subjid <subjid> --ribbon <ribbon-edit.nii.gz> --subcortical <subcortical.nii.gz>
