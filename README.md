# FreeSurfer Surface editing

[![DOI](https://zenodo.org/badge/665245427.svg)](https://zenodo.org/doi/10.5281/zenodo.12784724)

Rebuild FreeSurfer 7.4.1 white and pial surfaces (and the files that depend on them) from a **manually curated or predicted
segmentation**, instead of from FreeSurfer's own tissue boundaries. The segmentation is given as an aseg in FreeSurfer label space.

Typical workflow, for example with the output of the HOA-2 segmentation model (https://github.com/HOA-2/):

```
model segmentation  ->  label ids remapped to FreeSurfer's  ->  aseg_surface_script.sh  ->  FreeSurfer 7.4.1-style outputs
(HOA-2 labels)          (apply_gwmatter_crosswalk.py or          (surfaces, labels, stats)
                         merge_hoa_into_aseg.py)
```

## Layout

Top level (the current tools):

| File | Role |
|---|---|
| `aseg_surface_script.sh` | aseg (FreeSurfer labels) + T1 -> surfaces. Recomputes the recon-all stages that depend on the segmentation (steps `-t 0..11`). |
| `merge_hoa_into_aseg.py` | Replace the subcortical labels of an existing aseg with HOA labels (remap, carve, paint, refill). |
| `apply_gwmatter_crosswalk.py` | Remap a label volume to FreeSurfer ids with a TSV crosswalk. |

`legacy/` holds earlier tools, kept unchanged and not maintained: `ribbon_edit_script.sh` and `ribbon_edit_script_only_subcortical.sh`
(ribbon + subcortical files as input; the tool archived on Zenodo, documented in [`legacy/README.md`](legacy/README.md)),
`keep_main_components.py` and `mask_test_images.py` (aparc+aseg / nnUNet preprocessing).

## aseg_surface_script.sh

Uses an aseg in **FreeSurfer label space** and a T1 image to recompute the surfaces, reusing FreeSurfer 7.4.1 intermediate files.
The ribbon (labels 2, 3, 41, 42) is derived automatically from the aseg. The script writes a log `report.sh` and an `outputs/` folder
inside `$SUBJECTS_DIR/<subjid>/`.

### Prerequisites

- FreeSurfer 7.4.1, with `FREESURFER_HOME` set and `SUBJECTS_DIR` exported.
- Python with `nibabel` and `scipy` (`pip install nibabel scipy`, or
  `conda create --name env_ribbon_edit_script python=3.10 nibabel scipy -c conda-forge`).
- A recon-all folder for the subject, `$SUBJECTS_DIR/<subjid>/<subjid>_freesurfer`. Either let the script run recon-all (option `-i`), or run it
  first: `recon-all -s <subjid>_freesurfer -i <T1> -autorecon1 -autorecon2 -autorecon3 -hires -parallel -openmp 4`.
- An aseg file with FreeSurfer label values: white matter 2 (lh) / 41 (rh), cortex 3 / 42, subcortical structures
  (4, 5, 7, 8, 10-18, 24, 26, 28, 43-47, 49-54, 58, 60, 77, 85, ...) and midline labels such as 251-255.

### Examples

```bash
export SUBJECTS_DIR=`pwd`

# T1 + aseg, with recon-all run by the script
bash aseg_surface_script.sh -i <T1.nii.gz> -s <subjid> -a <aseg.mgz>

# recon-all folder already there
bash aseg_surface_script.sh -s <subjid> -a <aseg.mgz>

# restart from the orig surface (step 4), right hemisphere only
bash aseg_surface_script.sh -s <subjid> -t 4 -r

# whole folder: T1 detected by '*T1*', aseg by '*aseg*'
bash aseg_surface_script.sh -f <folder> -i '*' -a '*'
```

### Options

| Option | Meaning |
|---|---|
| `-s` | Subject id / subfolder name (required) |
| `-a` | Aseg file in FreeSurfer label space |
| `-i` | T1 image; also runs recon-all |
| `-t N` | Start at step N (see below) |
| `-p` | Value of `PIAL_BORDER_LOW` (default 5) |
| `-l` / `-r` | Left / right hemisphere only |
| `-f` | Run on every subfolder of the given folder (each must contain an aseg, and a T1 if `-i` is used) |
| `-k` | Reset the `outputs/` folder and `report.sh` |
| `-h` | Help |

### Steps (`-t`)

| Step | Name | What it does |
|---|---|---|
| 0 | aseg | Pad the aseg, extract the ribbon (labels 2, 3, 41, 42) |
| 1 | bmask | Brain mask from the aseg, mask the T1 |
| 2 | aseg.presurf | Create `aseg.presurf` and `aseg.presurf_wo_subc` |
| 3 | wm | White matter from the ribbon, pretess |
| 4 | orig | Orig surface tessellation and topology fix |
| 5 | bf-edit | Edit `brain.finalsurfs` with the cortex from the ribbon |
| 6 | stats | Stats, cortex label, sphere, aparc |
| 7 | pial | Pial surface (`mris_place_surface`) |
| 8 | non-hemi | Non-hemisphere files (copy, denoise, segment) |
| 9 | autorecon3 | Curvature and thickness stats |
| 10 | aparc+aseg | Cortical ribbon, aparc+aseg, wmparc, segstats |
| 11 | bonus | Extended parcellations (BA_exvivo, vpnl) |

## merge_hoa_into_aseg.py

Replaces FreeSurfer's automated subcortical labels in an aseg (or aparc+aseg) with an HOA label map: remap the HOA values with a TSV,
carve the replaced FreeSurfer labels, paint the HOA labels, refill the voxels that were carved but not repainted.

```bash
python merge_hoa_into_aseg.py --hoa <hoa.nii.gz> --aseg <aseg.mgz> --tsv <mapping.tsv> --out <out.mgz> \
    [--brainmask <mask>] [--reslice] [--refill surround|nearest|none] [--split-vdc] [--include-cerebellum]
```

The TSV has the columns `hoa_value` and `fs_id`. See the docstring of the script for the refill modes.

## apply_gwmatter_crosswalk.py

Remaps a label volume into FreeSurfer aseg label space with a TSV crosswalk (columns `gw_class` and `fs_id`). No dilation, erosion, carve or
refill: meant for dense whole-brain segmentations.

```bash
python apply_gwmatter_crosswalk.py --gwmatter <labels.nii.gz> --tsv <crosswalk.tsv> [--reference <T1.nii.gz>] --out <aseg.nii.gz>
```

Pass `--tsv` explicitly. `--reference` reslices the input onto that grid first (nearest neighbour); leave it out if the grids already match.
