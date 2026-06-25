#!/usr/bin/env python3
"""
keep_main_components.py

Post-processing pipeline for .nii.gz segmentation files.

Stage A — Main-component filtering
    For each non-zero label, keep only the largest 6-connected component.
    Voxels in smaller ("orphan") components are reassigned to the most common
    NON-ZERO label among their 26 immediate neighbors (using the original,
    unmodified volume as reference so the result is order-independent).
    When background (0) is the plurality winner, the runner-up non-zero label
    is used instead.  Orphans with an entirely-background 26-neighborhood
    become 0.

Stage A' — Reverse label remap  (applied immediately after Stage A)
    The compact integer labels produced by process_aparc_aseg.py Stage 3
    (1-13) are converted back to their original FreeSurfer label values so
    that the volume shares the same label space as the fuse folder.

Stage B — Label fusion  (optional, requires --fuse)
    A second segmentation folder is provided.  Files are matched by 6-digit
    subject ID found anywhere in the filename.
    Before fusion, FUSE_REMAP converts the fuse-folder labels (e.g. HOA
    atlas indices) into FreeSurfer label values so both volumes speak the
    same language.
    Merge rule: the Stage-A' result takes priority.  Where it is background
    (0), the corresponding remapped voxel from the fusion volume is inserted.

Stage C — Hole filling  (always runs)
    For each non-zero label (in sorted order), scipy binary_fill_holes fills
    any enclosed background void with that label.  Labels processed later
    overwrite earlier ones in the rare case of overlapping filled regions.

Stage D — T1 brain masking  (optional, requires --t1)
    The post-Stage-C volume is first saved as an unmasked intermediate to
    <output>/../unmasked/<filename>.  Then every label voxel whose
    corresponding T1 intensity is 0 is set to background (0), and the masked
    result is written to the normal output folder.
    T1 files are matched to segmentation files by the 6-digit subject ID
    found anywhere in the filename.
    When --t1 is omitted, the output folder still receives the unmasked
    result and the unmasked/ copy is written identically.

    Migration: on the first run of a script version that includes Stage D,
    any file already present in the output folder but absent from unmasked/ is
    automatically moved to unmasked/ (it was produced without T1 masking and
    will be reprocessed).

Usage
-----
    # Stages A + A' only, in-place:
    python keep_main_components.py /path/to/labels

    # Stages A + A' + B + C + D (all stages):
    python keep_main_components.py /path/to/labels \\
        --fuse /path/to/hoa_labels \\
        --t1   /path/to/t1_images \\
        --out  /path/to/output
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import nibabel as nib
from scipy.ndimage import label as nd_label        # connected-component labelling
from scipy.ndimage import binary_fill_holes        # morphological hole filling


# ── Connectivity constants ────────────────────────────────────────────────────

# 6-connectivity structure (face-adjacent only): diagonal touches do NOT merge.
_STRUCT_6 = np.zeros((3, 3, 3), dtype=bool)
_STRUCT_6[1, 1, 0] = _STRUCT_6[1, 1, 2] = True   # ±z
_STRUCT_6[1, 0, 1] = _STRUCT_6[1, 2, 1] = True   # ±y
_STRUCT_6[0, 1, 1] = _STRUCT_6[2, 1, 1] = True   # ±x

# 26-connectivity structure: face-, edge-, and corner-adjacent voxels.
_STRUCT_26 = np.ones((3, 3, 3), dtype=bool)

# 6 face-adjacent neighbor offsets — used for voting under 6-connectivity.
_OFFSETS_6 = [
    (1, 0, 0), (-1, 0, 0),
    (0, 1, 0), (0, -1, 0),
    (0, 0, 1), (0, 0, -1),
]

# All 26 neighbor offsets in {-1,0,1}³ \ {(0,0,0)} — used for 26-connectivity voting.
_OFFSETS_26 = [
    (di, dj, dk)
    for di in (-1, 0, 1)
    for dj in (-1, 0, 1)
    for dk in (-1, 0, 1)
    if (di, dj, dk) != (0, 0, 0)
]

# Map connectivity integer → (CC structure, neighbor offsets).
_CONNECTIVITY = {
    6:  (_STRUCT_6,  _OFFSETS_6),
    26: (_STRUCT_26, _OFFSETS_26),
}

# Regex matching the first 6-digit sequence in a filename (subject ID).
_ID_RE = re.compile(r"(\d{6})")


# ── Label-remap tables ────────────────────────────────────────────────────────

# Reverse of process_aparc_aseg.py Stage-3 LABEL_REMAP.
# Converts compact nnUNet integers (1-13) back to FreeSurfer label values.
REVERSE_REMAP = {
    1:  2,     # Left-Cerebral-White-Matter
    2:  41,    # Right-Cerebral-White-Matter
    3:  3,     # Left-Cerebral-Cortex
    4:  42,    # Right-Cerebral-Cortex
    5:  251,   # CC_Posterior
    6:  252,   # CC_Mid_Posterior
    7:  253,   # CC_Central
    8:  254,   # CC_Mid_Anterior
    9:  255,   # CC_Anterior
    10: 7,     # Left-Cerebellum-White-Matter
    11: 46,    # Right-Cerebellum-White-Matter
    12: 8,     # Left-Cerebellum-Cortex
    13: 47,    # Right-Cerebellum-Cortex
}

# Mapping from fuse-folder (HOA-style) label indices to FreeSurfer label values.
# None entries in the source list are skipped (no corresponding source label).
# Where a source label appears twice the value is identical — no conflict.
FUSE_REMAP = {
    0:  0,     # background → background
    1:  4,
    2:  43,
    3:  24,
    4:  14,
    5:  15,
    6:  63,
    7:  26,
    8:  58,
    9:  11,
    10: 50,
    11: 12,
    12: 51,
    13: 13,
    14: 52,
    15: 16,
    16: 10,
    17: 49,
    18: 5,
    19: 44,
    20: 17,
    21: 53,
    22: 18,
    23: 54,
    24: 85,
    25: 28,
    26: 60,
    27: 28,
    28: 60,
    29: 8,
    30: 47,
    31: 7,
    32: 46,
}


# ── Stage A helpers ───────────────────────────────────────────────────────────

def main_component_mask(binary_vol: np.ndarray,
                        structure: np.ndarray = _STRUCT_6) -> np.ndarray:
    """
    Return a boolean mask covering only the largest connected component
    of the non-zero region in *binary_vol*.

    Parameters
    ----------
    binary_vol : bool ndarray
    structure : ndarray (3,3,3)
        Connectivity kernel passed to scipy.ndimage.label.
        Use _STRUCT_6 (default) for 6-connectivity or _STRUCT_26 for 26.

    Returns an all-False mask when the input contains no True voxels.
    """
    if not binary_vol.any():
        return np.zeros_like(binary_vol, dtype=bool)

    labeled, n = nd_label(binary_vol, structure=structure)

    if n == 1:
        # Already a single component — nothing to filter.
        return binary_vol.astype(bool)

    # bincount index 0 is the background bin inside the CC-labeled array.
    sizes = np.bincount(labeled.ravel())
    sizes[0] = 0                        # ignore background
    main_idx = int(sizes.argmax())      # CC index with the most voxels
    return labeled == main_idx


def vote_neighbor_labels(vol: np.ndarray,
                         orphan_mask: np.ndarray,
                         offsets: list = _OFFSETS_6) -> np.ndarray:
    """
    For every True position in *orphan_mask*, return the most common
    NON-ZERO label among its immediate neighbors in *vol*.

    Background (0) is deliberately excluded from winning: if background
    has more neighbor votes than any non-zero label, the runner-up non-zero
    label is assigned instead.  Only when every neighbor is background does
    the assignment become 0.

    Implementation
    --------------
    For each offset in *offsets* the volume is roll-shifted so that
    shifted[x,y,z] holds the neighbor value of voxel (x,y,z) in that
    direction.  Wrapped edges are zeroed.  This gives a
    (n_orphans × len(offsets)) neighbor matrix with no Python loop over
    voxels.  A broadcasted equality comparison builds per-label vote counts.

    Parameters
    ----------
    vol : int ndarray (X, Y, Z)
        Original label volume — read-only, used for neighbor lookup.
    orphan_mask : bool ndarray (X, Y, Z)
        True at every voxel to be reassigned.
    offsets : list of (int, int, int)
        Neighbor offsets to consider.  Use _OFFSETS_6 (default, 6-connectivity)
        or _OFFSETS_26 (26-connectivity).

    Returns
    -------
    assignments : int ndarray (n_orphans,)
        Winner label for each orphan (same order as np.where(orphan_mask)).
        Zero only when all neighbors are background.
    """
    orphan_coords = np.where(orphan_mask)
    n_orphans = orphan_coords[0].size

    if n_orphans == 0:
        return np.array([], dtype=vol.dtype)

    # ── Gather neighbor labels for every orphan voxel at once ─────────────────
    # neighbor_matrix[i, k] = label of the k-th neighbor of orphan i.
    neighbor_matrix = np.zeros((n_orphans, len(offsets)), dtype=vol.dtype)

    for col, (di, dj, dk) in enumerate(offsets):
        # np.roll(vol, shift=(-di,-dj,-dk))[x,y,z] == vol[x+di, y+dj, z+dk]
        # (with modular wrap-around, which we zero out below).
        shifted = np.roll(vol, shift=(-di, -dj, -dk), axis=(0, 1, 2))

        # Nullify the wrapped-around edges so boundary voxels do not gain
        # phantom neighbors from the opposite face of the volume.
        if di > 0:
            shifted[-di:, :, :] = 0
        elif di < 0:
            shifted[:-di, :, :] = 0
        if dj > 0:
            shifted[:, -dj:, :] = 0
        elif dj < 0:
            shifted[:, :-dj, :] = 0
        if dk > 0:
            shifted[:, :, -dk:] = 0
        elif dk < 0:
            shifted[:, :, :-dk] = 0

        neighbor_matrix[:, col] = shifted[orphan_coords]

    # ── Vote: count non-zero neighbor labels, pick the plurality winner ───────
    # Excluding background (0) from the candidate set means: even when
    # background has more votes than any single non-zero label, we still assign
    # the best non-zero label (the "second most common neighbor overall").
    nonzero_labels = np.unique(vol)
    nonzero_labels = nonzero_labels[nonzero_labels != 0]

    if nonzero_labels.size == 0:
        return np.zeros(n_orphans, dtype=vol.dtype)

    # vote_matrix[i, k] = how many of orphan i's neighbors carry nonzero_labels[k].
    # Shape: (n_orphans, n_nonzero_labels).  Built via numpy broadcast; no voxel loop.
    vote_matrix = (
        neighbor_matrix[:, :, np.newaxis]             # (n_orphans, n_offsets, 1)
        == nonzero_labels[np.newaxis, np.newaxis, :]  # (1, 1, n_labels)
    ).sum(axis=1)                                     # sum over neighbors → (n_orphans, n_labels)

    winner_idx  = vote_matrix.argmax(axis=1)
    assignments = nonzero_labels[winner_idx]

    # Orphans whose entire 26-neighborhood is background → set to 0.
    all_background = (neighbor_matrix != 0).sum(axis=1) == 0
    assignments[all_background] = 0

    return assignments


# ── Shared remap helper ───────────────────────────────────────────────────────

def remap_volume(vol: np.ndarray, remap: dict) -> np.ndarray:
    """
    Apply *remap* to *vol* and return a new array.

    Labels present in *remap* are translated to their mapped value.
    Labels absent from *remap* become 0 (background).
    A fresh output array avoids in-place collisions when a source value
    equals a destination value of a different entry.

    Parameters
    ----------
    vol : int ndarray
        Input label volume (not modified).
    remap : dict {int → int}
        Label translation table.

    Returns
    -------
    out : int ndarray, same shape and dtype as *vol*.
    """
    out = np.zeros_like(vol)
    for src, dst in remap.items():
        out[vol == src] = dst
    return out


# ── Stage B helper ────────────────────────────────────────────────────────────

def fuse_labels(primary: np.ndarray, secondary: np.ndarray) -> np.ndarray:
    """
    Merge two label volumes.

    *primary* takes full priority: its non-zero labels are kept unchanged.
    *secondary* fills only the background (0) voxels of *primary*.

    Both arrays must have the same shape.

    Returns
    -------
    fused : int ndarray, same shape and dtype as *primary*.
    """
    if primary.shape != secondary.shape:
        raise ValueError(
            f"Shape mismatch: primary {primary.shape} vs secondary {secondary.shape}"
        )
    fused = primary.copy()
    background_mask = primary == 0
    fused[background_mask] = secondary[background_mask]
    return fused


# ── Stage C helper ────────────────────────────────────────────────────────────

def fill_label_holes(vol: np.ndarray) -> np.ndarray:
    """
    For each non-zero label, fill any enclosed background void (3-D hole)
    with that label.

    scipy.ndimage.binary_fill_holes identifies background regions that are
    completely surrounded in 3-D by the label's mask and sets them to True.
    New voxels (previously background) are assigned the filling label,
    overwriting whatever was there.

    Labels are processed in sorted order; in the unlikely event that filled
    regions of two labels overlap, the higher-valued label wins.

    Returns
    -------
    filled : int ndarray, same shape and dtype as *vol*.
    """
    filled = vol.copy()
    unique_labels = np.unique(vol)
    unique_labels = unique_labels[unique_labels != 0]   # skip background

    for lbl in unique_labels:
        mask        = vol == lbl
        mask_filled = binary_fill_holes(mask)           # True inside enclosed voids too
        new_voxels  = mask_filled & ~mask               # only the newly filled positions
        filled[new_voxels] = lbl

    return filled


# ── File-matching utility ─────────────────────────────────────────────────────

def extract_subject_id(path: Path) -> Optional[str]:
    """Return the first 6-digit sequence found in *path*'s filename, or None."""
    m = _ID_RE.search(path.name)
    return m.group(1) if m else None


def find_fuse_file(fuse_folder: Path, subject_id: str) -> Optional[Path]:
    """
    Return the .nii.gz file in *fuse_folder* whose name contains *subject_id*.
    Warns and returns the first match when multiple files share the same ID.
    Returns None when no match is found.
    """
    matches = [
        p for p in sorted(fuse_folder.glob("*.nii.gz"))
        if extract_subject_id(p) == subject_id
    ]
    if not matches:
        return None
    if len(matches) > 1:
        print(
            f"  WARNING: {len(matches)} fuse files share ID {subject_id} "
            f"— using {matches[0].name}"
        )
    return matches[0]


def find_t1_file(t1_folder: Path, subject_id: str) -> Optional[Path]:
    """Return the .nii.gz T1 file in *t1_folder* whose name contains *subject_id*."""
    matches = [
        p for p in sorted(t1_folder.glob("*.nii.gz"))
        if extract_subject_id(p) == subject_id
    ]
    if not matches:
        return None
    if len(matches) > 1:
        print(
            f"  WARNING: {len(matches)} T1 files share ID {subject_id} "
            f"— using {matches[0].name}"
        )
    return matches[0]


def _stage_d_mask_and_save(
    name_ref: Path,
    vol: np.ndarray,
    affine,
    header,
    dst: Path,
    t1_folder: Optional[Path],
) -> None:
    """Apply T1 brain mask to *vol* (in-place) and save to *dst*.

    *name_ref* is any path whose filename contains the 6-digit subject ID used
    to locate the matching T1.  When *t1_folder* is None the volume is saved
    unchanged.
    """
    if t1_folder is not None:
        sid     = extract_subject_id(name_ref)
        t1_path = find_t1_file(t1_folder, sid) if sid else None
        if t1_path is None:
            print(
                f"  WARNING: no T1 found for ID {sid} in {t1_folder} "
                f"— T1 masking skipped."
            )
        else:
            print(f"  Stage D — T1 masking with {t1_path.name} ...", flush=True)
            t1_img = nib.load(str(t1_path))
            t1_vol = np.asarray(t1_img.dataobj, dtype=t1_img.get_data_dtype())
            vol[t1_vol == 0] = 0
    nib.save(nib.Nifti1Image(vol, affine, header), str(dst))
    print(f"  Saved  →  {dst}", flush=True)


# ── Per-file pipeline ─────────────────────────────────────────────────────────

def process_file(
    src: Path,
    dst: Path,
    unmasked_dst: Path,
    fuse_folder: Optional[Path] = None,
    t1_folder: Optional[Path] = None,
    connectivity: int = 6,
    reassign_orphans: bool = False,
) -> None:
    """
    Run the pipeline on *src* and write the result to *dst*.

    Stage A  : main-component filtering; orphans → 0 by default, or
               reassigned to most common non-zero neighbor (--reassign-orphans)
    Stage A' : reverse remap back to FreeSurfer label values (REVERSE_REMAP)
    Stage B  : fuse-volume remap (FUSE_REMAP) + label fusion (skipped when fuse_folder is None)
    Stage C  : per-label hole filling (always runs)
    Stage D  : save unmasked intermediate to *unmasked_dst*, then apply T1
               brain mask (voxels where T1 == 0 → label 0) and save to *dst*.
               If *t1_folder* is None, *dst* receives the same data as
               *unmasked_dst*.

    Parameters
    ----------
    unmasked_dst : Path
        Destination for the post-Stage-C result before T1 masking.
    t1_folder : Path or None
        Folder containing T1 .nii.gz files matched by 6-digit subject ID.
        When None, T1 masking is skipped.
    connectivity : {6, 26}
        Connectivity for CC detection and neighbor voting.
        6 = face-adjacent only (default); 26 = face + edge + corner adjacent.
    reassign_orphans : bool
        If False (default), non-main-component voxels are set to 0.
        If True, they are reassigned to the most common non-zero neighbor label.
    """
    print(f"  Loading   {src.name} ...", flush=True)
    img = nib.load(str(src))
    vol = np.asarray(img.dataobj, dtype=np.int32)

    unique_labels = np.unique(vol)
    unique_labels = unique_labels[unique_labels != 0]

    if unique_labels.size == 0:
        print(f"  WARNING: no non-zero labels — skipping.")
        return

    # ── Stage A: keep only the main component per label ───────────────────────
    cc_struct, nb_offsets = _CONNECTIVITY[connectivity]
    print(f"  Stage A — main-component filtering (connectivity={connectivity}) ...",
          flush=True)

    # Build clean_vol: start from zeros, copy only the main-component voxels of
    # each label.  This guarantees that clean_vol contains exactly one connected
    # component per label — no orphan can survive.
    # orphan_mask is computed in parallel for logging and for --reassign-orphans.
    clean_vol    = np.zeros_like(vol)
    orphan_mask  = np.zeros(vol.shape, dtype=bool)
    for lbl in unique_labels:
        binary    = vol == lbl
        main_mask = main_component_mask(binary, structure=cc_struct)
        clean_vol[main_mask] = lbl          # only main component goes into output
        orphan_mask |= binary & ~main_mask  # non-main voxels of this label

    n_orphans = int(orphan_mask.sum())
    print(f"  Labels: {list(unique_labels)}  |  Orphan voxels: {n_orphans}", flush=True)

    if reassign_orphans and n_orphans > 0:
        # Reassign each orphan to the most common non-zero neighbor label.
        # Neighbor lookup uses the *original* vol (order-independent reference).
        assignments = vote_neighbor_labels(vol, orphan_mask, offsets=nb_offsets)
        clean_vol[orphan_mask] = assignments

    vol = clean_vol  # guaranteed: one main component per label (+ reassignments if any)

    # ── Stage A': reverse remap → FreeSurfer label space ──────────────────────
    # Converts compact nnUNet integers (1-13) produced by process_aparc_aseg.py
    # Stage 3 back to their original FreeSurfer values so that the volume shares
    # the same label space as the fuse folder.
    print(f"  Stage A' — reverse remap to FreeSurfer labels ...", flush=True)
    vol = remap_volume(vol, REVERSE_REMAP)

    # ── Stage B: remap fuse labels, then fuse ─────────────────────────────────
    if fuse_folder is not None:
        print(f"  Stage B — label fusion ...", flush=True)
        sid       = extract_subject_id(src)
        fuse_path = find_fuse_file(fuse_folder, sid) if sid else None

        if fuse_path is None:
            print(
                f"  WARNING: no fuse file found for ID {sid} in {fuse_folder} "
                f"— Stage B skipped."
            )
        else:
            fuse_img = nib.load(str(fuse_path))
            fuse_raw = np.asarray(fuse_img.dataobj, dtype=np.int32)
            # Convert fuse-folder labels (HOA indices) to FreeSurfer values
            # before merging so both volumes share the same label space.
            fuse_vol = remap_volume(fuse_raw, FUSE_REMAP)
            # Fill enclosed holes per label in the fuse volume so that no
            # internal voids remain before the two volumes are merged.
            fuse_vol = fill_label_holes(fuse_vol)
            vol = fuse_labels(vol, fuse_vol)
            print(f"  Fused with {fuse_path.name}", flush=True)

    # ── Stage C: fill holes per label ─────────────────────────────────────────
    print(f"  Stage C — hole filling ...", flush=True)
    vol = fill_label_holes(vol)

    # ── Stage D: save unmasked intermediate, then apply T1 brain mask ─────────
    nib.save(nib.Nifti1Image(vol, img.affine, img.header), str(unmasked_dst))
    print(f"  Saved (unmasked) → {unmasked_dst}", flush=True)
    _stage_d_mask_and_save(src, vol, img.affine, img.header, dst, t1_folder)


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Segmentation post-processing: "
            "(A) keep main CC per label + reassign orphans, "
            "(A') reverse remap to FreeSurfer label space, "
            "(B) optional HOA-label fusion, "
            "(C) optional per-label hole filling."
        )
    )
    parser.add_argument(
        "folder",
        type=Path,
        help="Folder containing .nii.gz segmentation files to process.",
    )
    parser.add_argument(
        "--fuse",
        type=Path,
        default=None,
        metavar="FUSE_FOLDER",
        help=(
            "Optional second label folder for Stage B fusion.  Files are "
            "matched by 6-digit subject ID.  The primary (Stage-A) labels "
            "take priority; the fuse labels fill background voxels only."
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        metavar="OUTPUT_FOLDER",
        help=(
            "Destination folder for the final files.  "
            "Defaults to in-place overwrite of the source folder."
        ),
    )
    parser.add_argument(
        "--reassign-orphans",
        action="store_true",
        default=False,
        help=(
            "Stage A: reassign non-main-component voxels to their most common "
            "non-zero neighbor label instead of setting them to 0 (default: off)."
        ),
    )
    parser.add_argument(
        "--connectivity",
        type=int,
        choices=[6, 26],
        default=6,
        metavar="{6,26}",
        help=(
            "Connectivity for Stage A CC detection and neighbor voting.  "
            "6 = face-adjacent only (default); 26 = face + edge + corner."
        ),
    )
    parser.add_argument(
        "--t1",
        type=Path,
        default=None,
        metavar="T1_FOLDER",
        help=(
            "Optional folder of T1 .nii.gz images matched by 6-digit subject ID.  "
            "Stage D zeros out label voxels where the T1 is 0 (brain mask).  "
            "The pre-masking volume is always saved to <output>/../unmasked/."
        ),
    )
    parser.add_argument(
        "--recompute",
        action="store_true",
        default=False,
        help=(
            "Recompute and overwrite output files that already exist.  "
            "By default, files whose destination already exists are skipped."
        ),
    )
    args = parser.parse_args()

    src_folder = args.folder.resolve()
    if not src_folder.is_dir():
        print(f"ERROR: folder not found: {src_folder}", file=sys.stderr)
        sys.exit(1)

    fuse_folder = args.fuse.resolve() if args.fuse else None
    if fuse_folder is not None and not fuse_folder.is_dir():
        print(f"ERROR: fuse folder not found: {fuse_folder}", file=sys.stderr)
        sys.exit(1)

    t1_folder = args.t1.resolve() if args.t1 else None
    if t1_folder is not None and not t1_folder.is_dir():
        print(f"ERROR: T1 folder not found: {t1_folder}", file=sys.stderr)
        sys.exit(1)

    out_folder = args.out.resolve() if args.out else src_folder
    if out_folder != src_folder:
        out_folder.mkdir(parents=True, exist_ok=True)

    unmasked_folder = out_folder.parent / "unmasked"
    unmasked_folder.mkdir(parents=True, exist_ok=True)

    files = sorted(src_folder.glob("*.nii.gz"))
    if not files:
        print(f"No .nii.gz files found in {src_folder}")
        sys.exit(0)

    print(f"Found {len(files)} file(s) in {src_folder}")
    if fuse_folder:
        print(f"Fuse    → {fuse_folder}")
    if t1_folder:
        print(f"T1      → {t1_folder}")
    print(f"Output  → {'in-place' if out_folder == src_folder else out_folder}")
    print(f"Unmasked → {unmasked_folder}\n")

    for i, src in enumerate(files, 1):
        dst          = out_folder     / src.name
        unmasked_dst = unmasked_folder / src.name
        print(f"[{i}/{len(files)}] {src.name}")

        # Migrate files produced without T1 masking (out_folder has the file but
        # unmasked_folder does not).  Move the old result to unmasked/, then apply
        # Stage D only — skips the expensive stages A–C.
        if dst.exists() and not unmasked_dst.exists():
            dst.rename(unmasked_dst)
            print(f"  Migrated pre-masking output → {unmasked_dst}")
            _img = nib.load(str(unmasked_dst))
            _vol = np.asarray(_img.dataobj, dtype=np.int32)
            _stage_d_mask_and_save(unmasked_dst, _vol, _img.affine, _img.header,
                                   dst, t1_folder)
            print()
            continue

        if dst.exists() and not args.recompute:
            print(f"  Skipped (output exists — use --recompute to overwrite)\n")
            continue
        process_file(src, dst, unmasked_dst,
                     fuse_folder=fuse_folder,
                     t1_folder=t1_folder,
                     connectivity=args.connectivity,
                     reassign_orphans=args.reassign_orphans)
        print()

    print("Done.")


if __name__ == "__main__":
    main()
