#!/usr/bin/env python3
"""
apply_gwmatter_crosswalk.py

Remap any GWmatter-scheme label volume (35 classes, 0-34, as produced by the
GWmatter model -- both _seg_pred_GWmatter.nii.gz and _seg_gt_GWmatter.nii.gz
share this same scheme) into FreeSurfer aseg label space, using the static,
already-validated crosswalk in
  ribbon_project/GWmatter_predictions/gwmatter35_to_fs.tsv

That crosswalk was derived empirically for subject 792654 (2026-09-16) via
spatial overlap against an independent HOA_Subcortical_Labels reference
(gwmatter_to_aseg.py, in 792654_edited_0925/) -- all 34 classes' derived FS
ids matched colortable_35class_GWmatter.txt's declared identity exactly, so
this static table is reused directly here rather than re-deriving it per
subject: these test6_withGT subjects have no HOA_Subcortical_Labels-style
reference file to derive against anyway (only pred/gt GWmatter-scheme
volumes), and the crosswalk is a property of the model's own output scheme,
not subject-specific anatomy.

Unlike gwmatter_to_aseg.py (which reslices onto a *different* reference
grid), this script assumes the input is already on the T1's own grid (both
_seg_pred_GWmatter.nii.gz and _seg_gt_GWmatter.nii.gz share the T1's grid
for these test6_withGT subjects) -- pass --reslice if that's ever not true.

No dilation/erosion, no carve/refill (same reasoning as gwmatter_to_aseg.py:
a dense whole-brain segmentation has no holes to fill).
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

DEFAULT_TSV = Path(__file__).parent.parent / "GWmatter_predictions" / "gwmatter35_to_fs.tsv"


def load_crosswalk(tsv_path: Path) -> dict[int, int]:
    mapping: dict[int, int] = {}
    with tsv_path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            mapping[int(row["gw_class"])] = int(row["fs_id"])
    if not mapping:
        sys.exit(f"ERROR: no mapping rows read from {tsv_path}")
    return mapping


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gwmatter", required=True, help="a _seg_pred_GWmatter.nii.gz or _seg_gt_GWmatter.nii.gz file")
    ap.add_argument("--tsv", default=str(DEFAULT_TSV))
    ap.add_argument("--reference", default=None, help="reslice onto this grid first (nearest); omit if already matching")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    crosswalk = load_crosswalk(Path(args.tsv))

    gw_img = nib.load(args.gwmatter)
    if args.reference is not None:
        from nibabel.processing import resample_from_to

        ref_img = nib.load(args.reference)
        if gw_img.shape != ref_img.shape or not np.allclose(gw_img.affine, ref_img.affine):
            gw_img = resample_from_to(gw_img, (ref_img.shape, ref_img.affine), order=0)
        out_affine, out_header = ref_img.affine, ref_img.header
    else:
        out_affine, out_header = gw_img.affine, gw_img.header

    gw = np.asarray(gw_img.dataobj).astype(np.int32)
    present = set(int(v) for v in np.unique(gw) if v != 0)
    unmapped = sorted(present - set(crosswalk.keys()))
    if unmapped:
        print(f"[warn] GWmatter classes present but not in crosswalk (left as background 0): {unmapped}")

    out = np.zeros_like(gw, dtype=np.int32)
    for g, fs_id in crosswalk.items():
        out[gw == g] = fs_id

    out_img = nib.Nifti1Image(out, out_affine, out_header)
    out_img.header.set_data_dtype(np.int32)
    nib.save(out_img, args.out)
    print(f"[ok] {args.gwmatter} -> {args.out} ({len(present & set(crosswalk.keys()))}/{len(present)} classes mapped)")


if __name__ == "__main__":
    main()
