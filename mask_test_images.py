#!/usr/bin/env python3
"""
mask_test_images.py

Stage 6 of process_aparc_aseg.py, extracted as a standalone script.

For each subject found in ./SubcorticalParcellations/dseg/, load the HOA
label mask (<id>_dseg.nii.gz) and the corresponding T1 from
./imagesTs_not_masked_yet/<id>_0000.nii.gz.

HOA cerebellum labels {29, 30, 31, 32} are neutralised (set to 0) before
masking so that cerebellum T1 voxels are preserved.
All remaining non-zero HOA voxels are zeroed out in the T1.
Output is saved to ./imagesTs/<id>_0000.nii.gz.

Output is skipped when it already exists (idempotent).
"""

import re
import sys
import numpy as np
import nibabel as nib
from pathlib import Path

# HOA cerebellum labels excluded from masking (T1 voxels are kept)
HOA_CEREB_LABELS = frozenset([29, 30, 31, 32])


def run_stage_6(subject_id: str,
                labels_hoa_dir: Path,
                t1_src_dir: Path,
                t1_dst_dir: Path) -> None:
    mask_path = labels_hoa_dir / f"{subject_id}_dseg.nii.gz"
    t1_path   = t1_src_dir     / f"{subject_id}_0000.nii.gz"
    out_path  = t1_dst_dir     / f"{subject_id}_0000.nii.gz"

    if out_path.exists():
        print(f"[{subject_id}] Already done — skipping.")
        return

    if not mask_path.exists():
        print(f"[{subject_id}] WARNING: HOA mask not found at {mask_path} — skipping.")
        return

    if not t1_path.exists():
        print(f"[{subject_id}] WARNING: T1 not found at {t1_path} — skipping.")
        return

    print(f"[{subject_id}] Masking test T1 ...")
    mask_img = nib.load(str(mask_path))
    t1_img   = nib.load(str(t1_path))

    mask = np.asarray(mask_img.dataobj, dtype=np.int32)
    t1   = np.asarray(t1_img.dataobj,   dtype=t1_img.get_data_dtype())

    # Treat cerebellum labels as background so those T1 voxels are preserved.
    effective_mask = mask.copy()
    effective_mask[np.isin(effective_mask, list(HOA_CEREB_LABELS))] = 0

    masked_t1 = np.where(effective_mask == 0, t1, 0).astype(t1_img.get_data_dtype())

    nib.save(nib.Nifti1Image(masked_t1, t1_img.affine, t1_img.header), str(out_path))
    print(f"[{subject_id}] Saved → {out_path}")


def main() -> None:
    script_dir    = Path(__file__).parent.resolve()
    labels_hoa_dir = script_dir / "SubcorticalParcellations" / "dseg"
    t1_src_dir     = script_dir / "imagesTs_not_masked_yet"
    t1_dst_dir     = script_dir / "imagesTs"

    if not labels_hoa_dir.is_dir():
        print(f"ERROR: HOA directory not found: {labels_hoa_dir}", file=sys.stderr)
        sys.exit(1)

    if not t1_src_dir.is_dir():
        print(f"ERROR: T1 source directory not found: {t1_src_dir}", file=sys.stderr)
        sys.exit(1)

    t1_dst_dir.mkdir(exist_ok=True)

    ts_files = sorted(labels_hoa_dir.glob("*.nii.gz"))
    subjects = []
    for f in ts_files:
        m = re.search(r'(\d{6})', f.name)
        if m:
            subjects.append(m.group(1))

    if not subjects:
        print(f"No subjects found in {labels_hoa_dir}")
        sys.exit(0)

    print(f"Found {len(subjects)} subject(s) in {labels_hoa_dir}")
    print(f"T1 source → {t1_src_dir}")
    print(f"Output    → {t1_dst_dir}\n")

    for sid in subjects:
        run_stage_6(sid, labels_hoa_dir, t1_src_dir, t1_dst_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
