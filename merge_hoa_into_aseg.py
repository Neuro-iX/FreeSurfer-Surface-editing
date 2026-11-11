#!/usr/bin/env python3
"""
merge_hoa_into_aseg.py

Replace FreeSurfer's automated subcortical segmentation with a manual HOA
labelmap inside a FreeSurfer aseg (or aparc+aseg) volume, using a
remap -> carve -> paint -> refill strategy (a true replace, not a union).

Resolution-agnostic: the HOA is resliced onto whatever grid the target is on.

Run on the two files that drive the surface stream (aseg.presurf.mgz,
aseg.auto_noCCseg.mgz), or post-hoc on a finished aparc+aseg.mgz to make the
volume match HOA exactly.

REFILL MODES (--refill), for voxels FS labelled as a replaced structure but HOA
does not (HOA smaller than FS):
  surround : fill with the nearest SURROUNDING (non-replaced) label -- WM, cortex,
             CSF, etc. The nucleus is NOT regrown, so the result matches HOA
             tightly (Dice ~1.0). DEFAULT.
  nearest  : fill with the nearest surviving label INCLUDING the nucleus itself,
             which regrows the structure back toward FS's larger size. (old behavior)
  none     : leave those voxels as 0 (unlabeled).

--split-vdc : keep HOA's VentralDC anterior/posterior distinction instead of
              collapsing both to FS 28/60. Anterior stays 28 (L)/60 (R);
              posterior becomes 8028 (L)/8060 (R). Add these to your LUT
              (see hoa_fs_custom_LUT.txt).

--include-cerebellum : also replace FS cerebellum with HOA cerebellum
              (29->8, 30->47, 31->7, 32->46). Default keeps FS cerebellum.
              (Cerebellum is not part of the cerebral ribbon; affects volumes only.)

Deps: numpy, scipy, nibabel.
"""
import argparse, sys, csv
import numpy as np
import nibabel as nib
from scipy.ndimage import distance_transform_edt

CEREBELLUM_MAP = {29: 8, 30: 47, 31: 7, 32: 46}
VDC_SPLIT_MAP  = {27: 8028, 28: 8060}   # HOA VDC-posterior L/R -> custom IDs


def load_mapping(tsv_path):
    mapping = {}
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            hv, fid = row["hoa_value"].strip(), row["fs_id"].strip()
            if hv in ("", "?") or fid in ("", "?"):
                sys.exit(f"ERROR: unresolved mapping row (hoa_value={hv!r}, fs_id={fid!r}).")
            mapping[int(hv)] = int(fid)
    if not mapping:
        sys.exit("ERROR: no mapping rows read from TSV.")
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hoa", required=True)
    ap.add_argument("--aseg", required=True)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--brainmask", default=None, help="confines the refill (recommended).")
    ap.add_argument("--reslice", action="store_true",
                    help="reslice HOA onto the aseg grid (nearest). = mri_vol2vol --regheader --nearest.")
    ap.add_argument("--refill", choices=["surround", "nearest", "none"], default="surround",
                    help="how to fill carved-but-unpainted voxels. Default: surround (tight match to HOA).")
    ap.add_argument("--split-vdc", action="store_true",
                    help="keep HOA VentralDC anterior/posterior (post -> 8028/8060).")
    ap.add_argument("--include-cerebellum", action="store_true",
                    help="replace FS cerebellum with HOA cerebellum (29->8,30->47,31->7,32->46).")
    args = ap.parse_args()

    mapping = load_mapping(args.tsv)
    if args.include_cerebellum:
        mapping.update(CEREBELLUM_MAP)
    if args.split_vdc:
        mapping.update(VDC_SPLIT_MAP)
    mapped_vals = sorted(mapping.keys())
    carve_ids = sorted(set(mapping.values()))

    aseg_img = nib.load(args.aseg)
    aseg = np.asarray(aseg_img.dataobj).astype(np.int32)
    hoa_img = nib.load(args.hoa)
    if args.reslice:
        from nibabel.processing import resample_from_to
        hoa_img = resample_from_to(hoa_img, (aseg_img.shape, aseg_img.affine), order=0)
    hoa = np.asarray(hoa_img.dataobj).astype(np.int32)
    if hoa.shape != aseg.shape:
        sys.exit(f"ERROR: HOA shape {hoa.shape} != aseg {aseg.shape}. Use --reslice or reslice first.")

    present = set(int(v) for v in np.unique(hoa) if v != 0)
    unmapped = sorted(present - set(mapped_vals))
    if unmapped:
        print(f"[warn] HOA values present but NOT in mapping (kept as FreeSurfer): {unmapped}")
    missing = sorted(set(mapped_vals) - present)
    if missing:
        print(f"[warn] mapping has HOA values not found in this labelmap: {missing}")
    print(f"[info] refill={args.refill} | split_vdc={args.split_vdc} | "
          f"cerebellum={'HOA' if args.include_cerebellum else 'FS'}")

    out = aseg.copy()
    carve_mask = np.isin(out, carve_ids)
    n_carved = int(carve_mask.sum())
    out[carve_mask] = 0

    painted = np.zeros_like(out, dtype=bool)
    for hv in mapped_vals:
        m = hoa == hv
        if m.any():
            out[m] = mapping[hv]
            painted |= m
    n_painted = int(painted.sum())

    hole = carve_mask & ~painted
    if args.brainmask is not None:
        hole &= (np.asarray(nib.load(args.brainmask).dataobj) > 0)
    n_holes = int(hole.sum())
    n_filled = 0
    if n_holes and args.refill != "none":
        if args.refill == "surround":
            src = (out != 0) & ~np.isin(out, carve_ids)   # surrounding tissue only
        else:  # nearest
            src = out != 0                                  # includes the nucleus (regrows)
        _, (ix, iy, iz) = distance_transform_edt(~src, return_indices=True)
        fill = out[ix, iy, iz]
        out[hole] = fill[hole]
        n_filled = n_holes

    nib.save(nib.MGHImage(out.astype(np.int32), aseg_img.affine, aseg_img.header), args.out)
    print(f"[ok] {args.aseg}")
    print(f"     resolution (mm)               : {tuple(round(float(z),3) for z in aseg_img.header.get_zooms()[:3])}")
    print(f"     carved / painted / refilled   : {n_carved} / {n_painted} / {n_filled}")
    if n_holes and args.refill == "none":
        print(f"     [warn] {n_holes} carved-but-unpainted voxels left as 0.")
    print(f"     wrote -> {args.out}")


if __name__ == "__main__":
    main()
