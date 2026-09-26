"""vorganoid command-line tool.

  vorganoid sizeaware TABLE --drug NAME [--control DMSO] [--bins 4] [--forskolin-col COL]
      Size-stratified CFTR-modulator effect from a per-organoid FIS table (CSV with donor, condition, A0, swelling).
  vorganoid matched-blocks TABLE --drug NAME [--small-max 722] [--large-min 1400]
      Strict donor/experiment/forskolin-dose matched size contrast; no clinical validation.
  vorganoid segment IMAGE --model unet.pt [--size 256] [--out labels.png]
      Instance-segment organoids in a brightfield image with the trained U-Net + watershed.
  vorganoid fidelity COUNTS --reference gtex_median_tpm.gct.gz [--organ-tissue "Liver"] [--top 5] [--purity hpa_single_cell.tsv.gz]
      Score an organoid gene-count matrix (Ensembl IDs, samples in columns) against reference tissues;
      reports top matches, rank of the intended tissue, and rank of GTEx cultured fibroblasts.
  vorganoid fidelity PROFILE --profile --reference gtex_median_tpm.gct.gz
      Score a precomputed, log-normalized expression profile (Ensembl ID, lcpm),
      e.g. data/geo/profiles/GSE186249.csv.gz. Do not treat this as a count matrix.
"""
from __future__ import annotations

import argparse, json, sys

import numpy as np
import pandas as pd

from .sizeaware import size_effects, attenuation


def cmd_sizeaware(a):
    df = pd.read_csv(a.table)
    if a.forskolin_col: df = df[df[a.forskolin_col] > 0]
    need = {"donor", "condition", "A0", "swelling"}
    if not need <= set(df.columns): sys.exit(f"table needs columns {sorted(need)}")
    per, edges = size_effects(df, a.drug, a.control, a.bins)
    res = attenuation(per) | {"bin_edges_A0": [float(e) for e in edges], "drug": a.drug, "control": a.control}
    print(per.round(3).to_string(index=False)); print(json.dumps(res, indent=1))
    return 0


def cmd_matched_blocks(a):
    from .matched_blocks import matched_block_effects
    try:
        result = matched_block_effects(pd.read_csv(a.table), a.drug, a.control,
                                       a.small_max, a.large_min, a.min_per_cell)
    except ValueError as exc:
        sys.exit(str(exc))
    print(json.dumps(result, indent=2))
    return 0


def cmd_segment(a):
    import torch, torch.nn.functional as F
    from PIL import Image
    from .seg import UNet, instances_from_probs
    if a.selection and (a.size != 512 or not a.tta_three):
        sys.exit("the supplied validation selection was frozen for --size 512 --tta-three")
    net = UNet()
    checkpoint = torch.load(a.model, map_location="cpu", weights_only=True)
    # Legacy weights are bare state dictionaries; train-only checkpoints include epoch and optimizer state.
    net.load_state_dict(checkpoint["net"] if "net" in checkpoint and "epoch" in checkpoint else checkpoint)
    net.eval()
    source = Image.open(a.image).convert("L")
    img = np.asarray(source, dtype=np.float32)
    x = np.asarray(source.resize((a.size, a.size), Image.BILINEAR), dtype=np.float32) / 255.0
    x = (x - x.mean()) / (x.std() + 1e-6)
    xt = torch.tensor(x)[None, None]
    with torch.no_grad():
        logits = [torch.softmax(net(xt), 1)]
        if a.tta_three:
            logits += [torch.softmax(net(xt.flip(-1)), 1).flip(-1),
                       torch.softmax(net(xt.flip(-2)), 1).flip(-2)]
        p = torch.stack(logits).mean(0)
    p = F.interpolate(p, size=img.shape, mode="bilinear")[0].numpy()
    if a.selection:
        from scipy import ndimage as ndi
        from skimage.segmentation import watershed
        selection = json.load(open(a.selection))
        if selection.get("eval_used") is not False:
            sys.exit("selection must be frozen before eval")
        if selection.get("selected_epoch") != checkpoint.get("epoch", -1) + 1:
            sys.exit("selected epoch does not match model checkpoint")
        fg_t, seed_t, bnd_t, minimum = selection["selected_params"]
        fg = (p[1] + p[2]) > fg_t
        seeds = (p[1] > seed_t) & (p[2] < bnd_t)
        markers, _ = ndi.label(seeds)
        lab = watershed(-p[1], markers, mask=fg)
        sizes = np.bincount(lab.ravel())
        minimum = int(minimum * (img.shape[0] * img.shape[1]) / (a.size * a.size))
        small = np.where(sizes < minimum)[0]
        lab[np.isin(lab, small[small > 0])] = 0
    else:
        lab = instances_from_probs(p, min_size=int(20 * (img.shape[0] * img.shape[1]) / (a.size * a.size)))
    ids = np.unique(lab)
    lab = np.searchsorted(ids, lab).astype(np.uint16)
    Image.fromarray(lab).save(a.out)
    areas = np.bincount(lab.ravel())[1:]
    print(json.dumps({"n_organoids": int((areas > 0).sum()), "median_area_px": float(np.median(areas[areas > 0])) if (areas > 0).any() else 0}))
    return 0

def cmd_fidelity(a):
    from .fidelity import mean_log_cpm, fidelity_scores, load_gtex
    sep = "," if ".csv" in a.counts else "\t"
    df = pd.read_csv(a.counts, sep=sep, index_col=0); df.index = df.index.astype(str).str.split(".").str[0]
    from .fidelity import parenchymal_genes
    genes = parenchymal_genes(a.purity) if a.purity else None
    if a.profile:
        if df.shape[1] != 1: sys.exit("--profile expects exactly one precomputed expression column")
        profile = pd.to_numeric(df.iloc[:, 0], errors="raise")
        if not np.isfinite(profile).all(): sys.exit("--profile contains non-finite values")
    else:
        try: profile = mean_log_cpm(df)
        except ValueError as exc: sys.exit(f"count input invalid: {exc}; use --profile only if input is already log-normalized")
    sc = fidelity_scores(profile, load_gtex(a.reference), genes=genes)
    rank = {t: int(list(sc.index).index(t)) + 1 for t in sc.index}
    fib = "Cells - Cultured fibroblasts"
    res = {"top": {t: round(float(v), 4) for t, v in sc.head(a.top).items()}, "fibroblast_rank": rank.get(fib),
           "culture_attractor": rank.get(fib) == 1}
    if a.organ_tissue: res["organ_tissue_rank"] = rank.get(a.organ_tissue)
    print(json.dumps(res, indent=1)); return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="vorganoid"); sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("sizeaware"); s.add_argument("table"); s.add_argument("--drug", required=True); s.add_argument("--control", default="DMSO")
    s.add_argument("--bins", type=int, default=4); s.add_argument("--forskolin-col"); s.set_defaults(f=cmd_sizeaware)
    b = sp.add_parser("matched-blocks", help="donor/plate/dose matched drug-size contrast; no clinical use")
    b.add_argument("table"); b.add_argument("--drug", required=True); b.add_argument("--control", default="DMSO")
    b.add_argument("--small-max", type=float, default=722.); b.add_argument("--large-min", type=float, default=1400.)
    b.add_argument("--min-per-cell", type=int, default=3); b.set_defaults(f=cmd_matched_blocks)
    g = sp.add_parser("segment"); g.add_argument("image"); g.add_argument("--model", required=True); g.add_argument("--size", type=int, default=256)
    g.add_argument("--out", default="labels.png"); g.add_argument("--selection", help="frozen validation-selected postprocessing JSON")
    g.add_argument("--tta-three", action="store_true", help="average original, horizontal and vertical flip predictions")
    g.set_defaults(f=cmd_segment)
    f = sp.add_parser("fidelity"); f.add_argument("counts"); f.add_argument("--reference", required=True)
    f.add_argument("--organ-tissue"); f.add_argument("--purity", metavar="HPA_TSV", help="restrict to HPA parenchymal-only genes (purity-corrected score)"); f.add_argument("--top", type=int, default=5)
    f.add_argument("--profile", action="store_true", help="input is a one-column, precomputed log expression profile, not raw counts")
    f.set_defaults(f=cmd_fidelity)
    a = ap.parse_args(argv); return a.f(a)


if __name__ == "__main__":
    sys.exit(main())
