"""CNN instance segmentation of organoids (U-Net, 3-class: background /
interior / boundary) with connected-component instance extraction, and the
COCO-style per-image Average Precision at an IoU threshold used by OrgaSegment.
"""
from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from PIL import Image
from scipy import ndimage as ndi
from skimage.segmentation import find_boundaries, watershed

DATA = Path(__file__).resolve().parents[2] / "data" / "raw" / "organoid_basic"


def list_pairs(split: str) -> list[tuple[str, str]]:
    pairs = []
    for img in sorted(glob.glob(str(DATA / split / "*_img.jpg"))):
        mask = img.replace("_img.jpg", "_masks_organoid.png")
        if Path(mask).exists():
            pairs.append((img, mask))
    return pairs


def load_pair(img_path: str, mask_path: str, size: int = 256):
    img = Image.open(img_path).convert("L")
    lab = np.array(Image.open(mask_path))
    x = np.array(img.resize((size, size), Image.BILINEAR), dtype=np.float32) / 255.0
    x = (x - x.mean()) / (x.std() + 1e-6)
    lab_s = np.array(Image.fromarray(lab).resize((size, size), Image.NEAREST))
    tgt = np.zeros(lab_s.shape, np.int64)
    tgt[lab_s > 0] = 1
    tgt[find_boundaries(lab_s, mode="inner") & (lab_s > 0)] = 2
    return x, tgt, lab


def block(ci, co):
    return nn.Sequential(nn.Conv2d(ci, co, 3, padding=1), nn.BatchNorm2d(co), nn.ReLU(inplace=True),
                         nn.Conv2d(co, co, 3, padding=1), nn.BatchNorm2d(co), nn.ReLU(inplace=True))


class UNet(nn.Module):
    def __init__(self, base: int = 16, n_cls: int = 3):
        super().__init__()
        c = [base, base * 2, base * 4, base * 8]
        self.d = nn.ModuleList([block(1, c[0]), block(c[0], c[1]), block(c[1], c[2]), block(c[2], c[3])])
        self.u = nn.ModuleList([nn.ConvTranspose2d(c[i + 1], c[i], 2, stride=2) for i in range(3)])
        self.ub = nn.ModuleList([block(c[i] * 2, c[i]) for i in range(3)])
        self.out = nn.Conv2d(c[0], n_cls, 1)

    def forward(self, x):
        skips = []
        for i, b in enumerate(self.d):
            x = b(x)
            if i < 3:
                skips.append(x); x = F.max_pool2d(x, 2)
        for i in reversed(range(3)):
            x = self.ub[i](torch.cat([self.u[i](x), skips[i]], 1))
        return self.out(x)


def instances_from_probs(prob: np.ndarray, min_size: int = 20) -> np.ndarray:
    """prob: 3 x H x W softmax. Seeds = confident interior; watershed on fg."""
    fg = (prob[1] + prob[2]) > 0.5
    seeds = (prob[1] > 0.6) & (prob[2] < 0.3)
    markers, _ = ndi.label(seeds)
    lab = watershed(-prob[1], markers, mask=fg)
    sizes = np.bincount(lab.ravel())
    small = np.where(sizes < min_size)[0]
    lab[np.isin(lab, small[small > 0])] = 0
    return lab


def iou_matrix(gt: np.ndarray, pr: np.ndarray) -> np.ndarray:
    g_ids = np.unique(gt); g_ids = g_ids[g_ids > 0]
    p_ids = np.unique(pr); p_ids = p_ids[p_ids > 0]
    if len(g_ids) == 0 or len(p_ids) == 0:
        return np.zeros((len(g_ids), len(p_ids)))
    gi = np.searchsorted(g_ids, gt[gt > 0]); pj = pr[gt > 0]
    # intersection via joint histogram
    pmap = {p: k for k, p in enumerate(p_ids)}
    inter = np.zeros((len(g_ids), len(p_ids)))
    m = pj > 0
    np.add.at(inter, (gi[m], np.array([pmap[v] for v in pj[m]], dtype=int)), 1)
    ga = np.array([(gt == g).sum() for g in g_ids])[:, None]
    pa = np.array([(pr == p).sum() for p in p_ids])[None, :]
    return inter / (ga + pa - inter)


def average_precision(gt: np.ndarray, pr: np.ndarray, thr: float = 0.5) -> float:
    """Stardist/Cellpose-style AP = TP / (TP + FP + FN) with greedy IoU matching."""
    iou = iou_matrix(gt, pr)
    ng, npr = iou.shape
    if ng == 0 and npr == 0:
        return 1.0
    tp = 0
    if ng and npr:
        used_g, used_p = set(), set()
        for g, p in sorted(zip(*np.where(iou >= thr)), key=lambda t: -iou[t]):
            if g not in used_g and p not in used_p:
                used_g.add(g); used_p.add(p); tp += 1
    return tp / (tp + (npr - tp) + (ng - tp))
