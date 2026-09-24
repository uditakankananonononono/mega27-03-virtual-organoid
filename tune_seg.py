"""Tune instance post-processing on the VAL split only, then score EVAL once.
Usage: python tune_seg.py results/unet_organoid_512.pt 512 _512
Grid: foreground prob threshold, seed interior threshold, boundary suppression, min size; optional flip TTA.
"""
import itertools, json, sys
import numpy as np, torch, torch.nn.functional as F
from scipy import ndimage as ndi
from skimage.segmentation import watershed
torch.set_num_threads(1)
sys.path.insert(0, "src")
from vorganoid.seg import UNet, list_pairs, load_pair, average_precision

ckpt, S, TAG = sys.argv[1], int(sys.argv[2]), sys.argv[3]
net = UNet(); net.load_state_dict(torch.load(ckpt)); net.eval()


def probs(img, msk, tta=True):
    x, _, lab = load_pair(img, msk, S)
    xt = torch.tensor(x)[None, None]
    with torch.no_grad():
        p = torch.softmax(net(xt), 1)
        if tta:
            p = p + torch.softmax(net(xt.flip(-1)), 1).flip(-1) + torch.softmax(net(xt.flip(-2)), 1).flip(-2)
            p = p / 3
    return p[0].numpy(), lab


def instances(p, fg_t, seed_t, bnd_t, min_size):
    fg = (p[1] + p[2]) > fg_t
    seeds = (p[1] > seed_t) & (p[2] < bnd_t)
    markers, _ = ndi.label(seeds)
    lab = watershed(-p[1], markers, mask=fg)
    sizes = np.bincount(lab.ravel()); small = np.where(sizes < min_size)[0]
    lab[np.isin(lab, small[small > 0])] = 0
    return lab


def score(pairs_probs, params):
    aps = []
    for p, lab in pairs_probs:
        # evaluate at model resolution upsampled to label size
        pu = F.interpolate(torch.tensor(p)[None], size=lab.shape, mode="bilinear")[0].numpy()
        scale = (lab.shape[0] / S) ** 2
        aps.append(average_precision(lab, instances(pu, *params[:3], int(params[3] * scale)), 0.5))
    return float(np.mean(aps)), float(np.std(aps)), aps


val = [probs(i, m) for i, m in list_pairs("val")[:20]]
grid = list(itertools.product([0.4, 0.5, 0.6], [0.5, 0.6, 0.7, 0.8], [0.2, 0.3, 0.5], [10, 20, 40]))
best = None
for g in grid:
    s = score(val, g)[0]
    if best is None or s > best[0]:
        best = (s, g)
print("best val", best, flush=True)
ev = [probs(i, m) for i, m in list_pairs("eval")]
mu, sd, aps = score(ev, best[1])
json.dump({"val_best_mAP50": best[0], "params": best[1], "eval_mAP50": mu, "eval_sd": sd, "per_image": aps,
           "published_mAP50": 0.76, "published_sd": 0.12, "tta": True}, open(f"results/seg_eval_tuned{TAG}.json", "w"), indent=1)
print("EVAL mAP50", mu, sd)
