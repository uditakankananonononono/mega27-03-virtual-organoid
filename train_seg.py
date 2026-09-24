import json, sys, time
import numpy as np, torch, torch.nn.functional as F
from PIL import Image
torch.set_num_threads(1)
sys.path.insert(0, "src")
from vorganoid.seg import UNet, list_pairs, load_pair, instances_from_probs, average_precision
S = int(sys.argv[2]) if len(sys.argv) > 2 else 256; EPOCHS = int(sys.argv[1]) if len(sys.argv) > 1 else 40; CROP = 256; TAG = sys.argv[3] if len(sys.argv) > 3 else ""
tr = [load_pair(i, m, S)[:2] for i, m in list_pairs("train") + list_pairs("val")]
X = torch.tensor(np.stack([t[0] for t in tr]))[:, None]; Y = torch.tensor(np.stack([t[1].astype(np.uint8) for t in tr])); del tr
print("train imgs", len(X), flush=True)
torch.manual_seed(0); net = UNet(); opt = torch.optim.Adam(net.parameters(), 2e-3)
w = torch.tensor([1.0, 1.0, 3.0]); t0 = time.time()
for ep in range(EPOCHS):
    net.train(); perm = torch.randperm(len(X)); tot = 0
    for b in range(0, len(X), 2):
        i = perm[b:b + 2]; x, y = X[i], Y[i]
        if S > CROP:
            oy, ox = np.random.randint(0, S - CROP + 1, 2); x, y = x[..., oy:oy + CROP, ox:ox + CROP], y[..., oy:oy + CROP, ox:ox + CROP]
        if np.random.rand() < 0.5: x, y = x.flip(-1), y.flip(-1)
        if np.random.rand() < 0.5: x, y = x.flip(-2), y.flip(-2)
        k = np.random.randint(4); x, y = torch.rot90(x, k, (-2, -1)), torch.rot90(y, k, (-2, -1))
        opt.zero_grad(); loss = F.cross_entropy(net(x), y.long(), weight=w); loss.backward(); opt.step(); tot += loss.item()
    print(f"ep {ep} loss {tot:.3f} t {time.time()-t0:.0f}s", flush=True)
torch.save(net.state_dict(), f"results/unet_organoid{TAG}.pt")
net.eval(); aps = []
for img, msk in list_pairs("eval"):
    x, _, lab = load_pair(img, msk, S)
    with torch.no_grad():
        p = torch.softmax(net(torch.tensor(x)[None, None]), 1)[0]
    p = F.interpolate(p[None], size=lab.shape, mode="bilinear")[0].numpy()
    pr = instances_from_probs(p, min_size=int(20 * (lab.shape[0] / S) ** 2))
    aps.append({"image": img.split("/")[-1], "ap50": average_precision(lab, pr, 0.5),
                "ap75": average_precision(lab, pr, 0.75), "n_gt": int(len(np.unique(lab)) - 1), "n_pred": int(len(np.unique(pr)) - 1)})
    print(aps[-1], flush=True)
a = np.array([r["ap50"] for r in aps])
json.dump({"per_image": aps, "mAP50": float(a.mean()), "sd": float(a.std()), "published_mAP50": 0.76, "published_sd": 0.12,
           "epochs": EPOCHS, "input_size": S}, open(f"results/seg_eval{TAG}.json", "w"), indent=1)
print("mAP50", a.mean(), a.std())
