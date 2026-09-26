"""VO3-R1 candidate v2: fine-tune sealed best.pt with stronger augmentation, low LR.
Levers used (locked in notes/prereg_ap50_beat.md): stronger augmentation + longer schedule.
Same train split only; val used ONLY for model selection (CE at epochs 5,10,...,30); eval untouched.
Incremental: --steps epochs per invocation; resume via results/trainonly_seg/resume_v2.pt."""
import argparse, json, sys, time
from pathlib import Path
import numpy as np, torch, torch.nn.functional as F
sys.path.insert(0, 'src')
from vorganoid.seg import UNet, list_pairs, load_pair

p = argparse.ArgumentParser(); p.add_argument('--steps', type=int, default=1); a = p.parse_args()
assert 1 <= a.steps <= 3
S = 512; CROP = 256; MAX = 30; LR = 5e-4
root = Path('results/trainonly_seg')
CK = root / 'resume_v2.pt'; BEST = root / 'best_v2.pt'; LOG = root / 'history_v2.json'
torch.set_num_threads(1)
train = list_pairs('train'); val = list_pairs('val'); ev = list_pairs('eval')
assert (len(train), len(val), len(ev)) == (184, 35, 12)
assert not ({i for i, _ in train} & {i for i, _ in val}) and not ({i for i, _ in train} & {i for i, _ in ev})
tr = [load_pair(i, m, S)[:2] for i, m in train]; va = [load_pair(i, m, S)[:2] for i, m in val]
X = torch.tensor(np.stack([t[0] for t in tr]))[:, None]; Y = torch.tensor(np.stack([t[1].astype(np.uint8) for t in tr])); del tr
VX = torch.tensor(np.stack([t[0] for t in va]))[:, None]; VY = torch.tensor(np.stack([t[1].astype(np.uint8) for t in va])); del va

def augment(x, y):
    oy, ox = np.random.randint(0, S - CROP + 1, 2)
    x, y = x[..., oy:oy + CROP, ox:ox + CROP], y[..., oy:oy + CROP, ox:ox + CROP]
    if np.random.rand() < 0.5: x, y = x.flip(-1), y.flip(-1)
    if np.random.rand() < 0.5: x, y = x.flip(-2), y.flip(-2)
    k = int(np.random.randint(0, 4))
    if k: x, y = torch.rot90(x, k, (-2, -1)), torch.rot90(y, k, (-2, -1))
    g = float(np.random.uniform(0.9, 1.1)); b = float(np.random.uniform(-0.05, 0.05))
    return x * g + b, y

net = UNet(); opt = torch.optim.Adam(net.parameters(), LR); w = torch.tensor([1., 1., 3.]); best = float('inf'); history = []
if CK.exists():
    ck = torch.load(CK, map_location='cpu', weights_only=False)
    net.load_state_dict(ck['net']); opt.load_state_dict(ck['opt']); torch.set_rng_state(ck['torch_rng']); np.random.set_state(ck['numpy_rng'])
    start = ck['epoch'] + 1; best = ck['best']; history = ck['history']; print('RESUME v2', start, flush=True)
else:
    src = torch.load(root / 'best.pt', map_location='cpu', weights_only=False); assert src['epoch'] == 39
    net.load_state_dict(src['net'])
    torch.manual_seed(2077); np.random.seed(2077); start = 0
    print('START v2 fine-tune from best.pt (epoch 40, val CE 0.075751) lr', LR, flush=True)
for ep in range(start, min(start + a.steps, MAX)):
    t0 = time.time(); net.train(); perm = torch.randperm(len(X)); tot = 0.
    for b in range(0, len(X), 2):
        ii = perm[b:b + 2]; x, y = X[ii], Y[ii]
        x, y = augment(x, y)
        logits = net(x); loss = F.cross_entropy(logits, y.long(), weight=w)
        opt.zero_grad(); loss.backward(); opt.step(); tot += float(loss)
    vl = None
    if (ep + 1) % 5 == 0 or ep + 1 == MAX:
        net.eval()
        with torch.no_grad():
            vl = float(sum(F.cross_entropy(net(VX[i:i + 4]), VY[i:i + 4].long(), weight=w) for i in range(0, len(VX), 4)) / 9)
        if vl < best:
            best = vl; torch.save({'net': net.state_dict(), 'epoch': ep, 'val_loss': vl}, BEST)
    history.append({'epoch': ep + 1, 'train_loss_sum': tot, 'val_loss': vl, 'best_val_loss': best, 'seconds': round(time.time() - t0, 1)})
    torch.save({'net': net.state_dict(), 'opt': opt.state_dict(), 'torch_rng': torch.get_rng_state(),
                'numpy_rng': np.random.get_state(), 'epoch': ep, 'best': best, 'history': history}, CK)
    LOG.write_text(json.dumps(history))
    print('v2 epoch', ep + 1, 'train', round(tot, 3), 'val', vl, 'best', round(best, 5), flush=True)
