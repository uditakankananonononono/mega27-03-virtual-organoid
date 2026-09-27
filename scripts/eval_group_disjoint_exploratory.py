"""One exploratory score on the frozen candidate-group-disjoint reassignment.

The default extraction parameters are those in src/vorganoid/seg.py, fixed here
before the first reassigned-eval score. Previously seen images are not a benchmark.
"""
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from scipy import ndimage as ndi
from skimage.segmentation import watershed

sys.path.insert(0, 'src')
from vorganoid.seg import UNet, average_precision, load_pair

ROOT = Path('results/group_disjoint_seg')
SHA = 'cc9f8ddfd07155854c8955ccd3089de777bd53fcf6131d0ef7aadb26142360f4'
MANIFEST_SHA = 'ba003cd9645067eb455f5a2d0166a53a0a41aa73378f938e848f2ca9c73a9b65'
AUTHOR_SHA = 'aeaec17c42fefde61a8b2215a616d419e8126e9f'
PARAMS = (0.5, 0.6, 0.3, 20)  # code-default, frozen before the reassigned eval
SIZE = 512


def instances(p, fg_t, seed_t, bnd_t, min_size):
    fg = (p[1] + p[2]) > fg_t
    seeds = (p[1] > seed_t) & (p[2] < bnd_t)
    markers, _ = ndi.label(seeds)
    lab = watershed(-p[1], markers, mask=fg)
    sizes = np.bincount(lab.ravel())
    small = np.where(sizes < min_size)[0]
    lab[np.isin(lab, small[small > 0])] = 0
    return lab


def main():
    assert len(sys.argv) == 2, 'pass path to pinned OrgaSegment author checkout'
    assert not (ROOT / 'exploratory_eval.json').exists(), 'eval already scored; no rerun'
    author_path = Path(sys.argv[1])
    assert subprocess.check_output(['git', '-C', str(author_path), 'rev-parse', 'HEAD'], text=True).strip() == AUTHOR_SHA
    spec = importlib.util.spec_from_file_location('author_metrics', author_path / 'lib' / 'metrics.py')
    author = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(author)
    manifest = json.loads(Path('results/group_disjoint_seg_split.json').read_text())
    assert manifest['source_manifest_sha256'] == MANIFEST_SHA
    assert manifest['threshold'] == .95 and manifest['split_counts'] == {'train': 180, 'val': 37, 'eval': 14}
    graph = json.loads(Path('results/scene_graph_forensics.json').read_text())
    assert graph['image_manifest_sha256'] == MANIFEST_SHA
    assert graph['threshold_sensitivity'][0]['threshold'] == .95
    pairs = {}
    for split, n in manifest['split_counts'].items():
        images = sorted((Path('data/derived/group_disjoint_seg') / split).glob('*_img.jpg'))
        assert len(images) == n
        pairs[split] = [(im, im.with_name(im.name.replace('_img.jpg', '_masks_organoid.png'))) for im in images]
        assert all(mask.exists() for _, mask in pairs[split])
    assignments = {im.name: split for split, rows in pairs.items() for im, _ in rows}
    assert len(assignments) == 231
    for comp in graph['threshold_sensitivity'][0]['multi_frame_components']:
        assert len({assignments[name.split('/')[-1]] for name in comp['members']}) == 1
    assert hashlib.sha256((ROOT / 'best.pt').read_bytes()).hexdigest() == SHA
    history = json.loads((ROOT / 'history.json').read_text())
    assert len(history) == 60 and history[-1]['epoch'] == 60
    assert history[-1]['val_loss'] == min(z['val_loss'] for z in history if z['val_loss'] is not None)
    ck = torch.load(ROOT / 'best.pt', map_location='cpu', weights_only=False)
    assert ck['epoch'] == 59 and ck['source_manifest_sha256'] == MANIFEST_SHA
    assert abs(ck['val_loss'] - history[-1]['val_loss']) < 1e-14
    torch.set_num_threads(1)
    net = UNet()
    net.load_state_dict(ck['net'])
    net.eval()
    out = []
    for img, msk in pairs['eval']:
        x, _, lab = load_pair(str(img), str(msk), SIZE)
        xt = torch.tensor(x)[None, None]
        with torch.no_grad():
            p = (torch.softmax(net(xt), 1) + torch.softmax(net(xt.flip(-1)), 1).flip(-1) + torch.softmax(net(xt.flip(-2)), 1).flip(-2)) / 3
        pu = F.interpolate(p, size=lab.shape, mode='bilinear')[0].numpy()
        pr = instances(pu, *PARAMS[:3], int(PARAMS[3] * (lab.shape[0] / SIZE) ** 2))
        ids = np.unique(pr)
        pr = np.searchsorted(ids, pr).astype(np.int32)
        ap_author, tp, fp, fn = author.average_precision(lab.astype(np.int32), pr, .5)
        ap_ours = average_precision(lab, pr, .5)
        row = {'image': img.name, 'ap_author': float(ap_author[0]), 'ap_ours': float(ap_ours),
               'tp_author': int(tp[0]), 'fp_author': int(fp[0]), 'fn_author': int(fn[0]),
               'n_gt': int(len(np.unique(lab)) - 1), 'n_pred': int(len(np.unique(pr)) - 1),
               'gt_id_gaps': int(lab.max() - (len(np.unique(lab)) - 1))}
        print(row, flush=True)
        out.append(row)
    a = np.array([r['ap_author'] for r in out])
    rng = np.random.default_rng(20260927)
    draws = rng.choice(a, size=(20000, len(a)), replace=True).mean(axis=1)
    result = {'status': 'EXPLORATORY; reassigned eval on known source images, not independent or prospective benchmark',
              'protocol': 'notes/prereg_group_disjoint_seg.md', 'source_manifest_sha256': MANIFEST_SHA,
              'model': str(ROOT / 'best.pt'), 'model_sha256': SHA, 'selected_epoch': 60,
              'selected_val_loss': ck['val_loss'], 'params': PARAMS,
              'selection': 'default in src/vorganoid/seg.py frozen in this script before new eval; no reassigned-val postprocessing grid',
              'author_release': 'https://github.com/Living-Technologies/OrgaSegment/tree/v1.0.1',
              'author_release_commit': AUTHOR_SHA, 'per_image': out,
              'mean_author': float(a.mean()), 'mean_ours': float(np.mean([r['ap_ours'] for r in out])),
              'image_bootstrap_ci95': list(map(float, np.quantile(draws, [.025, .975]))),
              'bootstrap_seed': 20260927, 'bootstrap_draws': 20000,
              'limitations': '14 exploratory eval files: 12 previously scored original eval files plus one original-train and one original-val file moved into eval; all 14 were known source images, and the moved files have original-eval scene-proxy links. Candidate groups do not prove donor/field independence; published author-model paired predictions unavailable; no benchmark-win claim'}
    (ROOT / 'exploratory_eval.json').write_text(json.dumps(result, indent=2) + '\n')
    print('mean author', result['mean_author'], 'ours', result['mean_ours'], 'CI', result['image_bootstrap_ci95'], flush=True)

if __name__ == '__main__':
    main()
