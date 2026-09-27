import json
from pathlib import Path

def test_frozen_candidate_group_split():
    x=json.loads(Path('results/group_disjoint_seg_split.json').read_text())
    g=json.loads(Path('results/scene_graph_forensics.json').read_text())
    assert x['source_manifest_sha256']==g['image_manifest_sha256']
    assert x['threshold']==0.95 and x['component_count']==215
    assert x['split_counts']=={'train':180,'val':37,'eval':14}
    assert x['n_moved']==len(x['moved'])==5
    assert sum(x['split_counts'].values())==231
    root=Path('data/derived/group_disjoint_seg')
    for split,n in x['split_counts'].items():
        ims=list((root/split).glob('*_img.jpg'))
        assert len(ims)==n
        assert all(i.is_symlink() and i.resolve().is_file() for i in ims)
        assert all(i.with_name(i.name.replace('_img.jpg','_masks_organoid.png')).resolve().is_file() for i in ims)
    assigned={}
    for split in ('train','val','eval'):
        for i in (root/split).glob('*_img.jpg'):
            assert i.name not in assigned
            assigned[i.name]=split
    for c in g['threshold_sensitivity'][0]['multi_frame_components']:
        assert len({assigned[k.split('/')[-1]] for k in c['members']})==1
