import json
from pathlib import Path
x=Path(__file__).resolve().parents[1]/'results/postresult_mouse_yolo_domain.json'
def test_mouse_domain_all_images_and_negative():
    r=json.loads(x.read_text());z=r['per_image'];s=r['summary']
    assert len(z)==r['n_expected_val']==84
    assert len({i['image'] for i in z})==84
    assert s['tp']==sum(i['tp'] for i in z)==850
    assert s['fp']==sum(i['fp'] for i in z)==911
    assert s['fn']==sum(i['fn'] for i in z)==1618
    assert sum(i['n_truth'] for i in z)==s['tp']+s['fn']==2468
    assert s['micro_f1']==2*s['tp']/(2*s['tp']+s['fp']+s['fn'])
    assert s['invalid_zero_area_annotations']==1
    assert all(i['tp']+i['fn']==i['n_truth'] and i['tp']+i['fp']==i['n_pred'] for i in z)
    assert s['micro_recall']<.35
