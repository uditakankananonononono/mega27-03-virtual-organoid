"""The delivered CLI must accept the actual train-only checkpoint and frozen selection."""
import json
from pathlib import Path
from PIL import Image
from vorganoid.cli import main


def test_delivered_model_and_selection(tmp_path, capsys):
    source=sorted(Path('data/raw/organoid_basic/train').glob('*_img.jpg'))[0]
    target=tmp_path/'labels.png'
    assert main(['segment',str(source),'--model','results/trainonly_seg/best.pt',
                 '--size','512','--selection','results/trainonly_seg/val_selection.json',
                 '--tta-three','--out',str(target)])==0
    result=json.loads(capsys.readouterr().out)
    assert target.exists() and Image.open(target).size==Image.open(source).size
    assert result['n_organoids']>0
