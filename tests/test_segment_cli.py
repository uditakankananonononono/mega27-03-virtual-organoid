"""The delivered CLI must accept the actual train-only checkpoint and frozen selection."""
import json
from pathlib import Path
from PIL import Image
import pytest
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


def test_frozen_selection_rejects_wrong_inference_protocol(tmp_path):
    source=sorted(Path('data/raw/organoid_basic/train').glob('*_img.jpg'))[0]
    base=['segment',str(source),'--model','results/trainonly_seg/best.pt',
          '--selection','results/trainonly_seg/val_selection.json','--out',str(tmp_path/'labels.png')]
    with pytest.raises(SystemExit, match='--size 512 --tta-three'):
        main(base+['--size','256','--tta-three'])
    with pytest.raises(SystemExit, match='--size 512 --tta-three'):
        main(base+['--size','512'])
    assert not (tmp_path/'labels.png').exists()
