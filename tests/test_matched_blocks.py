"""Strict CLI result must reproduce the registered negative, not pooled positives."""
import json
from pathlib import Path
import pandas as pd
import pytest
from vorganoid.matched_blocks import matched_block_effects
from vorganoid.cli import main

DATA = Path('data/raw/orgasegment/dis_merged_A0.csv')


def test_exact_registered_real_data(capsys):
    assert main(['matched-blocks', str(DATA), '--drug', 'VX445_VX661_VX770']) == 0
    out = json.loads(capsys.readouterr().out)
    frozen = json.load(open('results/blocked_size.json'))
    assert (out['n_donors'], out['n_positive_donors'], out['verdict']) == (12, 9, 'FAIL')
    assert out['one_sided_sign_p'] == frozen['sign_p_one_sided']
    assert len(out['blocks']) == frozen['n_eligible_blocks'] == 31
    expected = {z['donor']: z['median_effect'] for z in frozen['per_donor'] if z['n_blocks']}
    assert {z['donor'] for z in out['donors']} == set(expected)
    for z in out['donors']:
        assert z['median_effect'] == pytest.approx(expected[z['donor']], abs=1e-12)


def test_missing_or_invalid_inputs():
    with pytest.raises(ValueError, match='missing columns'):
        matched_block_effects(pd.DataFrame({'A0': [2]}), 'drug')
    with pytest.raises(ValueError, match='small_max'):
        matched_block_effects(pd.read_csv(DATA).head(), 'drug', small_max=2000, large_min=1000)


def test_reject_same_arm_and_nonfinite_inputs():
    base = pd.read_csv(DATA).head(20)
    with pytest.raises(ValueError, match="must differ"):
        matched_block_effects(base, "DMSO", "DMSO")
    with pytest.raises(ValueError, match="finite"):
        matched_block_effects(base, "drug", small_max=float("nan"))
    with pytest.raises(ValueError, match="numeric"):
        matched_block_effects(base.assign(A0="bad"), "drug")
    altered = base.copy()
    altered.loc[altered.index[0], "swelling"] = float("inf")
    assert matched_block_effects(altered, "VX445_VX661_VX770")["n_donors"] == 0
