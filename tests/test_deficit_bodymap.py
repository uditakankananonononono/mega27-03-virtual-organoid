import sys
sys.path.insert(0, "src")
from deficit_bodymap import parse_cell

def test_parse_cell_median():
    assert parse_cell("1,2,3,4,5") == 3.0
    assert parse_cell("7") == 7.0
