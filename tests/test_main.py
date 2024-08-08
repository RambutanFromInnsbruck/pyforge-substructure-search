import pytest
from src.main import substructure_search


@pytest.mark.parametrize(
    "substructure, expected",
    ("c1ccccc1", [{"molecule_id": 2, "molecule_structure": "c1ccccc1"},
                  {"molecule_id": 4, "molecule_structure": "CC(=O)Oc1ccccc1C(=O)O"}]))
def test_substructure_search(substructure, expected):
    result = substructure_search(substructure)
    assert result == expected
