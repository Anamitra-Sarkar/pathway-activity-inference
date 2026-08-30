import tempfile
from pathlib import Path

import pytest
from data_pipeline.gmt_parser import parse_gmt, parse_gmt_string


def test_parse_gmt_string_basic():
    content = "HALLMARK_IFN_ALPHA_RESPONSE\tInterferon\tIFIT1\tMX1\tISG15\nHALLMARK_P53\tP53 pathway\tTP53\tCDKN1A\tMDM2"
    result = parse_gmt_string(content)
    assert len(result) == 2
    assert result["HALLMARK_IFN_ALPHA_RESPONSE"]["genes"] == ["IFIT1", "MX1", "ISG15"]
    assert result["HALLMARK_P53"]["description"] == "P53 pathway"


def test_parse_gmt_file_and_dedup():
    content = "PW1\tdesc\tG1\tG2\tG1\tG3\nPW2\tdesc2\tG4\tG5"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gmt", delete=False) as f:
        f.write(content)
        path = f.name
    try:
        result = parse_gmt(path)
        assert result["PW1"]["genes"] == ["G1", "G2", "G3"]  # dedup preserve order
        assert result["PW2"]["genes"] == ["G4", "G5"]
    finally:
        Path(path).unlink()


def test_parse_gmt_empty_gene_filtered():
    content = "PW1\tdesc\tG1\t\tG2"
    result = parse_gmt_string(content)
    assert result["PW1"]["genes"] == ["G1", "G2"]


def test_parse_gmt_malformed_raises():
    with pytest.raises(ValueError):
        parse_gmt_string("OnlyOneField")

def test_parse_gmt_missing_file():
    with pytest.raises(FileNotFoundError):
        parse_gmt("/tmp/nonexistent_xyz.gmt")

def test_parse_gmt_duplicate_name_raises():
    content = "PW1\tdesc\tG1\nPW1\tdesc2\tG2"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gmt", delete=False) as f:
        f.write(content)
        path = f.name
    try:
        with pytest.raises(ValueError, match="Duplicate"):
            parse_gmt(path)
    finally:
        Path(path).unlink()
