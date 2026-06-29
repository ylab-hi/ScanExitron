"""Unit tests for pure, dependency-free functions."""
from scanexitron.constants import (
    CANONICAL_SPLICE_SITES,
    CHRMS,
    CHRMS_DICT,
    NON_MITO_CHRMS,
    REVERSE_CHRMS_DICT,
)
from scanexitron.core import BED_handler


# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------

def test_chrms_count():
    assert len(CHRMS) == 25  # chr1-22 + chrX + chrY + chrM


def test_chrms_contains_expected():
    for chrom in ("chr1", "chr22", "chrX", "chrY", "chrM"):
        assert chrom in CHRMS


def test_non_mito_excludes_chrm():
    assert "chrM" not in NON_MITO_CHRMS
    assert "chr1" in NON_MITO_CHRMS
    assert len(NON_MITO_CHRMS) == 24


def test_chrms_dict_b37_to_hg():
    assert CHRMS_DICT["1"] == "chr1"
    assert CHRMS_DICT["22"] == "chr22"
    assert CHRMS_DICT["X"] == "chrX"
    assert CHRMS_DICT["Y"] == "chrY"
    assert CHRMS_DICT["MT"] == "chrM"


def test_reverse_chrms_dict_roundtrip():
    for b37, hg in CHRMS_DICT.items():
        assert REVERSE_CHRMS_DICT[hg] == b37


def test_canonical_splice_sites():
    assert "GT-AG" in CANONICAL_SPLICE_SITES
    assert "GC-AG" in CANONICAL_SPLICE_SITES
    assert "AT-AC" in CANONICAL_SPLICE_SITES
    assert "AT-AG" not in CANONICAL_SPLICE_SITES


# ---------------------------------------------------------------------------
# BED_handler
# ---------------------------------------------------------------------------

def test_bed_handler_keeps_hg_chroms(tmp_path):
    bed = tmp_path / "input.bed"
    bed.write_text(
        "chr1\t100\t200\tjunc1\t5\t+\tGT-AG\n"
        "chr22\t300\t400\tjunc2\t3\t-\tGT-AG\n"
        "chrM\t500\t600\tjunc3\t1\t+\tGT-AG\n"  # mitochondrial — kept (CHRMS member)
    )
    result = BED_handler(bed, tmp_path)
    lines = result.read_text().splitlines()
    assert len(lines) == 3


def test_bed_handler_converts_b37(tmp_path):
    bed = tmp_path / "input.bed"
    bed.write_text("1\t100\t200\tjunc1\t5\t+\tGT-AG\n")
    BED_handler(bed, tmp_path)
    lines = bed.read_text().splitlines()
    assert lines[0].startswith("chr1\t")


def test_bed_handler_drops_unknown_chroms(tmp_path):
    bed = tmp_path / "input.bed"
    bed.write_text(
        "chr1\t100\t200\tjunc1\t5\t+\tGT-AG\n"
        "scaffold_1\t100\t200\tjunc2\t5\t+\tGT-AG\n"
    )
    BED_handler(bed, tmp_path)
    lines = bed.read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("chr1\t")


# ---------------------------------------------------------------------------
# version
# ---------------------------------------------------------------------------

def test_version_string():
    from scanexitron import __version__
    parts = __version__.split(".")
    assert len(parts) == 3
    assert all(p.isdigit() for p in parts)
