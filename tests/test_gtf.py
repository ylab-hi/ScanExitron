from pathlib import Path
import pytest
from scanexitron.gtf import extract_cds_bed


def test_extract_cds_bed_success(tmp_path):
    # Create a dummy GTF file
    gtf_file = tmp_path / "test_annotation.gtf"
    gtf_content = (
        "chr17\tYangLab\tCDS\t100\t200\t.\t+\t0\tgene_id \"ENSG00000005884\"; gene_name \"ITGA3\";\n"
        "chr17\tYangLab\tCDS\t300\t400\t.\t-\t0\tgene_id \"ENSG00000005379\";\n"
        "chr17\tYangLab\texon\t50\t150\t.\t+\t.\tgene_id \"ENSG00000005884\";\n"
    )
    gtf_file.write_text(gtf_content)

    # Run extraction
    bed_path = extract_cds_bed(gtf_file)

    # Verify bed file path is correctly formatted next to the GTF file
    assert bed_path == tmp_path / "test_annotation.CDS.bed"
    assert bed_path.exists()

    # Verify content
    lines = bed_path.read_text().splitlines()
    assert len(lines) == 2

    # HTSeq handles GFF 1-based start coordinate and converts it to 0-based half-open (start - 1)
    # GFF 100 -> 0-based 99; 200 -> 199
    # GFF 300 -> 0-based 299; 400 -> 399
    assert lines[0] == "chr17\t99\t199\tENSG00000005884\tITGA3\t+"
    # For ENSG00000005379, it falls back to gene_id as gene_name
    assert lines[1] == "chr17\t299\t399\tENSG00000005379\tENSG00000005379\t-"


def test_extract_cds_bed_caching(tmp_path, caplog):
    import logging
    # Create a dummy GTF file
    gtf_file = tmp_path / "test_cache.gtf"
    gtf_content = "chr17\tYangLab\tCDS\t100\t200\t.\t+\t0\tgene_id \"ENSG00000005884\"; gene_name \"ITGA3\";\n"
    gtf_file.write_text(gtf_content)

    # First run (extracts and caches)
    extract_cds_bed(gtf_file)

    # Second run (should use cache)
    with caplog.at_level(logging.INFO):
        extract_cds_bed(gtf_file)
        assert any("Using cached CDS BED" in record.message for record in caplog.records)

    # Third run with force=True (should regenerate)
    caplog.clear()
    with caplog.at_level(logging.INFO):
        extract_cds_bed(gtf_file, force=True)
        assert any("Extracting CDS features" in record.message for record in caplog.records)


def test_extract_cds_bed_missing_file():
    with pytest.raises(FileNotFoundError):
        extract_cds_bed("non_existent_file.gtf")
