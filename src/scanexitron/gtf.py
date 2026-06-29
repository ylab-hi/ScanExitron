"""GTF utilities for ScanExitron.

Provides :func:`extract_cds_bed`, which reads a GTF file with HTSeq and writes
a sorted, deduplicated 6-column BED of all CDS features.  The result is cached
alongside the GTF so repeated runs avoid re-parsing the (often large) annotation.
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def extract_cds_bed(gtf_path: str | Path, force: bool = False) -> Path:
    """Extract CDS features from *gtf_path* and write a 6-column BED file.

    The BED is written to ``<gtf_stem>.CDS.bed`` in the same directory as the
    GTF so it is reused across runs.  Pass ``force=True`` to regenerate it even
    when a cached copy exists.

    Output columns (0-based, tab-separated):
        chrom, start, end, gene_id, gene_name, strand

    ``start`` is 0-based (half-open), matching BED convention.

    Args:
        gtf_path: Path to a GTF annotation file (plain or gzip-compressed).
        force:    Overwrite the cached BED even if it already exists.

    Returns:
        :class:`~pathlib.Path` pointing to the written (or cached) CDS BED.

    Raises:
        ImportError: If the ``HTSeq`` package is not installed.
        FileNotFoundError: If *gtf_path* does not exist.
    """
    try:
        import HTSeq  # noqa: PLC0415
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "HTSeq is required for GTF parsing. Install it with: pip install htseq (or conda install -c bioconda htseq)"
        ) from exc

    gtf_path = Path(gtf_path)
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    # Cache lives next to the GTF: e.g. gencode.v38.annotation.CDS.bed
    # Strip one suffix at a time so .gtf.gz → .gtf → stem
    stem = gtf_path.name
    for ext in (".gz", ".gtf", ".gff", ".gff3"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    bed_path = gtf_path.parent / f"{stem}.CDS.bed"

    if bed_path.exists() and not force:
        logger.info("Using cached CDS BED: %s", bed_path)
        return bed_path

    logger.info("Extracting CDS features from %s -> %s", gtf_path, bed_path)

    rows: set[tuple[str, int, int, str, str, str]] = set()
    reader = HTSeq.GFF_Reader(str(gtf_path), end_included=False)
    for feature in reader:
        if feature.type != "CDS":
            continue
        iv = feature.iv
        gene_id = feature.attr.get("gene_id", ".")
        # gene_name is present in GENCODE; fall back to gene_id for other sources
        gene_name = feature.attr.get("gene_name", gene_id)
        rows.add((iv.chrom, iv.start, iv.end, gene_id, gene_name, iv.strand))

    sorted_rows = sorted(rows, key=lambda r: (r[0], r[1]))

    with open(bed_path, "w") as fh:
        for chrom, start, end, gene_id, gene_name, strand in sorted_rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t{gene_name}\t{strand}\n")

    logger.info("CDS BED written: %d regions -> %s", len(sorted_rows), bed_path)
    return bed_path
