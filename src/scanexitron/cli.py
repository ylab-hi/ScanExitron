import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional

import typer

from . import __version__
from .core import (
    MAPQ_filter,
    external_tool_checking,
    junction_caller,
    junction_overlap_CDS_to_position_BED,
    percent_spliced_out,
)

_HEADER = (
    "chrom\tstart\tend\tname\tao\tstrand\tgene_symbol\tlength\t"
    "splice_site\tgene_id\tpso\tpsi\tdp\ttotal_junctions\n"
)

app = typer.Typer(
    help="ScanExitron: detecting exitron splicing events from RNA-Seq data",
    epilog="Requires regtools, samtools, and bedtools on PATH.",
)
vcf_app = typer.Typer(help="Convert an exitron results table to VCF format")


@app.command()
def run(
    input: Path = typer.Argument(..., help="Input BAM/CRAM file (index must be present alongside it)"),
    fasta: Path = typer.Option(..., "--fasta", "-f", help="Reference FASTA file"),
    gtf: Path = typer.Option(..., "--gtf", "-g", help="Annotation GTF file"),
    ao: int = typer.Option(3, "--ao", "-a", help="Minimum reads supporting the exitron"),
    pso: float = typer.Option(0.05, "--pso", "-p", help="Minimum PSO value"),
    mapq: int = typer.Option(50, "--mapq", "-m", help="Minimum mapping quality"),
    strand: int = typer.Option(
        1, "--strand", "-s",
        help="RNA library strand specificity: 0=unstranded, 1=first-strand/RF, 2=second-strand/FR",
    ),
    threads: int = typer.Option(1, "--threads", "-t", help="Threads for samtools"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output prefix (default: input stem)"),
    verbose: bool = typer.Option(False, "--verbose", help="Enable verbose/debug logging"),
) -> None:
    """Detect exitron splicing events from RNA-Seq data."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    external_tool_checking()

    out_bam = MAPQ_filter(in_bam=input, threads=threads, mapq=mapq)
    if not out_bam:
        logging.error("MAPQ filtering failed. Exiting.")
        raise typer.Exit(code=1)

    prefix = output if output else input.stem
    outfile = prefix + ".exitron"

    with TemporaryDirectory() as tmp_dir, open(outfile, "w") as outstream:
        outstream.write(_HEADER)
        janno_file = junction_caller(
            bam_file=out_bam,
            fasta=fasta,
            gtf=gtf,
            strand=strand,
            tmp_dir=tmp_dir,
        )
        if janno_file:
            src_file, pos_bed = junction_overlap_CDS_to_position_BED(
                janno_file,
                fasta=fasta,
                gtf=gtf,
                ao_cutoff=ao,
                tmp_dir=tmp_dir,
            )
            if src_file and pos_bed:
                percent_spliced_out(
                    bam_file=input,
                    src_exitron_file=src_file,
                    position_bed_file=pos_bed,
                    ao_cutoff=ao,
                    pso_cutoff=pso,
                    mapq=mapq,
                    out=outstream,
                )

    logging.info("Results written to %s", outfile)


@vcf_app.command()
def convert(
    input: Path = typer.Argument(..., help="Input exitron results table (.exitron file)"),
    fasta: Path = typer.Option(..., "--fasta", "-f", help="Reference FASTA file"),
    output: str = typer.Option("output.vcf", "--output", "-o", help="Output VCF file"),
) -> None:
    """Convert an exitron results table to VCF format."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    from .vcf import exitron2vcf
    exitron2vcf(in_file=input, out_vcf=output, fasta=fasta)
    logging.info("VCF written to %s", output)


def main() -> None:
    app()


def vcf_main() -> None:
    vcf_app()
