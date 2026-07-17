import logging
import secrets
import shutil
import subprocess
import sys
from collections import OrderedDict
from io import BytesIO
from pathlib import Path
from typing import IO, Optional

import numpy as np
from pyfaidx import Fasta

from .constants import CANONICAL_SPLICE_SITES, CHRMS, CHRMS_DICT, NON_MITO_CHRMS
from .gtf import extract_cds_bed

logger = logging.getLogger(__name__)


def _remove(infile: str | Path) -> None:
    path = Path(infile)
    if path.is_file():
        path.unlink()


def run_cmd(cmd: str, msg: str = "") -> tuple[bool, Optional[bytes]]:
    logger.info(cmd)
    if "," in msg:
        begin, finish = msg.split(",", 1)
        logger.info(begin.strip())
    else:
        finish = msg
    try:
        result = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT, stdin=subprocess.PIPE
        )
    except subprocess.CalledProcessError as err:
        logger.error("Command failed: %s\n%s", err, err.output.decode(errors="replace"))
        return False, None
    logger.info(finish.strip())
    return True, result


def external_tool_checking() -> None:
    """Verify that regtools, bedtools, and samtools are on PATH."""
    for tool in ("regtools", "bedtools", "samtools"):
        try:
            path = subprocess.check_output(["which", tool], stderr=subprocess.STDOUT)
            logger.info("Found '%s': %s", tool, path.decode().strip())
        except subprocess.CalledProcessError:
            logger.error("Required tool '%s' not found. Please install it and retry.", tool)
            sys.exit(1)


def _done_file(name: str | Path) -> None:
    Path(f"{name}.done").write_text("done!")


def BED_handler(inbed: str | Path, tmp_dir: str | Path) -> Path:
    """Keep only canonical chromosomes; convert b37/GRCh notation to hg19/38."""
    inbed = Path(inbed)
    rnd_id = secrets.token_hex(16)
    tmp_file = Path(tmp_dir) / f"tmp.{rnd_id}.txt"
    with open(tmp_file, "w") as tmp, open(inbed) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in CHRMS:
                tmp.write("\t".join(parts) + "\n")
            elif parts[0] in CHRMS_DICT:
                parts[0] = CHRMS_DICT[parts[0]]
                tmp.write("\t".join(parts) + "\n")
    shutil.move(str(tmp_file), str(inbed))
    return inbed.resolve()


def MAPQ_filter(
    in_bam: str | Path, threads: int = 6, mapq: int = 50
) -> Optional[str]:
    """Filter BAM by mapping quality; returns path to filtered BAM or None on failure."""
    in_bam = Path(in_bam)
    prefix = in_bam.stem
    hq_bam = f"{prefix}.hq.bam"
    if Path(f"{hq_bam}.done").exists():
        logger.info("%s found, skipping MAPQ filtering.", hq_bam)
        return hq_bam
    cmd = (
        f"samtools view -q {mapq} -@ {threads} -O BAM -o {hq_bam} {in_bam} "
        f"&& samtools index {hq_bam}"
    )
    flag, _ = run_cmd(cmd, "BAM filtering begins,BAM filtering finished.")
    if flag:
        _done_file(hq_bam)
        return hq_bam
    return None


def seq_dict(fasta: Path) -> Fasta:
    return Fasta(str(fasta), sequence_always_upper=True)


def junction_caller(
    bam_file: str | Path,
    fasta: Path,
    gtf: Path,
    tmp_dir: str | Path,
    strand: int = 1,
    out_name: Optional[str] = None,
) -> Optional[str]:
    """Call splice junctions with regtools and annotate against the reference GTF."""
    bam_file = Path(bam_file)
    tmp_dir = Path(tmp_dir)
    prefix = bam_file.stem
    if not out_name:
        out_name = prefix

    if Path(f"{out_name}.janno.done").exists():
        logger.info("%s.janno found, skipping junction identification.", out_name)
        return f"{out_name}.janno"

    bed_path = tmp_dir / f"{prefix}.bed"
    cmd = (
        f"regtools junctions extract -s {strand} -i 5 -I 10000000 "
        f"{bam_file} -o {bed_path}"
    )
    bed_flag, _ = run_cmd(cmd, "Calling junctions start,Calling junctions finished!")
    if not bed_flag:
        return None

    bed = BED_handler(bed_path, tmp_dir)
    cmd = f"regtools junctions annotate {bed} {fasta} {gtf} -o {out_name}.janno"
    janno_flag, _ = run_cmd(cmd, f"{out_name}.janno generated!")
    if janno_flag:
        _remove(bed)
        _done_file(f"{out_name}.janno")
        return f"{out_name}.janno"
    return None


def junction_overlap_CDS_to_position_BED(
    janno: str,
    fasta: Path,
    gtf: Path,
    tmp_dir: str | Path,
    ao_cutoff: int = 3,
) -> tuple[Optional[str], Optional[str]]:
    """Intersect junctions with CDS annotation to identify exitron candidates.

    CDS regions are extracted on demand from *gtf* and cached alongside it as
    ``<gtf_stem>.CDS.bed`` for fast subsequent runs.
    """
    cds = extract_cds_bed(gtf)
    genome_seq = seq_dict(fasta)

    logger.info("Reading %s", janno)
    tmp_dir = Path(tmp_dir)
    rnd_id = secrets.token_hex(16)
    junction_bed = tmp_dir / f"{rnd_id}.junction.bed"
    total_junctions = 0

    with open(janno) as f, open(junction_bed, "w") as out:
        f.readline()  # skip header
        for line in f:
            l = line.rstrip().split("\t")
            total_junctions += int(l[4])
            chrm = l[0]
            start = int(l[1])
            end = int(l[2])
            stats = l[10]
            strand = l[5]
            if stats == "N" and strand != "?":
                if strand == "+":
                    left_site = genome_seq[chrm][start : start + 2].seq
                    right_site = genome_seq[chrm][end - 3 : end - 1].seq
                else:
                    left_site = genome_seq[chrm][end - 3 : end - 1].reverse.complement.seq
                    right_site = genome_seq[chrm][start : start + 2].reverse.complement.seq
                l[6] = f"{left_site}-{right_site}"
                if l[6].upper() in CANONICAL_SPLICE_SITES:
                    out.write("\t".join(l[:7]) + "\n")

    overlap_file = Path(tmp_dir) / f"{rnd_id}.overlap.bed"
    cmd = f"bedtools intersect -s -wo -a {junction_bed} -b {cds} > {overlap_file}"
    run_cmd(cmd, "Junctions intersect with CDS,Junctions intersect with CDS finished!")

    if not overlap_file.exists() or overlap_file.stat().st_size == 0:
        _remove(overlap_file)
        _remove(junction_bed)
        logger.info("No overlaps found in %s and gencode CDS", janno)
        return None, None

    tmp_dict: OrderedDict[str, str] = OrderedDict()
    with open(overlap_file) as f:
        for line in f:
            l = line.rstrip().split("\t")
            chrm = l[0]
            length = int(l[-1])
            junc_start = int(l[1])
            junc_end = int(l[2])
            junc_id = l[3]
            junc_read_no = l[4]
            strand = l[5]
            splice_site = l[6]
            ref_start = int(l[8])
            ref_end = int(l[9])
            gene_id = l[10]
            gene_name = l[11]
            pos_key = f"{chrm}:{junc_start}-{junc_end}"
            if (
                length == junc_end - junc_start
                and junc_start > ref_start
                and junc_end < ref_end
                and chrm in NON_MITO_CHRMS
                and int(junc_read_no) >= ao_cutoff
                and pos_key not in tmp_dict
            ):
                tmp_dict[pos_key] = (
                    f"{chrm}\t{junc_start}\t{junc_end}\t{junc_id}\t{junc_read_no}\t"
                    f"{strand}\t{gene_name}\t{length - 1}\t{splice_site}\t{gene_id}\t"
                    f"{total_junctions}"
                )

    _remove(overlap_file)
    _remove(junction_bed)

    if not tmp_dict:
        logger.info("No exitron found in %s", janno)
        return None, None

    base = Path(janno).stem
    position_bed_file = f"{base}.position.bed"
    src_exitron_file = f"{base}.src"

    position_set: set[str] = set()
    with open(src_exitron_file, "w") as src_out, open(position_bed_file, "w") as pos_out:
        for info in tmp_dict.values():
            src_out.write(info + "\n")
            parts = info.split("\t")
            chrm = parts[0]
            junc_start = int(parts[1])
            junc_end = int(parts[2])
            middle_point = int(np.median([junc_start, junc_end]))
            for pos in (junc_start, junc_end, middle_point):
                key = f"{chrm}\t{pos}"
                if key not in position_set:
                    pos_out.write(f"{chrm}\t{pos - 1}\t{pos}\n")
                    position_set.add(key)

    return src_exitron_file, position_bed_file


def percent_spliced_out(
    bam_file: str | Path,
    src_exitron_file: str,
    position_bed_file: str,
    ao_cutoff: int,
    pso_cutoff: float,
    mapq: int,
    out: IO[str],
) -> None:
    """Compute PSO/PSI for each exitron candidate and write passing records."""
    logger.info("Reading BAM file: %s", bam_file)
    depth_dict: dict[str, int] = {}
    cmd = f"samtools bedcov {position_bed_file} {bam_file} -Q {mapq}"
    depth_flag, result = run_cmd(cmd, "Calculate PSO and PSI.")

    if not depth_flag or result is None:
        raise RuntimeError("samtools bedcov failed; cannot compute PSO/PSI")

    for line in result.decode(errors="replace").splitlines():
        if line:
            chrm, _, pos, depth = line.rstrip().split()
            depth_dict[f"{chrm}\t{pos}"] = int(depth)
    with open(src_exitron_file) as f:
        for line in f:
            l = line.rstrip("\n").split("\t")
            chrm = l[0]
            start = int(l[1])
            end = int(l[2])
            ao = int(l[4])
            strand = l[5]
            middle_point = int(np.median([start, end]))

            if strand == "+":
                five_prime = depth_dict.get(f"{chrm}\t{start}", 0) - ao
                three_prime = depth_dict.get(f"{chrm}\t{end}", 0) - ao
                middle = depth_dict.get(f"{chrm}\t{middle_point}", 0) - ao
            else:
                five_prime = depth_dict.get(f"{chrm}\t{end}", 0) - ao
                three_prime = depth_dict.get(f"{chrm}\t{start}", 0) - ao
                middle = depth_dict.get(f"{chrm}\t{middle_point}", 0) - ao

            if five_prime < 0 or three_prime < 0:
                continue

            ave_dp = (five_prime + three_prime + middle) / 3.0
            try:
                pso = float(ao) / (ave_dp + ao)
            except ZeroDivisionError:
                logger.warning("ZeroDivisionError at %s %s %s", chrm, start, end)
                pso = 0.0

            psi = 1.0 - float(f"{pso:.3g}")
            dp = int(ao / pso) if pso > 0 else 0

            if ao >= ao_cutoff and pso >= pso_cutoff:
                out.write(
                    f"{chr(9).join(l[:-1])}\t{pso:.3g}\t{psi}\t{dp}\t{l[-1]}\n"
                )

    _remove(src_exitron_file)
    _remove(position_bed_file)
    logger.info("Finished reading BAM file: %s", bam_file)
