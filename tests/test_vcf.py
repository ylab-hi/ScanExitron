from pathlib import Path
from scanexitron.vcf import exitron2vcf


def test_exitron2vcf_conversion(tmp_path):
    # Create a dummy reference FASTA file
    fasta_file = tmp_path / "dummy.fa"
    # Sequence on chr17:
    # 0-based indices:
    # 012345678901234567890123456789
    # A C T G A T C G A T C G A T C G
    # Let's make it 100 bases long
    sequence = "ACTGATCGAT" * 10
    fasta_file.write_text(f">chr17\n{sequence}\n")

    # Create a dummy .exitron file
    exitron_file = tmp_path / "sample.exitron"
    exitron_content = (
        "chrom\tstart\tend\tname\tao\tstrand\tgene_symbol\tlength\tsplice_site\tgene_id\tpso\tpsi\tdp\ttotal_junctions\n"
        "chr17\t10\t20\tJUNC00000001\t15\t+\tITGA3\t9\tGT-AG\tENSG00000005884.18\t0.05\t0.95\t100\t1000\n"
    )
    exitron_file.write_text(exitron_content)

    # Output VCF file path
    vcf_file = tmp_path / "output.vcf"

    # Run conversion
    exitron2vcf(exitron_file, vcf_file, fasta_file)

    # Check VCF file exists
    assert vcf_file.exists()

    # Read and verify content
    lines = vcf_file.read_text().splitlines()

    # Find the data line (non-header)
    data_lines = [line for line in lines if not line.startswith("##") and not line.startswith("#")]
    assert len(data_lines) == 1

    # Fields: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample
    fields = data_lines[0].split("\t")
    assert len(fields) == 10
    assert fields[0] == "chr17"  # CHROM
    assert fields[1] == "10"     # POS (which is BED start = l[1])
    assert fields[2] == "JUNC00000001"  # ID

    # Let's verify REF and ALT slices:
    # l[1] = 10 (pos), l[2] = 20 (end).
    # genome_seq["chr17"][9:10] -> 10th base of sequence -> "T" (0-based index 9)
    # genome_seq["chr17"][9:19] -> 10th to 19th base -> "TACTGATCGA" (0-based indices 9 to 18)
    assert fields[3] == "TACTGATCGA"  # REF
    assert fields[4] == "T"          # ALT

    # Check INFO fields
    info = fields[7]
    assert "SVTYPE=DEL" in info
    assert "END=19" in info
    assert "AO=15" in info
    assert "DP=100" in info
    assert "STRAND=+" in info
    assert "SpliecedSite=GT-AG" in info
    assert "GeneName=ITGA3" in info
    assert "GeneID=ENSG00000005884" in info
    assert "SVLEN=-9" in info
    assert "PSO=0.05" in info
    assert "PSI=0.95" in info

    # Check FORMAT and SAMPLE fields
    assert fields[8] == "GT"
    assert fields[9] == "0/1"
