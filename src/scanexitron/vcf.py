import logging
from pathlib import Path
from pyfaidx import Fasta

logger = logging.getLogger(__name__)

_VCF_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##FILTER=<ID=LowQual,Description="PE/SR support below 3 or mapping quality below 20.">
##INFO=<ID=DP,Number=2,Type=Integer,Description="Total depth of junction">
##INFO=<ID=SpliecedSite,Number=2,Type=String,Description="Spliced site">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Junction strand">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=AO,Number=1,Type=Integer,Description="Reads support of the exitron">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=PSO,Number=1,Type=Integer,Description="Percent spliced-out">
##INFO=<ID=GeneName,Number=1,Type=String,Description="Gene name">
##INFO=<ID=GeneID,Number=1,Type=String,Description="Gene ID">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
##INFO=<ID=CE,Number=1,Type=Float,Description="Consensus sequence entropy">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Predicted length of the insertion">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Predicted microhomology length using a max. edit distance of 2">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def exitron2vcf(
    in_file: str | Path,
    out_vcf: str | Path,
    fasta: Path,
) -> None:
    """Convert an exitron results table to VCF format."""
    genome_seq = Fasta(str(fasta), sequence_always_upper=True, as_raw=True)
    sample = Path(in_file).stem

    with open(out_vcf, "w") as outfile:
        outfile.write(_VCF_HEADER.format(sample=sample))
        with open(in_file) as f:
            f.readline()  # skip header
            for line in f:
                l = line.rstrip().split("\t")
                chrom = l[0]
                pos = int(l[1])
                end = int(l[2])
                _alt = genome_seq[chrom][pos - 1 : pos]
                _ref = genome_seq[chrom][pos - 1 : end - 1]
                idx = l[3]
                ao = int(l[4])
                strand = l[5]
                gene_name = l[6]
                length = int(l[7])
                spliced_site = l[8]
                gene_id = l[9].split(".")[0]
                pso = float(l[10])
                psi = float(l[11])
                dp = int(l[12])
                outfile.write(
                    f"{chrom}\t{pos}\t{idx}\t{_ref}\t{_alt}\t.\t.\t"
                    f"SVTYPE=DEL;END={end - 1};AO={ao};DP={dp};STRAND={strand};"
                    f"SpliecedSite={spliced_site};GeneName={gene_name};GeneID={gene_id};"
                    f"SVLEN=-{length};PSO={pso};PSI={psi}\tGT\t0/1\n"
                )
