# ScanExitron

A computational workflow for exitron splicing identification from RNA-Seq data.

## Installation

### From PyPI

```bash
pip install scanexitron
```

### From Bioconda (recommended — includes external tool dependencies)

```bash
conda install -c bioconda -c conda-forge scanexitron
```

### From source

```bash
git clone https://github.com/ylab-hi/ScanExitron.git
cd ScanExitron
pip install -e ".[dev]"
```

## External Dependencies

ScanExitron requires three bioinformatics tools on your `PATH`:

| Tool | Version | Install |
|------|---------|---------|
| [regtools](https://github.com/griffithlab/regtools) | 0.5.0 | `conda install -c bioconda regtools=0.5.0` |
| [samtools](http://www.htslib.org/) | ≥1.10 | `conda install -c bioconda samtools` |
| [bedtools](https://bedtools.readthedocs.io/) | ≥2.26 | `conda install -c bioconda bedtools` |

## Reference Data Setup

Download genome FASTA and annotation files, then derive the CDS BED:

```bash
# hg38 genome
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz && gunzip hg38.fa.gz

# hg38 Gencode annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
gunzip gencode.v37.annotation.gtf.gz

# hg38 CDS BED (derived from GTF)
awk 'OFS="\t" {if ($3=="CDS") {print $1,$4-1,$5,$10,$16,$7}}' \
    gencode.v37.annotation.gtf | tr -d '";' > gencode.hg38.CDS.bed

# hg19 genome
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz && gunzip hg19.fa.gz

# hg19 Gencode annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz

# hg19 CDS BED
awk 'BEGIN{OFS="\t"} {if ($3=="CDS") {if ($13=="ccdsid") {print $1,$4-1,$5,$20,$22,$7} \
    else {print $1,$4-1,$5,$18,$20,$7}}}' \
    gencode.v19.annotation.gtf | tr -d '";' > gencode.hg19.CDS.bed
```

## Configuration

Copy `config.ini.example` to `config.ini` (in your working directory or `~/.config/scanexitron/config.ini`) and fill in the absolute paths:

```ini
[fasta]
hg38 = /abs/path/to/hg38.fa
hg19 = /abs/path/to/hg19.fa

[annotation]
hg38 = /abs/path/to/gencode.v37.annotation.gtf
hg19 = /abs/path/to/gencode.v19.annotation.gtf

[cds]
hg38 = /abs/path/to/gencode.hg38.CDS.bed
hg19 = /abs/path/to/gencode.hg19.CDS.bed
```

## Usage

### Exitron calling

```bash
scanexitron -i input.bam -r hg38
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | — | Input BAM/CRAM file (required; index must be alongside) |
| `-r / --ref` | `hg38` | Reference genome: `hg38` or `hg19` |
| `-a / --ao` | `3` | Minimum supporting reads for an exitron |
| `-p / --pso` | `0.05` | Minimum PSO value |
| `-m / --mapq` | `50` | Minimum mapping quality |
| `-s / --strand` | `1` | Strand: 0=unstranded, 1=first-strand/RF, 2=second-strand/FR |
| `-t / --threads` | `1` | Threads for samtools |
| `-c / --config` | `config.ini` | Path to config file |
| `--verbose` | off | Enable debug logging |

### Convert to VCF

```bash
exitron2vcf -i sample.exitron -o sample.vcf -r hg38
```

## Output

Results are written to `<input>.exitron` (tab-delimited):

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| start | Start position (0-based, half-open) |
| end | End position (0-based, half-open) |
| name | Junction identifier |
| ao | Supporting reads |
| strand | Strand |
| gene_symbol | Gene symbol |
| length | Exitron length (bp) |
| splice_site | Donor–acceptor dinucleotides (e.g. `GT-AG`) |
| gene_id | Ensembl gene ID |
| pso | Percent spliced out |
| psi | Percent spliced in |
| dp | Average depth |
| total_junctions | Total junctions in the sample |

RegTools intermediate results (`<input>.janno`) are retained for inspection.

## Development

```bash
pip install -e ".[dev]"
pytest
```

## License

[MIT](LICENSE)

## Citation

Please cite our papers at [Molecular Cell](https://www.sciencedirect.com/science/article/pii/S1097276521002239) and [STAR protocols](https://star-protocols.cell.com/protocols/977).

## Contact

Bug reports and feature requests: [GitHub Issues](https://github.com/ylab-hi/ScanExitron/issues)
