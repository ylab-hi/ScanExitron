# ScanExitron

<p align="center">
  <img src="https://raw.githubusercontent.com/ylab-hi/ScanExitron/main/docs/logo.png" alt="ScanExitron Logo" width="200" onerror="this.style.display='none'">
</p>

<p align="center">
  <strong>A computational workflow for exitron splicing identification from RNA-Seq data.</strong>
</p>

<p align="center">
  <a href="https://pypi.org/project/scanexitron/"><img src="https://img.shields.io/pypi/v/scanexitron.svg?logo=pypi&logoColor=white&color=blue" alt="PyPI Version"></a>
  <a href="https://anaconda.org/bioconda/scanexitron"><img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&logo=anaconda" alt="Bioconda"></a>
  <a href="https://github.com/ylab-hi/ScanExitron/actions"><img src="https://img.shields.io/github/actions/workflow/status/ylab-hi/ScanExitron/tests.yml?branch=main&logo=github&label=tests" alt="Tests Status"></a>
  <a href="https://github.com/ylab-hi/ScanExitron/blob/main/LICENSE"><img src="https://img.shields.io/github/license/ylab-hi/ScanExitron.svg?logo=open-source-initiative&logoColor=white&color=green" alt="License"></a>
</p>

---

## 🌟 Features

- **Accurate Identification**: High-performance detection of exitrons (exonic introns) directly from aligned RNA-Seq reads (BAM/CRAM).
- **VCF Conversion**: Built-in support to convert results into standard VCF files for downstream variant analysis.
- **Bioconda Support**: Easy deployment along with all underlying command-line dependencies.

---

## 🚀 Installation

### Option 1: Via Bioconda (Recommended)
This installs `scanexitron` along with all required external compiled tools automatically:
```bash
conda install -c bioconda scanexitron
```

### Option 2: Via PyPI
Installs the Python package. *Note: You must install the external dependencies separately (see below).*
```bash
pip install scanexitron
```

### Option 3: From Source
```bash
git clone https://github.com/ylab-hi/ScanExitron.git
cd ScanExitron
pip install -e ".[dev]"
```

---

## 🛠️ External Dependencies

If you did not install via `conda`, make sure the following bioinformatics tools are installed and available on your system `PATH`:

| Dependency | Required Version | Conda Install Command |
| :--- | :---: | :--- |
| 🛠️ **[regtools](https://github.com/griffithlab/regtools)** | `0.5.0` | `conda install -c bioconda regtools=0.5.0` |
| 🧬 **[samtools](http://www.htslib.org/)** | `≥ 1.10` | `conda install -c bioconda samtools` |
| 🗃️ **[bedtools](https://bedtools.readthedocs.io/)** | `≥ 2.26` | `conda install -c bioconda bedtools` |

---

## 📂 Reference Data Setup

Download your genome FASTA and matching transcript annotation files. For example, to setup `hg38` reference resources:

```bash
# Download and decompress the hg38 genome
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Download and decompress the Gencode release annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
gunzip gencode.v37.annotation.gtf.gz
```

---

## 💻 Usage

### 🔍 1. Exitron Calling
Run `scanexitron` by pointing it to your reference FASTA, transcript annotation GTF, and the input BAM/CRAM file (which must have a `.bai`/`.crai` index file alongside it):

```bash
scanexitron run -i input.bam -r hg38.fa -g gencode.v37.annotation.gtf
```

#### Calling Options
| Flag | Default | Description |
| :--- | :---: | :--- |
| **`-i, --input`** | *Required* | Path to the input BAM/CRAM file. |
| **`-r, --ref`** | *Required* | Path to the reference genome FASTA file. |
| **`-g, --gtf`** | *Required* | Path to the annotation GTF file. |
| **`-a, --ao`** | `3` | Minimum junction-spanning reads supporting the exitron. |
| **`-p, --pso`** | `0.05` | Minimum Percent Spliced Out (PSO) value. |
| **`-m, --mapq`** | `50` | Minimum mapping quality for alignment filtering. |
| **`-s, --strand`** | `1` | Strand specificity: `0` = unstranded, `1` = first-strand/RF, `2` = second-strand/FR. |
| **`-t, --threads`** | `1` | Number of threads to allocate for samtools. |
| **`-o, --output`** | *Input Stem* | Output prefix. |
| **`--verbose`** | *Off* | Enable verbose logging for debugging. |

---

### 🔄 2. Convert to VCF
Format your exitron tabular output into a standard VCF file using `scanexitron convert`:

```bash
scanexitron convert -i sample.exitron -r hg38.fa -o sample.vcf
```

#### VCF Converter Options
| Flag | Default | Description |
| :--- | :---: | :--- |
| **`-i, --input`** | *Required* | Path to the input tabular results (`.exitron` file). |
| **`-r, --ref`** | *Required* | Path to the reference genome FASTA file. |
| **`-o, --output`** | `output.vcf` | Path for the output VCF file. |

---

## 📊 Output Format

Results are written to `<input>.exitron` as a tab-delimited file with the following schema:

| Column | Type | Description |
| :--- | :---: | :--- |
| **`chrom`** | String | Chromosome name |
| **`start`** | Integer | Start position (0-based, half-open) |
| **`end`** | Integer | End position (0-based, half-open) |
| **`name`** | String | Junction identifier |
| **`ao`** | Integer | Number of supporting junction reads |
| **`strand`** | String | Genomic strand (`+` or `-`) |
| **`gene_symbol`** | String | HGNC gene symbol |
| **`length`** | Integer | Exitron length in base pairs |
| **`splice_site`** | String | Donor–acceptor dinucleotides (e.g., `GT-AG`) |
| **`gene_id`** | String | Ensembl gene ID |
| **`pso`** | Float | Percent Spliced Out (relative exitron abundance) |
| **`psi`** | Float | Percent Spliced In |
| **`dp`** | Float | Average local sequencing depth |
| **`total_junctions`** | Integer | Total junction-spanning reads detected in the sample |

*Note: RegTools intermediate files (`<input>.janno`) are preserved in the directory for manual inspection.*

---

## 👩‍💻 Development

Contributions and local development setups are welcome:

```bash
# Install development and testing dependencies
pip install -e ".[dev]"

# Run test suite
pytest
```

---

## 📖 Citation

If you use ScanExitron in your research, please cite our papers:
* **Molecular Cell**: [Molecular Cell Publication](https://www.sciencedirect.com/science/article/pii/S1097276521002239)
* **STAR Protocols**: [STAR Protocols Publication](https://star-protocols.cell.com/protocols/977)

---

## 📄 License

This project is licensed under the terms of the **[MIT License](file:///Users/twp7981/dev/ScanExitron/LICENSE)**.

---

## 📬 Contact & Support

- **Bug Reports & Feature Requests**: Please submit an issue on the [GitHub Issues tracker](https://github.com/ylab-hi/ScanExitron/issues).
