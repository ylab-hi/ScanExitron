# CLI Reference

Complete reference for the ScanExitron command-line interface.

## Command Overview

ScanExitron uses a command-line interface with two main subcommands:

- **`scanexitron run`**: Detects exitron splicing events from RNA-Seq BAM/CRAM files.
- **`scanexitron convert`**: Converts tabular exitron results into a standard VCF file.

---

## 1. `scanexitron run`

Detect exitron splicing events from RNA-Seq data.

### Usage

```bash
scanexitron run [OPTIONS]
```

### Options

| Flag | Type | Default | Description |
| :--- | :---: | :---: | :--- |
| **`-i, --input`** | Path | *Required* | Path to the input BAM/CRAM file. The coordinate-sorted index (`.bai`/`.crai`) file must be present alongside the BAM/CRAM file. |
| **`-r, --ref`** | Path | *Required* | Path to the reference genome FASTA file. |
| **`-g, --gtf`** | Path | *Required* | Path to the annotation GTF file. |
| **`-a, --ao`** | Integer | `3` | Minimum number of junction-spanning reads supporting the exitron. |
| **`-p, --pso`** | Float | `0.05` | Minimum Percent Spliced Out (PSO) value. |
| **`-m, --mapq`** | Integer | `50` | Minimum mapping quality (MAPQ) for alignment filtering. |
| **`-s, --strand`** | Integer | `1` | RNA library strand specificity:<br>`0` = unstranded<br>`1` = first-strand/RF (default)<br>`2` = second-strand/FR |
| **`-t, --threads`** | Integer | `1` | Number of threads to allocate for `samtools`. |
| **`-o, --output`** | String | *Input Stem* | Output prefix. By default, uses the input filename stem (e.g., `sample` for `sample.bam`). |
| **`--verbose`** | Flag | *Off* | Enable verbose/debug logging. |
| **`-h, --help`** | Flag | *Off* | Show the help message and exit. |

---

## 2. `scanexitron convert`

Convert an exitron results table into a standard VCF file.

### Usage

```bash
scanexitron convert [OPTIONS]
```

### Options

| Flag | Type | Default | Description |
| :--- | :---: | :---: | :--- |
| **`-i, --input`** | Path | *Required* | Path to the input tabular results (`.exitron` file). |
| **`-r, --ref`** | Path | *Required* | Path to the reference genome FASTA file. |
| **`-o, --output`** | String | `output.vcf` | Path for the output VCF file. |
| **`-h, --help`** | Flag | *Off* | Show the help message and exit. |
