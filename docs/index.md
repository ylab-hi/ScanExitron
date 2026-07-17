# ScanExitron

<div class="hero" markdown>

## ScanExitron: Identification of Exitron Splicing from RNA-Seq Data

A powerful computational workflow for detecting exitron (exonic intron) splicing events directly from aligned RNA-Seq reads (BAM/CRAM) at base resolution.

[Get Started](getting-started/quick-start.md){ .md-button .md-button--primary }
[View on GitHub](https://github.com/ylab-hi/ScanExitron){ .md-button }

</div>

______________________________________________________________________

## :material-star: Key Features

<div class="feature-grid" markdown>

<div class="feature-item" markdown>
### :material-target: High Accuracy
High-performance detection and quantification of exitrons directly from aligned reads (BAM/CRAM) using junction-spanning coverage.
</div>

<div class="feature-item" markdown>
### :material-sync: VCF Conversion
Convert tabular results into standard VCF format for seamless integration into downstream variant calling and annotation pipelines.
</div>

<div class="feature-item" markdown>
### :material-package-variant-closed: Bioconda Ready
Available on Bioconda with all complex underlying command-line tools (such as `regtools` and `samtools`) pre-packaged and auto-installed.
</div>

</div>

______________________________________________________________________

## Quick Start

Get up and running with ScanExitron in under 5 minutes:

```bash
# Install ScanExitron (via Bioconda)
conda install -c bioconda scanexitron

# Run ScanExitron calling
scanexitron run -i input.bam -r hg38.fa -g gencode.v37.annotation.gtf
```

Ready to dive in? Check out our [Quick Start Guide](getting-started/quick-start.md).

______________________________________________________________________

## What is an Exitron?

Exitrons (exonic introns) are internal regions of protein-coding exons that are spliced out. Although they exhibit typical splicing characteristics (canonical splice sites and splicing regulatory elements), their retention or splicing out can alter protein structure, affect function, or introduce premature stop codons leading to nonsense-mediated decay (NMD). ScanExitron provides an efficient workflow to detect, quantify, and analyze these events from standard transcriptomic datasets.

______________________________________________________________________

## Citation

If you use ScanExitron in your research, please cite our papers:

* **Molecular Cell**: Wang, Ting-You, et al. "Molecular Cell Publication" (2021). [DOI: 10.1016/j.molcel.2021.03.013](https://doi.org/10.1016/j.molcel.2021.03.013)
* **STAR Protocols**: Wang, Ting-You, et al. "STAR Protocols Publication" (2021). [DOI: 10.1016/j.xpro.2021.100580](https://doi.org/10.1016/j.xpro.2021.100580)

See [Citation](about/citation.md) for full citation formats and BibTeX entries.

______________________________________________________________________

## License

ScanExitron is licensed under the MIT License. See [License](about/license.md) for details.
