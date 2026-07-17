# Running ScanExitron

Learn how to perform exitron identification and convert your results into standard VCF files.

!!! info "Learning Objectives"
    By the end of this tutorial, you will be able to:
    - Run the exitron calling workflow with custom parameters
    - Convert calling results into VCF format
    - Interpret the results and quality logs

    **Prerequisites**:
    - ScanExitron installed
    - Coordinate-sorted and indexed BAM file (see [Reads Alignment](reads-alignment.md))
    - Reference genome FASTA and matching transcript GTF

---

## 1. Splicing Calling (`scanexitron run`)

The main command `scanexitron run` detects exitrons by comparing spliced alignments against transcript models.

### Basic command

```bash
scanexitron run -i sample.bam -r hg38.fa -g gencode.v37.annotation.gtf
```

### Customizing thresholds

Depending on your sequencing depth and library type, you may want to customize parameters:

- **Library Strand Specificity (`-s`)**: 
  - `0` = unstranded
  - `1` = first-strand/RF (default, typical for Illumina TruSeq stranded)
  - `2` = second-strand/FR (typical for dUTP-based protocols)
- **Minimum supporting reads (`-a`)**: By default, at least `3` junction-spanning reads must support an exitron. Increase this for deep sequencing to reduce false positives.
- **Minimum Percent Spliced Out (`-p`)**: By default, `0.05` (5%). Events with PSO below this threshold are filtered out.
- **Mapping Quality (`-m`)**: By default, reads with MAPQ < `50` are filtered out during alignment processing.

Example with custom settings:

```bash
scanexitron run -i sample.bam \
                -r hg38.fa \
                -g gencode.v37.annotation.gtf \
                -a 5 \
                -p 0.10 \
                -s 1 \
                -t 4 \
                -o sample_exitrons
```

This will run the pipeline using 4 threads, requiring 5 supporting reads and 10% minimum PSO, writing outputs to `sample_exitrons.exitron`.

---

## 2. Converting Results to VCF (`scanexitron convert`)

Once exitrons are called, you can format them into a standard VCF (Variant Call Format) file using `scanexitron convert`.

```bash
scanexitron convert -i sample_exitrons.exitron -r hg38.fa -o sample_exitrons.vcf
```

This will output a standard VCF file where each exitron splicing event is represented as a deletion. The VCF file will include annotations in the INFO field such as the PSO, PSI, and supporting reads (`AO`).

---

## Summary

You've learned how to:
- ✅ Run exitron identification with custom filters using `scanexitron run`
- ✅ Format tabular results to standard VCF format using `scanexitron convert`

!!! success "Exitron Splicing Data Ready!"
    Your exitron predictions are now ready for functional downstream analysis and variant annotations!
