Testing environment
---
* [Python](https://www.python.org) (v3.10 to v3.12, tested with v3.10.14)
* [SamTools](http://www.htslib.org/) (≥ v1.10, tested with v1.21)
* [BedTools](https://bedtools.readthedocs.io/en/latest/) (≥ v2.26.0, tested with v2.31.1)
* [pyfaidx](https://github.com/mdshw5/pyfaidx) (≥ v0.7.0, tested with v0.9.0.4)
* [RegTools](https://github.com/griffithlab/regtools) (v0.5.0, tested with v0.5.0)


Sample data for testing
---
This directory contains one sample dataset for testing purpose to make sure you have installed ScanExitron and its dependencies successfully. __test.bam__ is an RNA-seq dataset (in BAM format) that contains two exitrons.

```bash
scanexitron run -i test.bam -r hg38.fa -g gencode.gtf -a 3 -p 0.01 -s 0
```

Output file (test.exitron)
---
|chrom|start|end|name|ao|strand|gene_symbol|length|splice_site|gene_id|pso|psi|dp|total_junctions|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|chr17|50071324|50071418|JUNC00000010|16|+|ITGA3|93|GT-AG|ENSG00000005884.18|0.0267|0.9733|598|104499|
|chr17|58302312|58302408|JUNC00000893|3|-|TSPOAP1|95|GT-AG|ENSG00000005379.17|0.0337|0.9663|89|104499|
