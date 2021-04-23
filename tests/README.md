Testing environment
---
* [Python](https://www.python.org) (v3.7)
* [SamTools](http://www.htslib.org/) (v1.10)
* [BedTools](https://bedtools.readthedocs.io/en/latest/) (v2.26.0)
* [pyfaidx](https://github.com/mdshw5/pyfaidx) (v0.5.9)
* [pysam](https://github.com/pysam-developers/pysam) (v0.15.3)
* [bedops](https://bedops.readthedocs.io/en/latest/) (v2.4.38)
* [RegTools](https://github.com/griffithlab/regtools) (v0.4.2)


Sample data for testing
---
This directory contains one sample dataset for testing purpose to make sure you have installed ScanExitron and its dependencies sucessfully. __test.bam__ is a RNA-seq dataset (in BAM format) contains two exitrons. 

```
ScanExitron.py -i test.bam -r hg38
```

Output file (test.exitron)
---
|chrom|start|end|name|ao|strand|gene_symbol|length|splice_site|gene_id|pso|psi|dp|total_junctions|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|chr17|50071324|50071418|JUNC00000010|16|+|ITGA3|93|GT-AG|ENSG00000005884.15|0.02674|0.97326|598|104499|
|chr17|58302312|58302408|JUNC00000893|3|-|BZRAP1|95|GT-AG|ENSG00000005379.13|0.03371|0.96629|89|104499|
