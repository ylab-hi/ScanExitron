# ScanExitron
[![Build Status](https://travis-ci.org/ylab-hi/ScanExitron.svg?branch=master&status=passed)](https://travis-ci.org/ylab-hi/ScanExitron)

A computational workflow for exitron splicing identification


Prerequisites
----------------
You need Python 3.6 to run ScanExitron.

### install necessary python packages via anaconda
Install [anaconda](https://www.anaconda.com/download/) (python 3.6) firstly, then install dependent packages via conda in bioconda channel.
```
conda install -c bioconda samtools
conda install -c bioconda bedtools
conda install -c bioconda regtools
conda install -c bioconda pyfaidx
conda install -c bioconda pysam
conda install -c bioconda bedops
```

## Prepare the human genome FASTA sequences and annotation GTF file.

```
# hg38 genome
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# hg19 genome
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

# hg38 annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz
gunzip gencode.v21.annotation.gtf.gz

# hg19 annotation
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz

# hg38 CDS
grep -P "\tCDS\t" gencode.v21.annotation.gtf | gtf2bed  > gencode.hg38.CDS.bed
# hg19 CDS
grep -P "\tCDS\t" gencode.v19.annotation.gtf | gtf2bed  > gencode.hg19.CDS.bed
```

### configure config.ini file
```
[fasta]

# reference genome file in FASTA format (absolute path)

hg38=/abs/path/to/hg38.fa
hg19=/abs/path/to/hg19.fa

[annotation]

# gene annotation file in GTF format (absolute path)

hg38=/abs/path/to/gencode.v21.annotation.gtf
hg19=/abs/path/to/gencode.v19.annotation.gtf

[cds]

# CDS annotation in BED format (absolute path)

hg38=/abs/path/to/gencode.hg38.CDS.bed
hg19=/abs/path/to/gencode.hg19.CDS.bed
```

Usage
-------------------------
#### Exitron calling using RNA-seq data
```
ScanExitron.py -i input_rna_seq_bam_file -r [hg38/hg19] -m mapping_quality
```

#### Options:

```	
-h, --help            show this help message and exit
-i INPUT, --input INPUT
                        RNA-seq alignment file (BAM/CRAM)
--mapq                  consider reads with MAPQ >= cutoff (default: 0)		
-r {hg19,hg38}, --ref {hg19,hg38}
                        reference genome (default: hg38)
```

#### Input:
```	
input_bam_file      :input RNA-seq BAM/CRAM file. (e.g., rna-seq.bam)
reference_genome    :specify reference genome (hg19 or hg38)
```

#### Output:
```
exitron_file			:Reported exitrons in a TAB-delimited file. (rna-seq.exitron)
```

__Report Columns__

|Column Name | Description |
| ---------- | ----------- |
|chrom       | The chromosome of this exitron|
|start       | The start position of this exitron in the zero-based, half-open coordinate system |
|end	       | The stop position of this exitron in the zero-based, half-open coordinate system |
|name        | Identifier for the junction |
|ao          | Observed supporting reads for exitron |
|strand      | The strand the exitron is identified |
|gene_symbol | The Gene symbol of the affected gene |
|length      | Length of the exitron |
|splice_site | The two basepairs at the donor and acceptor sites separated by a hyphen |
|gene_id     | The Ensembl ID of the affected gene |
|pso         | The percent spliced out (PSO) index  |
|psi         | The percent spliced in (PSI) index |
|dp          | The average depth of the exitron   |
|total_junctions | The total number of junctions in the sample |


We also keep [RegTools](https://github.com/griffithlab/regtools) interim results (rna-seq.janno) for developers.

For a detailed explanation, please refer to [The Documentation of RegTools](https://regtools.readthedocs.org/en/latest/)

License
----------------
The project is licensed under the [MIT license](https://opensource.org/licenses/MIT).

Contact
-----------------
Bug reports or feature requests can be submitted on the <a href="https://github.com/ylab-hi/ScanExitron/issues">ScanExitron Github page</a>.
