# ScanExitron
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

License
----------------
The project is licensed under the [MIT license](https://opensource.org/licenses/MIT).

Contact
-----------------
Bug reports or feature requests can be submitted on the <a href="https://github.com/ylab-hi/ScanExitron/issues">ScanExitron Github page</a>.
