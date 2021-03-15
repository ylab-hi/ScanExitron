# Applications of ScanExitron

We applied ScanExitron on RNA-Seq data of TCGA project.

### Data preparation
RNA-Seq data in BAM format from TCGA cohort are available on [Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov).

### Install ScanExitron and its dependences
Following the instructions on [ScanExitron GitHub](https://github.com/ylab-hi/ScanExitron) to install ScanExitron.

### Use ScanExitron to analyze RNA-Seq BAM file
We used __genCmds.py__ to generate commands to run ScanExitron


## Additional analysis
Besides RNA-Seq data, some TCGA cohorts have corresponding Mass spectra data available at [CPTAC portal](https://cptac-data-portal.georgetown.edu/cptacPublic/).
CPTAC proteomic data can use to confirm the expression of peptides derived from Exitron splicing.

### Exitron-derived polypeptides identification using OpenMS
1. Create decoy sequences database.
```
DecoyDatabase -in <in.fasta> -out <db.fasta>
```
2. Search CPTAC mass spectra.
```
MSGFPlusAdapter -ini <config.ini> -in <spectra.mzML> -out <out.idXML> -database <db.fasta> -executable <MSGFPlus.jar> -java_memory 20000 -threads 16
```
3. Refresh the mapping of peptides to proteins and add target/decoy information.
```
PeptideIndexer -in <out.idXML> -fasta <db.fasta> -out <pi_out.idXML> -allow_unmatched -enzyme:specificity 'semi'
```
4. Merge peptide identification files from multiple runs.
```
IDMerger -in <pi_out.idXML files> -out <merged.idXML>
```
5. Control for false discovery rate
```
FalseDiscoveryRate -in <merged.idXML> -out <fdr_out.idXML>
IDFilter -in <fdr_out.idXML> -out <fdr_filtered.idXML> -score:pep 0.05
```
We have packed the OpenMS commonds to a simple Python script [cptac.py](https://github.com/ylab-hi/ScanExitron/blob/master/Applications/cptac.py), and the corresponding OpenMS configure file [config.ini](https://github.com/ylab-hi/ScanExitron/blob/master/Applications/config/config.ini) is also available.


