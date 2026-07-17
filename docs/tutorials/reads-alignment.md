# Reads Alignment

Learn how to perform read alignment and prepare BAM files for use with ScanExitron.

!!! info "Learning Objectives"
    By the end of this tutorial, you will be able to:
    - Prepare reference genomes and annotations
    - Align short-read RNA-seq data using STAR
    - Align long-read RNA-seq data using Minimap2
    - Sort and index output BAM files

    **Prerequisites**:
    - SAMtools installed
    - Aligners installed (STAR or Minimap2)

    **Time**: ~30 minutes to several hours, depending on dataset size.

---

## 1. Reference Data Preparation

Download the reference genome FASTA and the annotation GTF matching your species. For example, to setup `hg38` reference resources:

```bash
# Genome FASTA
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Transcript Annotation GTF (Gencode)
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
gunzip gencode.v37.annotation.gtf.gz
```

---

## 2. Short-read RNA-seq Alignment (STAR)

For Illumina short-read RNA-seq data, we recommend using **STAR** (Spliced Transcripts Alignment to a Reference) as it accurately identifies splicing junctions.

### Step 2.1: Generate STAR Genome Index

```bash
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir hg38_star_index \
     --genomeFastaFiles hg38.fa \
     --sjdbGTFfile gencode.v37.annotation.gtf \
     --sjdbOverhang 99  # (Read length - 1)
```

### Step 2.2: Run Spliced Alignment

```bash
STAR --runThreadN 8 \
     --genomeDir hg38_star_index \
     --readFilesIn read1.fastq.gz read2.fastq.gz \
     --readFilesCommand gunzip -c \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix sample_star_
```

This generates `sample_star_Aligned.sortedByCoord.out.bam`.

---

## 3. Long-read RNA-seq Alignment (Minimap2)

For third-generation long-read RNA-seq data (PacBio Iso-Seq or ONT Nanopore), **Minimap2** is the standard aligner.

### Step 3.1: Run Alignment

```bash
# For PacBio Iso-Seq (high accuracy)
minimap2 -ax splice -t 8 hg38.fa reads.fastq.gz | samtools sort -@ 4 -o sample_minimap2.bam -

# For Oxford Nanopore (ONT) Direct RNA-seq
minimap2 -ax splice -uf -k14 -t 8 hg38.fa reads.fastq.gz | samtools sort -@ 4 -o sample_minimap2.bam -
```

---

## 4. Sorting and Indexing

ScanExitron requires coordinate-sorted BAM files with an index (`.bai`) located alongside them.

If your BAM file is not sorted and indexed:

```bash
# Sort by coordinate
samtools sort -@ 4 -o sorted_output.bam input.bam

# Generate index (.bai)
samtools index sorted_output.bam
```

Verify that the index file `sorted_output.bam.bai` is present in the same directory.

---

## Summary

You've learned how to:
- ✅ Prepare reference datasets (`hg38.fa` and annotation `gtf`)
- ✅ Generate spliced alignments for short reads (STAR) or long reads (Minimap2)
- ✅ Coordinate-sort and index output BAMs

!!! success "Alignment Files Ready!"
    Your BAM and BAI files are now ready to be processed with ScanExitron!
