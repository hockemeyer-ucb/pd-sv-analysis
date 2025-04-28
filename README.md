# iSCORE-PD SNV analysis

## Pre-requesite

### Software

- [R](https://cran.r-project.org/)
- [Java Runtime (JRE)](https://www.oracle.com/java/technologies/downloads/?er=221886)
- [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- [wget](https://www.gnu.org/software/wget/) or [curl](https://curl.se/)
- C compiler with make if running from Linux

### R packages

#### Bioconductor

**From R terminal:**

```{R}
install.packages("BiocManager")
```

#### For analyses:

- httpgd
- VariantAnnotation
- data.table
- TxDb.Hsapiens.UCSC.hg38.knownGene
- BSgenome.Hsapiens.UCSC.hg38
- org.Hs.eg.db
- parallel
- DT
- htmltools
- ggplot2
- fastreeR
- ape
- gridExtra

**From R terminal:**

```{R}
BiocManager::install(c("httpgd",
  "VariantAnnotation",
  "data.table",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "BSgenome.Hsapiens.UCSC.hg38",
  "org.Hs.eg.db",
  "parallel",
  "DT",
  "htmltools",
  "ggplot2",
  "fastreeR",
  "ape",
  "gridExtra"))
```

#### For kniting the RMarkdown document

- rmarkdown

**From R terminal:**

```{R}
install.packages("rmarkdown")
```

## Installation

Clone the repo:

```{bash}
git clone https://github.com/hockemeyer-ucb/pd-sv-analysis.git
```

## Download test data

To test your installation, you can download the joint genotype for the iSCORE-PD cell lines for chromosome 22
in a directory named SNV within the pd-sv-analysis repo.

**wget**

```{bash}
wget -P pd-sv-analysis/SNV 'https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr22.vcf.gz'
```

**curl**

```{bash}
curl --create-dirs -O --output-dir pd-sv-analysis/SNV 'https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr22.vcf.gz'
```

## Knit
