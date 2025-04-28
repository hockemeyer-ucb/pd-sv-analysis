# iSCORE-PD SNV analysis

## Pre-requesite

### Software

- [R](https://cran.r-project.org/)
- [Java Runtime (JRE)](https://www.oracle.com/java/technologies/downloads/?er=221886)
- [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- [wget](https://www.gnu.org/software/wget/) or [curl](https://curl.se/)
- C compiler with make if running from Linux

### Additional libraries if installing from vanilla ubuntu

```{bash}
## Update apt packages
sudo apt update -y

## Install Ubuntu libraries and dependencies
sudo apt install -y libssl-dev /
  libcurl4-openssl-dev /
  libpng-dev /
  libxml2-dev /
  liblzma-dev /
  libbz2-dev /
  libblas-dev /
  liblapack-dev /
  gfortran /
  default-jre /
  default-jdk /
  libpcre2-dev /
  libdeflate-dev /
  libzstd-dev /
  libtirpc-dev /
  pandoc /
  libfontconfig1-dev /
  libtiff-dev /
  libcairo2-dev

## Update location of java for R deployment
sudo R CMD javareconf
```

### R packages

#### Bioconductor

**From R terminal:**

```{R}
install.packages("BiocManager")
```

#### For analyses:

- remotes
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
options(Ncpus = parallel::detectCores())

BiocManager::install(c("remotes",
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

remotes::install_github("nx10/httpgd")
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

To test your installation, you can download the iSCORE-PD cell lines joint genotypes for chromosome 22
in a directory named SNV within the pd-sv-analysis repo.

**wget**

```{bash}
wget -P pd-sv-analysis/SNV 'https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr22.vcf.gz'
```

**curl**

```{bash}
curl --create-dirs -O --output-dir pd-sv-analysis/SNV 'https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr22.vcf.gz'
```

## Suggested hardware

- OS: Ubuntu 24.04 LTS
- Architecture: ARM64
- Diskspace: 25 Gb
- CPU: 16
- RAM: 32 Gb
- Test dataset runtime: ~5 mins
- Full dataset runtime: ~2 hrs

I have successfully run the test and the full data set on a MacBook Air with 8 CPUs and 16Gb of RAM.

## Knit the RMakdown documents (with .Rmd files)

The 00_main_document_Run1and2.Rmd is the main file referencing all the other children files. To knit
(i.e to generate an html report), run the following commands from an R terminal:

```{R}
library(rmarkdown)
render("./00_main_document_Run1and2.Rmd", output_dir = "html_output")
```

This will generate an standalone html report in the html_output directory.

If loading the .Rmd file from [RStudio](https://posit.co/download/rstudio-desktop/) or [VScode](https://code.visualstudio.com/),
the kniting command will automatically create the document in the final html_output directory.

## Analyzing the full dataset

### Cleaning the slate

Prior to dowloading the full dataset, make sure that the SNV and the data directory are emptied or deleted. Some of the larger files
that are compute intensive are cached in the data directory when ran the first time. There's no mechanism for linking the content of the
SNV directory with the files in the data directory.

From within the pd-sv-analysis repo, run:

```{bash}
rm -rf SNV data
```

### Dowloading the full dataset

The joint genotypes for all cell lines for each chromosome can be found here:

- [study-chr1.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr1.vcf.gz)
- [study-chr2.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr2.vcf.gz)
- [study-chr3.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr3.vcf.gz)
- [study-chr4.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr4.vcf.gz)
- [study-chr5.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr5.vcf.gz)
- [study-chr6.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr6.vcf.gz)
- [study-chr7.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr7.vcf.gz)
- [study-chr8.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr8.vcf.gz)
- [study-chr9.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr9.vcf.gz)
- [study-chry10.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr10.vcf.gz)
- [study-chry11.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr11.vcf.gz)
- [study-chry12.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr12.vcf.gz)
- [study-chry13.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr13.vcf.gz)
- [study-chry14.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr14.vcf.gz)
- [study-chry15.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr15.vcf.gz)
- [study-chry16.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr16.vcf.gz)
- [study-chry17.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr17.vcf.gz)
- [study-chry18.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr18.vcf.gz)
- [study-chry19.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr19.vcf.gz)
- [study-chry20.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr20.vcf.gz)
- [study-chry21.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr21.vcf.gz)
- [study-chry22.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr22.vcf.gz)
- [study-chrX.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chrX.vcf.gz)
- [study-chrY.vcf.gz](https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chrY.vcf.gz)

Here's a simple one liner to download all the files, run in bash from the pd-sv-analysis directory

```{bash}
for x in {1..22} X Y;do wget -P SNV  "https://pd-cell-lines-data.s3.us-west-2.amazonaws.com/joint-genotyping/study-chr${x}.vcf.gz"; done
```

### Analyse the full dataset

The scripts will pick up eachindivuals gunzipped VCF files in the SNV directory and run with them. From an R terminal, simply run the
folowing command and be patient. On a 16 CPU linux machine, it takes ~2 hrs.

```{R}
library(rmarkdown)
render("./00_main_document_Run1and2.Rmd", output_dir = "html_output")
```

## Supporting Files

This analysis uses files specific for the Parkinson Disease CRISPR engineered cell lines. If you would like to run this analysis, the following files
in the supporting_files directory will need to be updated to your design.

- [sgRNA.txt](supporting_files/sgRNA.txt)
- supporting_files/iSCORE-PD_cells_grouped_with_guides.csv

## Creating you own supplemental CAS-Offinder supporing files
