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

- CPU: 16
- RAM: 32 Gb
- Test Run Time: ~5 mins
- Full dataset runtime: ~3 hrs

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
