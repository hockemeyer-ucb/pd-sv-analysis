# iSCORE-PD SNV analysis

## Pre-requesite

### Software

- [R 4.5](https://cran.r-project.org/)
- [Java Runtime (JRE) 21.0.1](https://www.oracle.com/java/technologies/downloads/?er=221886)
- [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- [wget](https://www.gnu.org/software/wget/) or [curl](https://curl.se/)
- C compiler with make if running from Linux

### Additional libraries if installing from vanilla Ubuntu 24.04 LTS

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
options(Ncpus = parallel::detectCores())
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

The script will read the VCF files in parallel up to the max number of CPU and will consume copius amount of RAM during the serialization process. It has been successully develop and tested using the following parameters. Please note that in it's current version, the script will crash an instance with 8 CPU and 64 Gb of RAM.

- OS: Ubuntu 24.04 LTS
- Architecture: ARM64
- Diskspace: 25 Gb
- CPU: 16
- RAM: 124 Gb
- Test dataset runtime: ~5 mins
- Full dataset runtime: ~2 hrs

I have successfully run the test data on a MacBook Air with Apple M3 AMD64 chipset with 8 CPUs and 24Gb of RAM in less than 5 minutes.

## Kniting the RMarkdown documents (with .Rmd files)

The 00_main_document_Run1and2.Rmd is the main file referencing all the other children files. To knit
(i.e to generate an html report), run the following commands from an R terminal:

```{R}
library(rmarkdown)
render("./00_main_document_Run1and2.Rmd", output_dir = "html_output")
```

This will generate an standalone html report in the html_output directory.

If loading the .Rmd file from [RStudio](https://posit.co/download/rstudio-desktop/) or [VScode](https://code.visualstudio.com/),
their kniting command will automatically create the document in the final html_output directory.

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

The scripts will pick up each indivuals gunzipped VCF files in the SNV directory and run with them. From an R terminal, simply run the
folowing command and be patient. On a 16 CPU linux machine, it takes ~2 hrs.

```{R}
library(rmarkdown)
render("./00_main_document_Run1and2.Rmd", output_dir = "html_output")
```

## Supporting Files

This analysis uses files specific for the Parkinson Disease CRISPR engineered cell lines. If you would like to run this analysis, the following files
in the supporting_files directory will need to be updated to your design.

#### [iSCORE-PD_design.csv](supporting_files/iSCORE-PD_design.csv)

the iSCORE-PD_design.csv file is a comma seperated text file starting with a header and with one cell line per line with the following column header:

the iSCORE-PD_design.csv

#### [iSCORE-PD_cells_grouped_by_editing_methods.csv](supporting_files/iSCORE-PD_cells_grouped_by_editing_methods.csv)

The iSCORE-PD_cells_grouped_with_guides.csv file is a comma seperated text file starting with a header and with one cell line per line with the following
column header:

- samples: Sample ID found in VCF header of the joint genotyping output file
- group: Group ID linking the sample in group
- meta: Additional group relation. Not used anymore

It is used to establish the sample in the VCF to the different cell line edited group. In our analysis, each CRISPR edit has serveral cell line clones.

#### [iSCORE-PD_cells_grouped_with_guides.csv](supporting_files/iSCORE-PD_cells_grouped_with_guides.csv)

The iSCORE-PD_cells_grouped_by_editing_methods.csv is an extention of the previous file with extra columns with the Id of the RNA guide(s)
used to edit each cell lines. The two files could be consolidated but where kept seperate for the ability to change the samples analyzed
in one or the other files. The extra columns represent the maximal number of guides used in any of the cell lines. In our case, some edits
used up to 3 guides, so we add 3 extra columns and for each lines with list the Id of the guides used (either one, two or three guide IDs)

- samples: Sample ID
- group: Group ID
- meta: Additional group relation
- editing_group: Type of Cas use for CRISPR (Cas9, TALEN, PE)
- guide1
- guide2
- guide3

#### [sgRNA.txt](supporting_files/sgRNA.txt)

This file is used by (Cas-OFFinder)[https://github.com/snugel/cas-offinder] to predict the putative location of CRISPR Off-target sites.
It starts with the location of the chromsome FASTA files for the [human Hg38 chromosomes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
on your local hard drive. The second line correspond to the type of RNA used for the search with the size of the RNA and DNA bulge
(Details can be found with the (Cas-OFFinder)[https://github.com/snugel/cas-offinder] documentation). Then the following lines are the RNA sequences for each guide
with the number of missmatch to be considered.

For examples:

```
/home/ubuntu/working/genomes/hg38.fullAnalysisSet.chroms
NNNNNNNNNNNNNNNNNNNNNRG
GGAGGGAGTGGTGCATGGTGNNN	5	SNCA_A53T_peg
TCATAGGAATCTTGAATACTNNN	5	SNCA_A53T_nc
CAGGGTGTGGCAGAAGCAGCNNN	5	SNCA_A30P_peg
```

#### [cas-offinder-out.txt](supporting_files/cas-offinder-out.txt)

The cas-offinder-out.txt file is the output of running (Cas-OFFinder)[https://github.com/snugel/cas-offinder]. It was run once using
[sgRNA.txt](supporting_files/sgRNA.txt) as input like this on a GPU instance on AWS:

```{bash}
cas-offinder sgRNA.txt G cas-offinder-out.txt
```
