# GBS-pipeline

This set of wrapper scripts are intended to run the whole GBS pipeline, from raw reads to the final VCF file.
They use [NGSEP](https://sourceforge.net/projects/ngsep/files/Library/) (tested with v.3.0.2), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (tested with v2.2.9), [Picard](http://broadinstitute.github.io/picard/index.html) (tested with v.1.140), [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) (tested with v.0.36) and [htslib](http://www.htslib.org/download/) for [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html). Make sure you have them installed and then provide the location for their executables in these scripts accordingly.
Please follow these steps:

1) Create the following directory and put your raw reads there: 
    `${WD}/reads/lane/`
2) Create an `INDEXFILE` containing information about sequencing lane, flowcell, barcode and sample-name. In addition, create a `FILES2DECONV` which contains information about the sequencing lane, flowcell and raw fastq filename(s). Give the location and name of these files in the `runPlate.sh` script properly. Take a look at the provided template files in this repository or check NGSEP-Deconvolute for more details.
3) Modify the running parameters of `runPlate.sh`. Then run Deconvolution, Trimming and Mapping:

`./runPlate.sh 'DTM'`

4) Run `calculateReadPos.sh` from the directory `${WD}/mapping`
