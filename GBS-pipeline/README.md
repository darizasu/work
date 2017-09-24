# GBS-pipeline

This set of wrapper scripts are intended to run the whole GBS pipeline, from raw reads to the final VCF file.
They use [NGSEP](https://sourceforge.net/projects/ngsep/files/Library/) (tested with v.3.0.2), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (tested with v2.2.9), [Picard](http://broadinstitute.github.io/picard/index.html) (tested with v.1.140), [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) (tested with v.0.36) and [htslib](http://www.htslib.org/download/) for [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html). Make sure you have them installed and then provide the location for their executables in these scripts accordingly.
## runPlate.sh
Please follow these steps to run the first stage of the pipeline, called "runPlate":

0) Select your working directory, where all the temporary and final results will be created. From now on (on this small vignette), I will refer to your selected location as `${WD}`. Notice you will have to provide this location in the parameter "WD" in `runPlate.sh` and `runPopulation.sh` scripts.
1) Create the following directory and put your raw reads there: 
    `${WD}/reads/lane/`
2) Create an [`INDEXFILE`](https://github.com/darizasu/work/blob/master/GBS-pipeline/INDEXFILE.txt) containing information about sequencing lane, flowcell, barcode and sample-name.
Create a [`FILES2DECONV`](https://github.com/darizasu/work/blob/master/GBS-pipeline/FILES2DECONV.txt) file which contains information about the sequencing lane, flowcell and raw fastq filename(s).
If you have [adapter sequences](https://github.com/darizasu/work/blob/master/GBS-pipeline/adaptersGBS.fa) to be removed from your reads, put them in a fasta file in `${WD}/reads/`.
Give the location and name of these files in the `runPlate.sh` script properly. Take a look at the provided template files in this repository or check [NGSEP - Deconvolute](https://sourceforge.net/projects/ngsep/files/Library/) for more details.
3) Modify the running parameters of `runPlate.sh`, which are all the lines in the script before the "DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD" warning.
Then run Deconvolution, Trimming and Mapping, it could take several hours:

`./runPlate.sh 'DTM' &`

   Once the mapping step is done, a file called `readPosStats_XX.txt` is created in `${WD}`. It contains the total sequencing error bias along the read position for unique-mapping reads (first column) and multi-mapping reads (second column), calculated for all the samples in the plate. Plot these results along the read position. You have to decide the i5 and i3 parameters by looking at how many consecutive positions in the 5' end and the 3' end have a high sequencing error bias (guide yourself by the behavior of the slope). i5 is to ignore this many base pairs from the 5' end of the reads. i3 is to ignore this many base pairs from the 3' end of the reads, check [NGSEP - FindVariants](https://sourceforge.net/projects/ngsep/files/Library/) manual page. Once you have decided the value for these parameters, specify them in the script `runPlate.sh` accordingly.

4) Run Variant discovery, it could take 1 -2h:

`./runPlate.sh 'V' &`

   When Variant Discovery is done, you'll be half way to get the final population VCF file. Check the individual BAM and VCF files for every one of your samples were produced successfully in the directory `${WD}/mapping`

5) Create a [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt) file, which contains information about the BAM and VCF files for every sample you want to include in your population's VCF file. It also contains some of the parameters to call variants using [NGSEP - FindVariants](https://sourceforge.net/projects/ngsep/files/Library/). Its full path must be specified in the parameter `samples2population` either in the [`compareRepeatedSamples.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/compareRepeatedSamples.sh) (optional), and the [`runPopulation.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/runPopulation.sh) scripts, which will be executed in the next stages of the pipeline. This file may contain comment lines starting with '#'. This [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt) file contains 4 tab-separated columns:

        /path/to/sampleID[tab]sample_name[tab]ignore5[tab]ignore3

  * `/path/to/sample` gives the full location + sample-prefix to a sample's VCF and BAM files.
  In other words, you must be able to find these files in the specified directory (check parameters `BAMext` and `VCFext` in the `runPopulation.sh` script for more info):
  /path/to/sample_${VCFext}
  /path/to/sample_${BAMext}
  These files are produced after running the first stage of this pipeline with `runPlate.sh` correctly.

  * `sample_name` is the name that every sample will take in the final VCF file.
  Beware to avoid repeated names in this second column. In case of repeated samples, you must name every sample as (e.g.): 
  `sample_name-p01F12`, the '-p' in lowercase is mandatory (MANDATORY) after 'samplename' ('p' stands for 
  'plate'), and p01F12 may indicate the plate number and well that contained that sample.

  * ignore5 and ignore3 are to ignore this many base pairs from the 5' and 3' ends in [NGSEP - FindVariants](https://sourceforge.net/projects/ngsep/files/Library/).
 
See the example below:

    /bioinfo2/projects/GBSplates/01/mapping/ALB_213	ALB_213-p01H10	4	10
    
## compareRepeatedSamples.sh
Please follow these steps to run the optional stage of the pipeline, called "compareRepeatedSamples". This stage should be run just in case you have one or more samples that have been sequenced more than once and you want to verify that they are actually the same:
## runPopulation.sh
Please follow these steps to run the second stage of the pipeline, called "runPopulation":

