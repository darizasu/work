# GBS-pipeline

This set of wrapper scripts are intended to run the whole GBS pipeline, from raw reads to the final VCF file.
They use [NGSEP](https://sourceforge.net/projects/ngsep/files/Library/) (tested with v.3.0.2), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (tested with v2.2.9), [Picard](http://broadinstitute.github.io/picard/index.html) (tested with v.1.140), [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) (tested with v.0.36) and [htslib](http://www.htslib.org/download/) for [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html). Make sure you have them installed and then provide the location for their executables in these scripts accordingly.
## runPlate.sh
Please follow these steps to run the first stage of the pipeline, called "runPlate":

0) Select your working directory, where all the temporary and final results will be created. From now on (on this small vignette), I will refer to your selected location as `${WD}`. Notice you will have to provide this location in the parameter `WD` in [`runPlate.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/runPlate.sh) and [`runPopulation.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/runPopulation.sh) scripts.
1) Create the following directory and put your raw reads there: 
    `${WD}/reads/lane/`
2) Create an [`INDEXFILE`](https://github.com/darizasu/work/blob/master/GBS-pipeline/INDEXFILE.txt) containing information about sequencing lane, flowcell, barcode and sample-name.
Create a [`FILES2DECONV`](https://github.com/darizasu/work/blob/master/GBS-pipeline/FILES2DECONV.txt) file which contains information about the sequencing lane, flowcell and raw fastq filename(s).
If you have [adapter sequences](https://github.com/darizasu/work/blob/master/GBS-pipeline/adaptersGBS.fa) to be removed from your reads, create a fasta file with them and leave it in `${WD}/reads/`.
Give the location and name of these files in the [`runPlate.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/runPlate.sh) script properly. Take a look at the provided template files in this repository or check [NGSEP - Deconvolute](https://sourceforge.net/projects/ngsep/files/Library/) for more details.
3) Modify the running parameters of [`runPlate.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/runPlate.sh), which are all the lines in the script before the "DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD" warning.
Then run Deconvolution, Trimming and Mapping, it could take several hours:

`./runPlate.sh 'DTM' &`

   Once the mapping step is done, a file called `readPosStats_XX.txt` is created in `${WD}`. It contains the total sequencing error bias along the read position for unique-mapping reads (first column) and multi-mapping reads (second column), calculated for all the samples in the plate. Plot these results along the read position. You have to decide the i5 and i3 parameters by looking at how many consecutive positions in the 5' end and the 3' end have a high sequencing error bias (guide yourself by the behavior of the slope). i5 is to ignore this many base pairs from the 5' end of the reads. i3 is to ignore this many base pairs from the 3' end of the reads, check [NGSEP - FindVariants](https://sourceforge.net/projects/ngsep/files/Library/) manual page. Once you have decided the value for these parameters, specify them in the script [`runPlate.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/runPlate.sh) accordingly.

4) Run Variant discovery, it could take 1 -2h:

`./runPlate.sh 'V' &`

   When Variant Discovery is done, you'll be half way to get the final population VCF file. Check the individual BAM and VCF files for every one of your samples were produced successfully in the directory `${WD}/mapping`

5) Create a [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt) file, which contains information about the BAM and VCF files for every sample you want to include in your population's VCF file. It also contains some of the parameters to call variants using [NGSEP - FindVariants](https://sourceforge.net/projects/ngsep/files/Library/). Its full path must be specified in the parameter `samples2population` either in the [`compareRepeatedSamples.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/compareRepeatedSamples.sh) (optional), and the [`runPopulation.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/runPopulation.sh) scripts, which will be executed in the next stages of the pipeline. This file may contain comment lines starting with '#'. This [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt) file contains 4 tab-separated columns:

        /path/to/sampleID[tab]sample_name[tab]ignore5[tab]ignore3

  * `/path/to/sample` gives the full location + sample-prefix to a sample's VCF and BAM files.
  In other words, you must be able to find these files in the specified directory (check parameters `BAMext` and `VCFext` in the [`runPopulation.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/runPopulation.sh) script for more info):
  /path/to/sample_${VCFext}
  /path/to/sample_${BAMext}
  These files are produced after running the first stage of this pipeline with [`runPlate.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/runPlate.sh) correctly.

  * `sample_name` is the name that every sample will take in the final VCF file.
  Beware to avoid repeated names in this second column. In case of repeated samples, you must name every sample as (e.g.): 
  `sample_name-p01F12`, the '-p' in lowercase is mandatory (MANDATORY) after 'samplename' ('p' stands for 
  'plate'), and p01F12 may indicate the plate number and well that contained that sample.

  * `ignore5` and `ignore3` are to ignore this many base pairs from the 5' and 3' ends in [NGSEP - FindVariants](https://sourceforge.net/projects/ngsep/files/Library/).
 
See the example below:

    /bioinfo2/projects/GBSplates/01/mapping/ALB_213	ALB_213-p01H10	4	10
    
## compareRepeatedSamples.sh
Please follow these steps to run the optional stage of the pipeline, called "compareRepeatedSamples". This stage should be run just in case you have one or more samples that have been sequenced more than once and you want to verify that they are actually the same. First of all, you need a dense list of variants for your species in VCF format. This list can be obtained by running [NGSEP - MergeVariants](https://sourceforge.net/projects/ngsep/files/Library/), check it out for more details.
The samples to be compared in this stage will be genotyped with the provided list of variants.

1) After creating a [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt) file, make sure that the repeated samples you want to compare have the same `sample_name` prefix followed by `-p` in lowercase (MANDATORY) and their corresponding unique identifiers, all of this in the second tab-separated value. Please refer to the 5th step in the previous stage 'runPlate'. You can check the provided template [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt) file. In the following example, the sample BAT_093 was sequenced three different times, and it will be compared to see if all the three samples are exactly the same genotype. Notice the second tab-separated value has the same prefix followed by unique `-pXXX` identifiers:

		#/bionas1/bean/GBSplates/21/mapping/v2.1/BAT_093	BAT_093-p21	6	12
		#/bionas1/bean/GBSplates/23/mapping/v2.1/BAT_093-D02	BAT_093-p23D02	6	12
		#/bionas1/bean/GBSplates/23/mapping/v2.1/BAT_093-E10	BAT_093-p23E10	6	12

The complete file path for [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt) should be specified in the `samples2population` parameter in the [`compareRepeatedSamples.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/compareRepeatedSamples.sh) script.

2) Specify the list of variants for your species in VCF format in the parameter `myVariants` in the script [`compareRepeatedSamples.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/compareRepeatedSamples.sh).

3) Modify the running parameters of [`compareRepeatedSamples.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/compareRepeatedSamples.sh) properly, which are all the lines in the script before the "DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD" warning.
Then run the script from the `${WD}` directory. It could take several hours depending on the number of repeated samples you want to compare and the density of variants provided in the `myVariants` parameter:

`./compareRepeatedSamples.sh &`

Once [`compareRepeatedSamples.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/compareRepeatedSamples.sh) finishes, `${WD}/genotyping` contains a single directory for every repeated sample identified in [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt). Taking back the example above, the directory `${WD}/genotyping/BAT_093` contains the following files:

Symlinks to the original BAM file for every repeated sample identified in [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt).
* `BAT_093-p21_bowtie2_sorted.bam`
* `BAT_093-p23D02_bowtie2_sorted.bam`
* `BAT_093-p23E10_bowtie2_sorted.bam`

A text file with the output of [NGSEP - CompareVCF](https://sourceforge.net/projects/ngsep/files/Library/). Check carefully the details of the comparison, as they will tell you if the tested samples are the same or not.
* `BAT_093_CompareVCF_q60.txt`

A merged BAM file that contains alignments for every BAM file with the prefix `BAT_093` identified in [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt).
* `BAT_093_bowtie2_sorted.bam`

A VCF file with variants discovered from the merged BAM file and its corresponding log file.
* `BAT_093_bowtie2_NGSEP.log`
* `BAT_093_bowtie2_NGSEP.vcf.gz`

In addition, notice that new lines have been added at the bottom part of the file [`samples2population`](https://github.com/darizasu/work/blob/master/GBS-pipeline/samples2population.txt). These new lines contain some information about the newly merged samples, their location, their name and their corresponding [NGSEP - FindVariants](https://sourceforge.net/projects/ngsep/files/Library/) parameters.

Now, based on the analysis you've done for every `sample_name_CompareVCF_q60.txt` file, you decide if you want to include the merged BAM and VCF file in the final population VCF. All you need to do is to comment (with '#') the samples you want to avoid using for the final population VCF. Going back to the example above, we found that all individual entries for BAT_093 were the same. Then, we will comment the lines that mark individual entries an leave uncommented the line with the merged files produced after running [`compareRepeatedSamples.sh`](https://github.com/darizasu/work/blob/master/GBS-pipeline/compareRepeatedSamples.sh):

	#/bionas1/bean/GBSplates/21/mapping/v2.1/BAT_093	BAT_093-p21	6	12
	#/bionas1/bean/GBSplates/23/mapping/v2.1/BAT_093-D02	BAT_093-p23D02	6	12
	#/bionas1/bean/GBSplates/23/mapping/v2.1/BAT_093-E10	BAT_093-p23E10	6	12
	/bioas1/bean/myPopulation/genotyping/BAT_093/BAT_093	BAT_093	6	12



## runPopulation.sh
Please follow these steps to run the second stage of the pipeline, called "runPopulation":

