# WGS-pipeline
This set of scripts are intended to run WGS pipeline from raw reads to the final VCF file
## runSamplesWGS.sh
This wrapper script is intended to run a test to get Insert Length parameters (using CollectInsertSizeMetrics from Picard) Trimming (using Trimmomatic), Mapping (using Bowtie2) and FindVariants (using NGSEP) for a list of WGS samples.
Please follow these instructions to run this script properly.
First of all, you need to have installed [NGSEP](https://sourceforge.net/projects/ngsep/files/Library/) (tested with v.3.0.2), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (tested with v2.2.9), [Picard](http://broadinstitute.github.io/picard/index.html) (tested with v.1.140) and [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) (tested with v.0.36).

The variable `samples2test` specifies the full path and filename to the file that contains parameters for Mapping and FindVariants.
From the begining, this file is created, used, then modified the 1st time, used, then modified the 2nd time, and used for the last time.
In the 1st step (see below), `samples2test` contains one line for each sample_id
In the 2nd step (see below), `samples2test` contains three tab-separated columns: sample_id[tab]minins[tab]maxins.
In the 3rd step (see below), `samples2test` contains five tab-separated columns: sample_id[tab]minins[tab]maxins[tab]I5prime[tab]I3prime
Please follow these steps.

  1) Create a txt file containing a sample_id per line. Comment lines with '#' are accepted.
  Specify the name of this file in the `samples2test` variable as described above.
  Execute the 1st run of this script:
  
  `./runMapping_WGS.sh 'I'` (The 'I' runs the InsertLength test only). 
  
  It takes the first 1.000.000 reads for each sample's reads and maps them to the reference genome.
  Then it runs Picard's CollectInsertSizeMetrics module to get the insert size characteristics for each sample's sequencing run. After running this module, a pdf histogram and a txt files are produced showing the insert size distribution. You have to decide the minimum and maximum insert size length for each sample by checking the points where the slope starts to increase and decrease rapidly (minimum and maximum insert size length respectively).
  minins is the minimum fragment length for a valid paired-end alignment.
  maxins is the maximum fragment length for a valid paired-end alignment, check Bowtie2 manual page.
  Once you have decided the minins and maxins parameters, specify them in `samples2test` file for each sample_id as the second and third tab-separated columns respectively.

  2) Execute the 2nd run of this script: 
  
  `./runMapping_WGS.sh 'TM'` (The 'TM' runs trimming and mapping only). 
  
 Once correctly completed, the files `*_bowtie2_readpos.stats` are produced by using NGSEP-QualStats. Plot these results. You have to decide the I5prime and I3prime parameters by looking at how many consecutive positions in the 5' end and the 3' end have a high sequencing error bias (guide yourself by the behavior of the slope).
  I5prime is to ignore this many base pairs from the 5' end of the reads.
  I3prime is to ignore this many base pairs from the 3' end of the reads, check NGSEP manual page.
  Once you have decided the I5prime and I3prime parameters, specify them in 'samples2test' for each sample_id as the fourth and fifth tab-separated columns respectively.

  3) Exectute the 3rd run of this script:
  
  `./runMapping_WGS.sh 'V'` (The 'V' runs VariantDiscovery only).
