# GBS-pipeline

This set of wrapper scripts are intended to run the whole GBS pipeline, from raw reads to the final VCF file.
They use [NGSEP](https://sourceforge.net/projects/ngsep/files/Library/) (tested with v.3.0.2), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (tested with v2.2.9), [Picard](http://broadinstitute.github.io/picard/index.html) (tested with v.1.140), [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) (tested with v.0.36) and [htslib](http://www.htslib.org/download/) for [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html). Make sure you have them installed and then provide the location for their executables in these scripts accordingly.
Please follow these steps:

0) Select your working directory, where all the temporary and final results will be created. From now on (on this small vignette), I will refer to your selected location as `${WD}`. Notice you will have to provide this location in the parameter "WD" in `runPlate.sh` and `runPopulation.sh` scripts.
1) Create the following directory and put your raw reads there: 
    `${WD}/reads/lane/`
2) Create an `INDEXFILE` containing information about sequencing lane, flowcell, barcode and sample-name.
Create a `FILES2DECONV` file which contains information about the sequencing lane, flowcell and raw fastq filename(s).
If you have adapter sequences to be removed from your reads, put them in a fasta file in `${WD}/reads/`.
Give the location and name of these files in the `runPlate.sh` script properly. Take a look at the provided template files in this repository or check NGSEP-Deconvolute for more details.
3) Modify the running parameters of `runPlate.sh`, which are all the lines in the script before the "DO NOT MODIFY ANYTHING FROM THIS POINT FORWARD" warning.
Then run Deconvolution, Trimming and Mapping, it could take several hours:

`./runPlate.sh 'DTM' &`

4) Run the script `calculateReadPos.sh` providing the `${WD}` and "plateID" information, it shouldn't take more than 10 seconds:

`calculateReadPos.sh ${WD} myPlateX`

After running it, a file called "readPosStats_myPlateX.txt" is created in `${WD}`. It contains the total sequencing error bias along the read position for unique-mapping reads (first column) and multi-mapping reads (second column), calculated for all the samples in the plate. Plot these results along the read position. You have to decide the i5 and i3 parameters by looking at how many consecutive positions in the 5' end and the 3' end have a high sequencing error bias (guide yourself by the behavior of the slope). i5 is to ignore this many base pairs from the 5' end of the reads. i3 is to ignore this many base pairs from the 3' end of the reads, check NGSEP-FindVariants manual page. Once you have decided the value for these parameters, specify them in the script `runPlate.sh` accordingly.

5) Run Variant discovery, it could take 1 -2h:

`./runPlate.sh 'V' &`

6) That's all folks ! Just kidding, you're half way to get the final population VCF file. Check the individual BAM and VCF files for every one of your samples were produced successfully in the directory `${WD}/mapping`
