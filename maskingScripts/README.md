This directory contains a brief description of the steps one could take to mask variant calls prior to calculating summary statistics for input to FILET.

## Files included in this directory:

makeMaskArmFilesFromQualAndRM.py
makeFilteredSnpVcfsForEachArm.py
phasedVcfsToFastas.py
miscFuncs.py

There are certain parts in the first three files where there is hard-coded information relevant only to the simulans-sechellia analysis conducted in the FILET paper. I have tried to mark these places with #TODO comments explaining how these would need to be modified for a different data set. The fourth file contains a few utility functions that you shouldn't have to touch.

As is the case for the main FILET pipeline, all python scripts were written for Python v2 and would need some tweaking to get working with Python3.

## Masking steps for real data

1) `python makeMaskArmFilesFromQualAndRM.py gVcfFileName1 gVcfFileName2 repeatGffFileName refFaFileName qualThresh sampSizePerc maskedRefFaFileName`

Arguments: `gVcfFileName1` and `gVcfFileName2` are GVCF files for populations 1 and 2, respectively; `repeatGffFileName` is a GFF file containing the locations of annotated repetitive elements (e.g. via RepeatMasker); `refFaFileName` is the reference genome in fasta format; `qualThresh` is our genotype quality cutoff (20 in the FILET paper), and `sampSizePerc` is the fraction of all samples that must survive this cutoff for the site (whether polymorphic or monomorphic) to be retained in the analysis; `maskedRefFaFileName`, which is written by this script, is a modified version of the reference genome in fasta format which will contain an 'N' for each site either residing within repetitive sequence or failing to pass our quality cutoff.

2) `python makeFilteredSnpVcfsForEachArm.py vcfFileName1 vcfFileName2 gVcfFileName1 gVcfFileName2 maskedRefFaFileName qualThresh filteredVcfPathPrefix_`

Arguments: `vcfFileName1` and `vcfFileName2` are VCF files for populations 1 and 2 (with variant sites only), while `gVcfFileName1` and `gVcfFileName2` are the same GVCFs used as input in the previous step; `maskedRefFaFileName` is the masked reference fasta file produced by the previous step; `qualThresh` is the quality threshold used to mask individual genotypes (whereas the previous step was masking entire sites/columns of the alignment); `filteredVcfPathPrefix` is the beginning of the path to the output files that will be created by this step--one for each chromosome or chromosome arm in the reference genome--and each file name will be of the form `filteredVcfPathPrefix_arm.vcf` where "arm" is the name of the chromosome or chromosome arm.

3) Step three is to phase the vcfs produced by the previous step (e.g. using shapeit). This is done separately for each chromosome arm, and the output for this step should be a phased vcf file for each chrom/arm.

4) `phasedVcfsToFastas.py vcfFileName maskedRefFaFileName fastaFileName1 fastaFileName2`

This step is done separately for each chromosome arm. vcfFileName is the output for a given arm from step 3. `maskedRefFaFileName` is once again our masked reference genome. `fastaFileName1` and `fastaFileName1` are our resulting fasta inputs for `pgStatsBedSubpop_forML` (see `examplePipeline.sh`).

At this point we are now ready to proceed with calculating FILET's feature vector, as described in Step 5 of examplePipeline.sh. One thing to keep in mind is that that step will require a separate reference genome file for each chromosome arm (while the scripts here use one fasta file containing the whole reference genome).

## Masking steps for simulated data

Once you have completed the above steps, you also have all the necessary information to mask your simulated training/test data if desired. Here is the pipeline:

1) `python makeMaskFilesFromMaskedRef.py maskedRefFaFileName maskFileDir winSize`

This step uses the masking information from our real data to write out masking files to apply to our simulated data. Here, `maskedRefFileName` is the name of our masked reference genome created by the above pipeline, `maskFileDir` is a directory where the masking files will be written during this step, and `winSize` is the size of windows that will be fed as input to FILET, given in bp. Note that this step may create a large number of files, one for each window in the genome.

2) `cat simulationInputFile | msMaskAllRows maskFileName | python removeNedOutColumnsFromMsFile.py stdin > outputFile"

This step will take one of the mask files created by the previous step and applies the masking to a single simulation replicate. Here `simulationInputFile` is the path your .ms-style formatted simulation input, and `maskFileName` is the path to your mask file (one of many created in the previous step). The masked version will be written to `outputFile`. To run this command on ms-style data with multiple replicates, the mask file will also have to have data for the corresponding number of regions. I.e. if `simulationInputFile` has 100 reps, then you will want `maskFileName` to have data from 100 windows from the masked reference genome---this can be achieved by simply concatenating 100 randomly selected regions together (ideally in a random order).
