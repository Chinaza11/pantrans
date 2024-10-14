#! /bin/bash

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID : $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""


### Set a shortcut to the program
multimap=/projects/cooper_research/Programs/multi-mapping-counter

### Set a shortcut to the gtf annotation file
gtf=/projects/cooper_research/Ref_Genomes/Rio/version_2.1/annotation/SbicolorRio_468_v2.1.gene.gtf

### Define the path to the BAM files
src=/projects/cooper_research2/chinaza/Iron_RNASeq/rio-ref-genome/aligned_files/*sortedByCoord.out.bam

### Run the counter on all files at once to get a matrix
$multimap/mmquant -a $gtf -r $src -o rawCounts20.txt \
	-n 069-con1 069-con2 069-no1 069-n02 \
	972-con1 972-con2 972-no1 972-no2 \
	Gra-con1 Gra-con2 Gra-con3 Gra-no1 Gra-no2 Gra-no3 \
	Leo-con1 Leo-con2 Leo-no1 Leo-no2 Rio-con1 Rio-no1 \
	-f BAM

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
