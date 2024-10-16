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
gtf=/users/cnnamdi/rna_seq/pantranscriptome/using_cdna/result/pseudo-annotation.gtf

### Define the path to the alignment files
sam_files=/projects/cooper_research2/chinaza/Iron_RNASeq/aligned_files/*.sam

: << 'COMMENT'
#double_checking files
for file in $sam_files
do
       echo $file
done
COMMENT

### Run the counter on all files at once to get a matrix
$multimap/mmquant -a $gtf -r $sam_files \
        -o rawCounts20_pantranscriptome.txt \
        -n 069-con1 069-con2 069-no1 069-n02 \
        972-con1 972-con2 972-no1 972-no2 \
        Gra-con1 Gra-con2 Gra-con3 Gra-no1 Gra-no2 Gra-no3 \
        Leo-con1 Leo-con2 Leo-no1 Leo-no2 Rio-con1 Rio-no1 \
        -t 20 -e N -f SAM -p

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
