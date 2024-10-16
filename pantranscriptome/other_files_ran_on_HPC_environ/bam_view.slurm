#! /bin/bash

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID : $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

module load samtools/1.11

cd $SLURM_SUBMIT_DIR

echo ""
samtools view Rio-no1Aligned.sortedByCoord.out.bam | head -n6
echo ""

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
