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

out=/projects/cooper_research2/chinaza/Iron_RNASeq/aligned_files

for file in /projects/cooper_research2/chinaza/Iron_RNASeq/aligned_files/*.bam
do

#	echo $file
	filename=`basename $file .bam`
#	echo $out/$filename

	samtools view --threads 24 \
	-h $file -o $out/$filename.sam
done

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
