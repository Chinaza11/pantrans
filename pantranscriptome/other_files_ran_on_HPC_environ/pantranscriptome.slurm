#! /bin/bash

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID : $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

module load anaconda3/2022.10
#module load muscle/3.8.31
#module load clustal-omega/1.2.4

cd $SLURM_SUBMIT_DIR

#python making_pav_table_for_cluster.py

#python sequence_alignment.py

#python making_pantranscriptome.py

#python making_pantranscriptome_2.py

#python making_pseudo-gff.py

python making_pantranscriptome_for_star.py

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
