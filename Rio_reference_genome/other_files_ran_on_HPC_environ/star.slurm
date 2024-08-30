#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=star
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --mail-user=cnnamdi@charlotte.edu
#SBATCH --mail-type=ALL

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID : $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

module load star

cd $SLURM_SUBMIT_DIR

### Use the previously created index file for the Rio Reference Genome
genomeDir=/projects/cooper_research/Ref_Genomes/Rio/version_2.1/star

### Define a shortcut variable for your output folder
out=/projects/cooper_research2/chinaza/Iron_RNASeq/rio-ref-genome/aligned_files

### Loop through all L001_R1 files in the raw data directory
### Get a prefix for the output bam, and run STAR on all 4 files associated with each sample
for file in /projects/cooper_research2/chinaza/Iron_RNASeq/Trimmed_Fastq/*L001_R1.paired.fastq
do
        filename=`basename $file .paired.fastq`
        prefix=`echo $filename | cut -d'_' -f1`
        file2=`echo $file | sed 's/R1/R2/g'`
        file3=`echo $file | sed 's/L001/L002/g'`
        file4=`echo $file3 | sed 's/R1/R2/g'`

        #### Run the STAR aligner
        STAR --genomeDir $genomeDir \
                --runThreadN 16 \
                --readFilesIn $file,$file3 $file2,$file4 \
                --outFileNamePrefix $out/$prefix \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard
done


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
