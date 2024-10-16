#! /bin/bash

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

: << 'COMMENT'
# creating index file
genomeDir=/projects/cooper_research2/chinaza/Iron_RNASeq/star
genomeFastaFiles=/users/cnnamdi/rna_seq/pantranscriptome/using_cdna/result/pantranscriptome_for_star.fasta

STAR --runThreadN 24 \
        --runMode genomeGenerate \
        --genomeDir $genomeDir \
        --genomeFastaFiles $genomeFastaFiles \
	--limitGenomeGenerateRAM 43371899829 \
	--genomeSAindexNbases 12
COMMENT


# mapping read to genome
### Loop through all L001_R1 files
### Get a prefix for the output bam, and run star on all 4 files associated with each sample

out=/projects/cooper_research2/chinaza/Iron_RNASeq/aligned_files

for file in /projects/cooper_research2/chinaza/Iron_RNASeq/Trimmed_Fastq/*L001_R1.paired.fastq
do
        filename=`basename $file .paired.fastq`
        prefix=`echo $filename | cut -d'_' -f1`
        file2=`echo $file | sed 's/R1/R2/g'`
        file3=`echo $file | sed 's/L001/L002/g'`
        file4=`echo $file3 | sed 's/R1/R2/g'`

	#### running STAR
	STAR --genomeDir /projects/cooper_research2/chinaza/Iron_RNASeq/star \
		--runThreadN 100 \
        --readFilesIn $file,$file3 $file2,$file4 \
        --outFileNamePrefix $out/$prefix \
        --alignEndsType EndToEnd \
		--alignIntronMax 1 \
		--outSAMunmapped Within \
		--outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard
done

#        echo $filename
#        echo $prefix
#        echo $file
#        echo $file2
#        echo $file3
#        echo $file4


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
