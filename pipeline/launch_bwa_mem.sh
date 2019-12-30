
module load bwa/0.7.12
module load samtools/1.3

OUTPUT_DIR=/nl/umw_job_dekker/users/fc85w/pacbio_aligner/bwa_outputs

DIGESTED_READS_DIR=/nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads

for f in $DIGESTED_READS_DIR/*fastq; do
	echo "$f"
	dataset=$(echo "$f" | cut -d/ -f 8 | cut -d. -f 1)
	echo "$dataset.bam"
	bsub -q short -n 32 -R "rusage[mem=32000]" -W 240 \
		"bwa mem -t 32 $OUTPUT_DIR/hg19.fa $f | samtools view -bS - > $OUTPUT_DIR/$dataset.bam
         	 samtools sort -n $OUTPUT_DIR/$dataset.bam > $OUTPUT_DIR/$dataset.sorted.bam
		 samtools flagstat $OUTPUT_DIR/$dataset.sorted.bam > $OUTPUT_DIR/$dataset.flagstat_out"
done


