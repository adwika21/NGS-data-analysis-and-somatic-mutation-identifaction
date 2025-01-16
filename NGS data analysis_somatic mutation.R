fastqc PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz \
PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz \
PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz \
PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz \
-o qc_reports/

multiqc qc_reports/ -o qc_summary/

bowtie2-build reference.fasta reference_index
# Align paired-end reads using Bowtie2
bowtie2 -x reference_index \
        -1 PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz \
        -2 PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz \
        -S tumor_aligned.sam

bowtie2 -x reference_index \
        -1 PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz \
        -2 PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz \
        -S normal_aligned.sam

# Convert SAM to BAM
samtools view -bS tumor_aligned.sam > tumor_aligned.bam
samtools view -bS normal_aligned.sam > normal_aligned.bam

# Sort the BAM files
samtools sort tumor_aligned.bam -o tumor_aligned_sorted.bam
samtools sort normal_aligned.bam -o normal_aligned_sorted.bam

# Index the sorted BAM files
samtools index tumor_aligned_sorted.bam
samtools index normal_aligned_sorted.bam


# Generate mpileup for tumor and normal BAM files
samtools mpileup -f reference.fasta tumor_dedup.bam > tumor.mpileup
samtools mpileup -f reference.fasta normal_dedup.bam > normal.mpileup

# Run VarScan to detect somatic mutations
java -jar VarScan.v2.4.3.jar somatic tumor.mpileup normal.mpileup somatic_output --min-var-freq 0.05 --p-value 0.01 --output-vcf 1

# View the VCF results
cat somatic_output.Somatic.VCF

# filter variants further with bcftools
bcftools view somatic_output.Somatic.VCF -i 'QUAL>20' > filtered_somatic.vcf

