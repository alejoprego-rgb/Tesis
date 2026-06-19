#for x in 1 2 3; do
  for y in 0 1 5 10; do
    cutadapt -j 10 \
      -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
      -A GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG \
      -q 20,20 \
      -m 36 \
      -o /u01/home/galonso/RNAseq_Alejo/trimmed/A${x}D${y}_1_trimmed.fastq.gz \
      -p /u01/home/galonso/RNAseq_Alejo/trimmed/A${x}D${y}_2_trimmed.fastq.gz \
      /u01/home/galonso/RNAseq_Alejo/datos_crudos/Amastigotes_axenicos/A${x}D${y}_1.fastq.gz \
      /u01/home/galonso/RNAseq_Alejo/datos_crudos/Amastigotes_axenicos/A${x}D${y}_2.fastq.gz
  done
done

for x in $(seq 1 4)
do 
bowtie -S -v 2 -p 4 --best -X 1000 -I 0 ~/RNAseq_Alejo/Genomas/dual-seq/dual-seq_index -1 ~/RNAseq_Alejo/trimmed/AI${x}_1_trimmed.fastq.gz -2 ~/RNAseq_Alejo/trimmed/AI${x}_2_trimmed.fastq.gz > ~/RNAseq_Alejo/alineamiento/amas/AI${x}-bowtie.sam && \
samtools view -bS ~/RNAseq_Alejo/alineamiento/amas/AI${x}-bowtie.sam > ~/RNAseq_Alejo/alineamiento/amas/AI${x}-bowtie.bam && \
samtools sort ~/RNAseq_Alejo/alineamiento/amas/AI${x}-bowtie.bam -o ~/RNAseq_Alejo/alineamiento/amas/AI${x}-bowtie_sorted.bam && \
samtools index ~/RNAseq_Alejo/alineamiento/amas/AI${x}-bowtie_sorted.bam
done

for x in $(seq 1 3)
do
featureCounts -t mRNA -g ID -a ~/RNAseq_Alejo/Genomas/TcDm28cT2T_manualCurated.gff -o ~/RNAseq_Alejo/conteos/epis-bowtie-counts${x}.txt ~/RNAseq_Alejo/alineamiento/epis/E${x}-bowtie_sorted.bam
done
