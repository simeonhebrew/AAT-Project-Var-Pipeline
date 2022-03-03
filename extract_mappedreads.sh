
#Indexing the reference genome
bwa index TriTrypDB-52_TcongolenseIL3000_Genome.fasta

#Mapping reads against the reference genome
bwa mem -t 8 -M TriTrypDB-52_TcongolenseIL3000_Genome.fasta ${Read1} ${Read2} > file.bam

#Retrieving mapping statistics
samtools flagstat file.bam

#Extracting mapped reads
samtools view -b -F 4 file.bam > mapped.bam

#Sorting mapped reads
samtools sort -n -o sortedaligned.bam aligned.bam

#Converting mapped reads from bam to fastq format
bedtools bamtofastq -i sortedaligned.bam \
                      -fq aln.end1.fq \
                      -fq2 aln.end2.fq
