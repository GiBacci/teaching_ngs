# reads quality check using FastQC
fastqc B01_aligned_1.fastq.gz

# check the first 4 lines of fastq file
zcat B01_aligned_1.fastq.gz | head -n 4

# count the number of reads in the fastq files
zcat B01_aligned_1.fastq.gz | awk 'END {print NR/4}'
zcat B01_aligned_2.fastq.gz | awk 'END {print NR/4}'

# calculate the average reads length for each file
zcat B01_aligned_1.fastq.gz | awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}'
zcat B01_aligned_2.fastq.gz | awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}'

# calculate the length of the reference genome
awk '!/^>/{seq=seq$1}END{print length(seq)}' assembly.fasta

# for loops for mapping the reads to the reference genome
for i in B01 B02; do bowtie2 -x assembly -1 ${i}_aligned_1.fastq.gz -2 ${i}_aligned_2.fastq.gz --very-fast -S ${i}.sam; done

# displaying the first five lines of a sam file
head -n 5 B01.sam

# Removal of PCR duplicates
for i in B01 B02; do samtools sort -o ${i}_sorted.sam ${i}.sam; done
for i in B01 B02; do samtools rmdup ${i}_sorted.sam ${i}_rmdup.sam; done

# Generation of the coverage table with bedtools
for i in B01 B02;do samtools view -b -S ${i}_rmdup.sam | bedtools genomecov -pc -ibam stdin > ${i}.coverage; done

#visualize the first 5 rows of a GFF file
head -n 5 assembly.gff

#remove non-CDS lines from a GFF file
awk -F'\t' '{if($3 == "CDS") print $0}' assembly.gff > assembly_cds.gff

#intersect bam and gff files
for i in B01 B02;do samtools view -b -S ${i}_rmdup.sam | bedtools intersect -c -bed -a stdin -b assembly_cds.gff > ${i}_intersect.bed; done

#show the content of the bed file
head -n 3 B01_intersect.bed

#count the reads that map outside the coding region
cut -f13 B01_intersect.bed | sort | uniq -c | sort -nr

#for loop for counting the reads that map outside the coding region. 
for i in B01 B02;do echo ${i}; cut -f13 ${i}_intersect.bed | sort | uniq -c | sort -nr; done

#using awk to generate a bed file with non ambiguously mapping reads
awk '{if($NF == 1) print $0}' B01_intersect.bed > B01_unambiguous.bed

#embedding the awk command to generate a bed file with non ambiguously mapped reads in a for loop
for i in B01 B02;do awk '{if($NF == 1) print $0}' ${i}_intersect.bed > ${i}_unambiguous.bed; done

#generating non-ambiguous reads mapping bed files
for i in B01 B02;do bedtools coverage -a assembly_cds.gff -b ${i}_unambiguous.bed > ${i}_genecov.bed; done

