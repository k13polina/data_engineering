ctrl1=/disk1/COURSE/Data/Homework/Reads/RNA-seq/Student_13/st13_ctrl_1.fastq
ctrl2=/disk1/COURSE/Data/Homework/Reads/RNA-seq/Student_13/st13_ctrl_2.fastq
exp1=/disk1/COURSE/Data/Homework/Reads/RNA-seq/Student_13/st13_exp_1.fastq
exp2=/disk1/COURSE/Data/Homework/Reads/RNA-seq/Student_13/st13_exp_2.fastq

athaliana_genome=/disk1/COURSE/Data/2-trimming_practice/references/Athaliana.fasta
athaliana_annotation=/disk1/COURSE/Data/2-trimming_practice/references/Athaliana.gtf

/disk1/COURSE/Tools/STAR-2.7.10a/bin/Linux_x86_64/STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles $athaliana_genome \
--sjdbGTFfile $athaliana_annotation \
--sjdbOverhang 59

/disk1/COURSE/Tools/STAR-2.7.10a/bin/Linux_x86_64/STAR \
--runThreadN 10 \
--genomeDir ./ \
--readFilesIn $exp1 \
--outFileNamePrefix STAR_mapping_exp1

samtools view -c STAR_mapping_exp2Aligned.out.sam
samtools view -c -F 0x4 STAR_mapping_exp2Aligned.out.sam
samtools view -c -f 0x800 STAR_mapping_exp2Aligned.out.sam

samtools view -F 0x4 -F 0x800 STAR_mapping_exp2Aligned.out.sam -b -h -@ 10 | samtools sort -@ 10 - > STAR_mapping_exp2.sort.bam && samtools index STAR_mapping_exp2.sort.bam
samtools coverage STAR_mapping_ctrl2.sort.bam -> STAR_mappin.txt
samtools depth STAR_mapping_ctrl1.sort.bam -> depth_ctrl1.txt

samtools depth -b /disk1/COURSE/Data/2-trimming_practice/references/Athaliana.gff STAR_mapping_ctrl1.sort.bam > depth.txt

htseq-count -f sam --stranded no --nonunique none STAR_mapping_ctrl1Aligned.out.sam /disk1/COURSE/Data/2-trimming_practice/references/Athaliana.gtf > output_ctrl1.txt