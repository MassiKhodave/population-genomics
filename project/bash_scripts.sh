#!/bin/bash
 
# running trimmomatic
 
cd /data/project_data/brachy/uncleaned
 
for f1 in *.fastq.gz

do

java -classpath /data/popgen/Trimmomatic-0.39/trimmomatic-0.39.jar org.usadellab.trimmomatic.TrimmomaticSE \
        -threads 10 \
        -phred33 \
         /data/project_data/brachy/uncleaned/"$f1" \
         /data/project_data/brachy/cleaned/"$f1".fastq.gz \
        ILLUMINACLIP:/data/popgen/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:6:20 \
        MINLEN:36 
>> log.txt
 
done

## running Salmon

cd /data/project_data/brachy/cleaned
conda activate salmon

salmon index -t Bdistachyon_556_v3.2.transcript_primaryTranscriptOnly.fa -i brachy_index -p 8
conda activate salmon

for i in $(ls /data/project_data/brachy/cleaned | grep '.fastq.gz' | cut -f 1-3 -d "_"| uniq);
do

    echo "starting sample ${i}"
	read1=$(ls /data/project_data/brachy/cleaned | grep ${i} | grep '.fastq.gz')
	salmon quant -i /data/project_data/brachy/cleaned/brachy_index \
        -l A \
         -r /data/project_data/brachy/cleaned/${read1} \
		    -p 8  \
         --softclip \
         --seqBias \
         --gcBias \
         -o /data/project_data/brachy/transcripts_quant/${i}

    echo "sample ${i} done"

done

