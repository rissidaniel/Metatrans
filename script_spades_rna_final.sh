############ rename fastq reads first

        echo "starting"

# for metagenomics remove the --nobins option from SqueezeMeta and SAmple_file.txt
# change all the paths for the current folder

mkdir raw
cp *fastq.* raw

mkdir raw/repaired
mkdir raw/repaired/fastqc

cd raw
#gunzip *.gz;

echo "modifing fastq header name"

gunzip *.fastq.gz


echo  "reparing"

for f in *_L001_R1_001.fastq;
        do SAMPLE=`basename ${f%%_L001*}`;
                echo "$SAMPLE";
                echo "reorient sequences that have become disordered using bbmap: repair.sh script";
                repair.sh in="$SAMPLE"_L001_R1_001.fastq in2="$SAMPLE"_L001_R2_001.fastq out=repaired/"$SAMPLE"_L001_R1_001_repaired.fastq  out2=repaired/"$SAMPLE"_L001_R2_001_repaired.fastq outs="$SAMPLE"_L001_001_singles.fastq;
done;

cd repaired/

gzip *.fastq

for f in *_L001_R1_001_repaired.fastq.gz;
        do SAMPLE=`basename ${f%%_L001*}`;
                echo "$SAMPLE";
                echo "trimmommatic";
                trimmomatic PE "$SAMPLE"_L001_R1_001_repaired.fastq.gz "$SAMPLE"_L001_R2_001_repaired.fastq.gz "$SAMPLE"_forward_paired.fq.gz "$SAMPLE"_forward_unpaired.fq.gz "$SAMPLE"_reverse_paired.fq.gz "$SAMPLE"_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15;
                                fastqc "$SAMPLE"_forward_paired.fq.gz "$SAMPLE"_reverse_paired.fq.gz -o fastqc;

done;

mkdir after_trim/
mv *_reverse_paired.fq.gz after_trim/
mv *_forward_paired.fq.gz after_trim/

cd after_trim

       echo "spades"

#for f in *_forward_paired.fq.gz;
#        do SAMPLE=`basename ${f%%_forward*}`;
#                echo "$SAMPLE";
#               spades --rna -1 "$SAMPLE"_forward_paired.fq.gz -2 "$SAMPLE"_reverse_paired.fq.gz   -o "$SAMPLE"_spades  -t 32;
#done;

SqueezeMeta.pl -a spades  -m coassembly -p /home/drissi/Metat_complete/results_spades -s /home/drissi/Metat_complete/sample_file.txt -f /home/drissi/Metat_complete/raw/repaired/after_trim/ --nobin --nopfam -t 40 --euk

	echo "finished"
