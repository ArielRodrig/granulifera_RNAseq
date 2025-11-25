# --------------- Preprocessing --------------#
#

mkdir /data/biolinux/Projects/Oophaga_DFG/granulifera/
mkdir /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues
mkdir /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly
mkdir /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/pre-processed_reads

RAWREADS="/data/g-zoologie/LINUX_VM/Projects/Oophaga_DFG/RNASeq_reads/granulifera"
OUTPUTDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/pre-processed_reads"

cd /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/pre-processed_reads
for f1 in $RAWREADS/*_R1_001.fastq.gz; do
    f2=${f1%%_R1_001.fastq.gz}"_R2_001.fastq.gz"
    o1="${f1##*/}"
    o2="${f2##*/}"
    fastp --thread 16 --detect_adapter_for_pe --in1 $f1 --in2 $f2 --out1 $OUTPUTDIR/${o1%%.fastq.gz}".trim.fastq.gz" --out2 $OUTPUTDIR/${o2%%.fastq.gz}".trim.fastq.gz" --html ${o1%%_R1_001.fastq.gz}".fastp.html" --json ${o1%%_R1_001.fastq.gz}".fastp.json"
done


# R-correct reads 
for f1 in $OUTPUTDIR/*_R1_001.trim.fastq.gz; do
    f2=${f1%%_R1_001.trim.fastq.gz}"_R2_001.trim.fastq.gz"
    perl /path/to/run_rcorrector.pl \
    -t 16 \
    -1 $f1 \
    -2 $f2 \
    -verbose \
    -od "$OUTPUTDIR"
done 

# Contamination screening with Kraken2
# krakendb was configured with the standard database plus 
# the mouse, human, plant and fungal sequences 

KRAKENDB="${KRAKENDB:-/path/to/kraken2_db}" 
for f1 in $OUTPUTDIR/*_R1_001.trim.fastq.cor.gz; do
    f2=${f1%%_R1_001.trim.fastq.cor.gz}"_R2_001.trim.fastq.cor.gz"
	o1="${f1##*/}"
    o2="${f2##*/}"
kraken2 \
  --db "$KRAKENDB" \
  --threads "$OMP_NUM_THREADS" \
  --paired  \
  --classified-out ${f1%%_R1_001.trim.fastq.cor.gz}".trim.cor.contam_R#.fq.gz" \
  --unclassified-out ${f1%%_R1_001.trim.fastq.cor.gz}".trim.cor.clean_R#.fq.gz" \
  --gzip-compressed \
  --confidence 0.5 \
  --output "${o1}.kraken.out" \
  --report "${o1}.kraken.report.txt" \
  --use-names $f1 $f2
done 


# Pool by locality and normalize 

localities="DAM ESQ PAL SAN"

for loc in $localities; do
cat $OUTPUTDIR/*$loc*".trim.cor.clean_R1.fq.gz" > $OUTPUTDIR/$loc".trim.cor.clean_R1.fq.gz"
cat $OUTPUTDIR/*$loc*".trim.cor.clean_R2.fq.gz" > $OUTPUTDIR/$loc".trim.cor.clean_R2.fq.gz"
    f1=$OUTPUTDIR/$loc".trim.cor.clean_R1.fq.gz"
    f2=$OUTPUTDIR/$loc".trim.cor.clean_R2.fq.gz"
bbnorm.sh -Xmx100g in1=$f1 in2=$f2 out1=${f1%%_R1.fq.gz}".norm_R1.fastq.gz" out2=${f2%%_R2.fq.gz}".norm_R2.fastq.gz" target=50 mindepth=2 prefilter=t;
rm $f1
rm $f2
done


#	remove ribosomal RNA with SortMeRNA

WORKDIR=$OUTPUTDIR"/sortmerna"
cd $OUTPUTDIR
for f1 in *.trim.cor.clean.norm_R1.fastq.gz; do
    f2=${f1%%.trim.cor.clean.norm_R1.fastq.gz}".trim.cor.clean.norm_R2.fastq.gz"
mkdir $WORKDIR
/data/biolinux/Software/sortmerna/bin/sortmerna --ref /data/biolinux/Software/sortmerna/database/smr_v4.3_fast_db.fasta \
--reads $f1 \
--reads $f2 \
--threads 32 -workdir $WORKDIR -blast 1 \
-v 1 -fastx 1 -paired_in 1 --zip-out --out2 1 \
--aligned ${f1%%.trim.cor.clean.norm_R1.fastq.gz}".trim.cor.norm.rRNAcontaminated" \
--other ${f1%%.trim.cor.clean.norm_R1.fastq.gz}".trim.cor.norm.rRNAfiltered"
rm -r $WORKDIR
done
