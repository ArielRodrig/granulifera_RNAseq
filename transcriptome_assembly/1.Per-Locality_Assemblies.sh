
########################
#       ASSEMBLIES     #
########################

localities="DAM ESQ PAL SAN"
for loc in $localities; do
READSDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/pre-processed_reads"
WORKDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/assemblies"
mkdir $WORKDIR
READ1=$READSDIR/$loc".trim.cor.norm.rRNAfiltered_fwd.fq.gz"
READ2=$READSDIR/$loc".trim.cor.norm.rRNAfiltered_rev.fq.gz"
TRINITYDIR=$WORKDIR"/trinity_temp"
mkdir $TRINITYDIR
cd $TRINITYDIR
/data/biolinux/Software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left $READ1  --right $READ2 --min_kmer_cov 2 --SS_lib_type RF --CPU 30 --max_memory 100G --bflyCalculateCPU --full_cleanup --no_normalize_reads --output "trinity.$loc"
cp "trinity.$loc.Trinity.fasta" $WORKDIR
rm -r $TRINITYDIR
done




# rnaSPADES 
localities="DAM ESQ PAL SAN"
for loc in $localities; do
READSDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/pre-processed_reads"
WORKDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/assemblies"
GZIP_READ1=$READSDIR/$loc".trim.cor.norm.rRNAfiltered_fwd.fq.gz"
GZIP_READ2=$READSDIR/$loc".trim.cor.norm.rRNAfiltered_rev.fq.gz"
gunzip -k $GZIP_READ1
gunzip -k $GZIP_READ2
READ1=${GZIP_READ1%%.fq.gz}".fq"
READ2=${GZIP_READ2%%.fq.gz}".fq"
SPADESDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/assemblies/spades_temp"
mkdir $SPADESDIR
cd $WORKDIR
rnaspades.py -k auto -o $SPADESDIR --ss rf -1 $READ1 -2 $READ2 -m 124 -t 24 --disable-gzip-output  --checkpoints "last"
cp $SPADESDIR/"transcripts.fasta" $WORKDIR"/spades_"$loc".transcripts.fasta"
cp $SPADESDIR/"spades.log" $WORKDIR"/spades_"$loc".log"
rm -r $SPADESDIR
rm $READ1
rm $READ2
done

# OASES

localities="ESQ PAL SAN DAM"
for loc in $localities; do
READSDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/pre-processed_reads"
WORKDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/assemblies"
GZIP_READ1=$READSDIR/$loc".trim.cor.norm.rRNAfiltered_fwd.fq.gz"
GZIP_READ2=$READSDIR/$loc".trim.cor.norm.rRNAfiltered_rev.fq.gz"
gunzip -k $GZIP_READ1
gunzip -k $GZIP_READ2
READ1=${GZIP_READ1%%.fq.gz}".fq"
READ2=${GZIP_READ2%%.fq.gz}".fq"
TEMPDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/assemblies/oases_temp"
mkdir $TEMPDIR
kmerValues="19 37 55 73"
for kmer in $kmerValues; do
OASESDIR="$TEMPDIR/k$kmer"
mkdir $OASESDIR
cd $OASESDIR
echo VELVETH_START; velveth $OASESDIR $kmer -shortPaired -fastq -separate -strand_specific $READ1 $READ2; echo VELVETG_START; velvetg $OASESDIR -read_trkg yes -min_contig_lgth 150 -cov_cutoff 4 -ins_length 200 -clean yes; echo OASES_START; oases $OASESDIR -cov_cutoff 4 -ins_length 200 -min_trans_lgth 200
cp $OASESDIR/"transcripts.fa" $WORKDIR"/oases_k"$kmer"_"$loc".transcrips.fasta"
cp $OASESDIR/"Log" $WORKDIR"/oases_k"$kmer"_"$loc".Log"
cp $OASESDIR/"stats.txt" $WORKDIR"/oases_k"$kmer"_"$loc".stats.txt"
rm -r $OASESDIR
done
rm $READ1
rm $READ2
done


# SOAPdenovo-trans

localities="DAM ESQ PAL SAN"
for loc in $localities; do
READSDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/pre-processed_reads"
WORKDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/assemblies"
GZIP_READ1=$READSDIR/$loc".trim.cor.norm.rRNAfiltered_fwd.fq.gz"
GZIP_READ2=$READSDIR/$loc".trim.cor.norm.rRNAfiltered_rev.fq.gz"
gunzip -k $GZIP_READ1
gunzip -k $GZIP_READ2
READ1=${GZIP_READ1%%.fq.gz}".fq"
READ2=${GZIP_READ2%%.fq.gz}".fq"
TEMPDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/assemblies/soap_temp"
mkdir $TEMPDIR
kmerValues="19 37 55 73"
for kmer in $kmerValues; do
SOAPDIR="$TEMPDIR/k$kmer"
mkdir $SOAPDIR
cd $SOAPDIR
echo -e "#maximal read length\nmax_rd_len=150\n[LIB]\n#maximal read length in this lib\nrd_len_cutof=150\n#average insert size\navg_ins=200\n#if sequence needs to be reversed \nreverse_seq=0\n#in which part(s) the reads are used\nasm_flags=3\n#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)\nmap_len=32\n#fastq file for read 1\nq1=$READ1\n#fastq file for read 2\nq2=$READ2\n" > SOAP_Conf.txt

SOAPdenovo-Trans-127mer all -s SOAP_Conf.txt -K $kmer -o SOAP_k$kmer.$loc -L 200 -F -p 30
cp *.scafSeq $WORKDIR"/SOAP_k"$kmer"."$loc".transcrips.fasta"
cp *.scafStatistics $WORKDIR"/SOAP_k"$kmer"."$loc".scafStatistics"
rm -r $SOAPDIR
done
rm $READ1
rm $READ2
done

# Idba_Trans 
#  --no_correct option was used as the pre-processing steps included r-corrector per sample

localities="DAM ESQ PAL SAN"
for loc in $localities; do
READSDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/pre-processed_reads"
WORKDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/assemblies"
GZIP_READ1=$READSDIR/$loc".trim.cor.norm.rRNAfiltered_fwd.fq.gz"
GZIP_READ2=$READSDIR/$loc".trim.cor.norm.rRNAfiltered_rev.fq.gz"
gunzip -k -f $GZIP_READ1
gunzip -k -f $GZIP_READ2
READ1=${GZIP_READ1%%.fq.gz}".fq"
READ2=${GZIP_READ2%%.fq.gz}".fq"
TEMPDIR="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/assemblies/idba_temp"
mkdir $TEMPDIR
MERGEDREADS=${READ1%%.rRNAfiltered_fwd.fq}".merged.fq"
/data/biolinux/Software/ConSemble3+/idba-1.1.3/bin/fq2fa --merge --filter $READ1 $READ2 $MERGEDREADS
/data/biolinux/Software/ConSemble3+/idba-1.1.3/bin/idba_tran -o $TEMPDIR -r $MERGEDREADS --num_threads 18 --mink 19 --maxk 73 --step 18 --max_isoforms 50 --no_correct
cp $TEMPDIR/transcript-73.fa $WORKDIR/"idba_"$loc"_transcripts.fasta"
cp $TEMPDIR/log $WORKDIR/"idba_"$loc"log.txt"
rm -r $TEMPDIR
rm $READ1
rm $READ2
rm $MERGEDREADS
done



