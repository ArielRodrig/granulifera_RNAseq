# estimate % heterozygositty using publiched VCF file
vcftools --vcf /data/biolinux/Projects/Oophaga/BayeScan/Oophaga.filtered.bialelic.DP10.missing50.recode.vcf --het --out /data/biolinux/Projects/Oophaga_DFG/granulifera/eyes/assembly/evidentialgene/vcftools_oophaga_output.het
# vcftools will output the mumber of mozygotic positions and sites.
#Open the table in excel and calculate the hetherozygosity: =((N_SITES-O(HOM))/N_SITES)*100
# In this case we have 21-22% transcriptomic hetherozygosity across the 15 samples.

evigene="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/evidentialgene"
mkdir $evigene
cd $evigene
myspecies="granulifera"
assemblies="/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/assemblies"

cd $evigene
# SOAP-denovo
# remove extra info from headers (which contain semicolons in Oases and SOAP)
/data/biolinux/Software/evigene_new/evigene/scripts/rnaseq/trformat.pl -pre granulifera -format=soapt -out $assemblies/$myspecies"_SOAP_transcripts.fasta" -log -in $assemblies/"SOAP_k"*".transcrips.fasta"
sed 's/nt\=*.\; cf\=..*\; //g' $assemblies/$myspecies"_SOAP_transcripts.fasta" > $assemblies/$myspecies"_SOAP_transcripts.ok.fasta"
seqkit rmdup -s $assemblies/$myspecies"_SOAP_transcripts.ok.fasta" > $myspecies"_SOAP_transcripts.fa" 
[INFO] 768242 duplicated records removed

# OASES
# remove extra info from headers (which contain semicolons in Oases and SOAP)
/data/biolinux/Software/evigene_new/evigene/scripts/rnaseq/trformat.pl -pre granulifera -format=velvet -out $assemblies/$myspecies"_OASES_transcripts.fasta" -log -in $assemblies/"oases_k"*"transcrips.fasta"
sed 's/nt\=*.\; cf\=..*\; //g' $assemblies/$myspecies"_OASES_transcripts.fasta" > $assemblies/$myspecies"_OASES_transcripts.ok.fasta"
# Now remove duplicates and sequences shorter than 200
seqkit rmdup -s $assemblies/$myspecies"_OASES_transcripts.ok.fasta" > $myspecies"_OASES_transcripts.fa" 
[INFO] 211505 duplicated records removed

# SPADES
/data/biolinux/Software/evigene_new/evigene/scripts/rnaseq/trformat.pl -pre granuliferaspad -out $myspecies"_SPADES_transcripts.fa" -log -in $assemblies/"spades"*"transcripts.fasta"
# Trinity
/data/biolinux/Software/evigene_new/evigene/scripts/rnaseq/trformat.pl -pre granulifera -format=trinity -out $myspecies"_Trinity_transcripts.fa" -log -in $assemblies/"trinity."*".fasta"

# Idba_tran
/data/biolinux/Software/evigene_new/evigene/scripts/rnaseq/trformat.pl -pre granulifera -format=idba -MINTR=200 -out $myspecies"_Idba_transcripts.fa" -log -in $assemblies/"idba_"*"transcripts.fasta"

seqkit stats $myspecies"_SOAP_transcripts.fa" 
seqkit stats $myspecies"_OASES_transcripts.fa" 
seqkit stats $myspecies"_SPADES_transcripts.fa"
seqkit stats $myspecies"_Idba_transcripts.fa"
seqkit stats $myspecies"_Trinity_transcripts.fa"

############################# QC of input transcriptomes ##########################

## BUSCO assesment
cd /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/BUSCO
for fasta in /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/evidentialgene/granulifera_*_transcripts.fa; do
tx="${fasta##*/}"
/data/biolinux/Software/busco/bin/busco --config /data/biolinux/Software/busco/config/config.ini \
-m tran -i $fasta -o $tx".BUSCO" -l vertebrata_odb10 -c 30
done

# rnaQUAST assesment
export PATH="/data/biolinux/Software/GeneMark/gmst_linux_64:$PATH"
export PATH="/data/biolinux/Software/metaeuk/bin:$PATH"
python /data/biolinux/Software/rnaQUAST-2.2.2/rnaQUAST.py --transcripts $evigene/*.fa --output_dir /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/QUAST --threads 24 --labels Idba_trans OASES SOAP_denovo rnaSPADES Trinity



#########################################################################################
#											
#											
#			Evidential Gene							
#											
#########################################################################################
cat transcript_collection.fa | seqkit rmdup -n -d transcript_collection.dupIds -D transcript_collection.dupIds.details -o transcript_collection.dedupedIds.fa

ls $evigene/*.fa
cat $evigene/*.fa > $evigene/transcript_collection.fa
cd $evigene
/data/biolinux/Software/evigene_new/evigene/scripts/prot/tr2aacds.pl -tidy -NCPU 24 -MAXMEM 122880 -MINAA 70 -reorient -species Oophaga_granulifera -pHeterozygosity=2 -log -cdna transcript_collection.fa

# copy and rename the output 
cp /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/evidentialgene/okayset/transcript_collection.okay.aa /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/evidentialgene/okayset/transcript_collection.okay.aa.fasta
/data/biolinux/Software/bbmap/reformat.sh in=/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/evidentialgene/okayset/transcript_collection.okay.aa.fasta out=/data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/evidentialgene/okayset/transcript_collection.okay.pep.fasta trd=t


#########################################################################################
#											
#											
#			completeness and contiguity assesment							
#											
#########################################################################################

# BUSCO completeness assesment
cd /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/BUSCO
/data/biolinux/Software/busco/bin/busco --config /data/biolinux/Software/busco/config/config.ini -m prot -i /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/evidentialgene/okayset/transcript_collection.okay.pep.fasta -o EvidentialGene.okayset.BUSCO -l vertebrata_odb10 -c 30 -f

# rnaQUAST contiguity assesment
mkdir /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/QUAST
cd /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/QUAST
export PATH="/data/biolinux/Software/GeneMark/gmst_linux_64:$PATH"
export PATH="/data/biolinux/Software/metaeuk/bin:$PATH"
python /data/biolinux/Software/rnaQUAST-2.2.2/rnaQUAST.py --transcripts $evigene/okayset/transcript_collection.okay.mrna --output_dir /data/biolinux/Projects/Oophaga_DFG/granulifera/all_tissues/assembly/QUAST --threads 24 --labels EvidentialGene


# Now lets make a BUSCO plot comparing all the input assemblies and the Evigene consensus:
mkdir /data/biolinux/Projects/Oophaga_DFG/granulifera/eyes/assembly/BUSCO/busco_plots
# Copy all files BUSCO.txt files to the busco_plots directory and run:
python3 /data/biolinux/Software/busco/scripts/generate_plot.py -wd "/data/biolinux/Projects/Oophaga_DFG/granulifera/eyes/assembly/BUSCO/busco_plots"


