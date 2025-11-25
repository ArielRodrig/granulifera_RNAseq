# run EnTAP annotation on the reference transcriptome

# After updating diamond, download the eggnog databases and configure each search databases 
diamond makedb --in /data/biolinux/GenBank/refseq/RefSeq_Amphibia.aa.fasta -d /data/biolinux/GenBank/refseq/refseq_amphibia
diamond makedb --in /data/biolinux/GenBank/uniprot/uniprot_sprot.fasta -d /data/biolinux/GenBank/uniprot/uniprot_sprot.dmnd
diamond makedb --in /data/biolinux/GenBank/uniprot/uniref90.fasta -d /data/biolinux/GenBank/uniprot/uniref90.dmnd
diamond makedb --in /data/biolinux/GenBank/uniprot/TrEMBL_Amphibia.fasta -d /data/biolinux/GenBank/uniprot/TrEMBL_Amphibia.dmnd
diamond dbinfo -d /data/biolinux/Software/EnTAP/databases/eggnog_proteins.dmnd
diamond dbinfo -d /data/biolinux/Software/EnTAP/bin/eggnog_proteins.dmnd

# Edit the entap_config.ini with deired settings and configure EnTAP
cd /data/biolinux/Software/EnTAP
EnTAP --config -t 8 --ini /data/biolinux/Software/EnTAP/entap_config.ini

WORKDIR=$HOME"/evidentialgene/okayset"
cd /data/biolinux/Software/EnTAP
# run EnTAP with settings speciffied in the entap_config.ini file

EnTAP --runP --ini /data/biolinux/Software/EnTAP/entap_config.ini -i \
$WORKDIR/transcript_collection.okay.mrna --out-dir $HOME/EnTAP_annotation/entap_wo_expression_filter \
-d /data/biolinux/GenBank/uniprot/uniprot_sprot.dmnd \
-d /data/biolinux/GenBank/uniprot/TrEMBL_Amphibia.dmnd \
-d /data/biolinux/GenBank/refseq/refseq_vert.dmnd \
-d /data/biolinux/GenBank/uniprot/uniref90.dmnd -t 24

# Extract the transcripts with hits on reference databases from the original EvidentialGene file:
# 1 Select transcripts with hits:
awk -F: '$2 != ""' $HOME/EnTAP_annotation/entap_outfiles/final_results/final_annotations_no_contam_lvl0.tsv > $HOME/EnTAP_annotation/entap_outfiles/final_results/final_annotations_no_contam_lvl0.hits.tsv
# USe the above file to create a list of transcripts and extract the sequences:
seqtk subseq $WORKDIR/transcript_collection.okay.mrna $HOME/EnTAP_annotation/entap_outfiles/final_results/final_annotations_no_contam_lvl0.hits.ids.txt > $WORKDIR/transcript_collection.okay.mrna.ORFs.annotations-no-contams.fasta
# remove extra spaces
cut -d" " -f1 $WORKDIR/transcript_collection.okay.mrna.ORFs.annotations-no-contams.fasta > $WORKDIR/transcript_collection.okay.mrna.ORFs.annotations-no-contams.shortheaders.fasta
# Index the coding 
samtools faidx $WORKDIR/transcript_collection.okay.mrna.ORFs.annotations-no-contams.shortheaders.fasta



