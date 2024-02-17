## PREPARE GENE ANNOTATION FOR REFERENCE GENOME  ------------------------------------

# download annotations and cds sequences from lepbase.org ---------------------------
wget http://download.lepbase.org/v4/features/Heliconius_melpomene_melpomene_Hmel2.5.gff3.gz ## for annotations (including gene annotations)
gunzip Heliconius_melpomene_melpomene_Hmel2.5.gff3.gz


# extract gene names and location from .gff3 file -----------------------------------
awk '{if ($3 == "gene") print }' Heliconius_melpomene_melpomene_Hmel2.5.gff3 | sed 's/=/\t/g' | awk '{print $10,$1,$4,$5,$7}' | sed 's/ /\t/g' > hmelv25_genes_coord.tab 


# create file with gene sequences ---------------------------------------------------
perl /n/home12/fseixas/code/heliconius_seixas/1.pseudo_references/2.3.gene_order/1.1.get_gene_seqs.pl Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa hmel2.5_genes_coord.tab 500 hmel2.5_genes.buf500.bp.fa

