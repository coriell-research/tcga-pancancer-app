#!/usr/bin/env bash
#
# Download PANCAN datasets from UCSC Xena
#
# -----------------------------------------------------------------------------
OUT=/home/gcalendo/data/tcga-pancan-app/raw_data

# Gene expression
wget -nc -O $OUT/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
wget -nc -O $OUT/tcga_RSEM_gene_fpkm.gz https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_fpkm.gz
wget -nc -O $OUT/tcga_RSEM_gene_tpm.gz https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz
wget -nc -O $OUT/tcga_RSEM_Hugo_norm_count.gz https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_Hugo_norm_count.gz

# Mutations
wget -nc -O $OUT/mc3.v0.2.8.PUBLIC.xena.gz https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/mc3.v0.2.8.PUBLIC.xena.gz

# Methylation
wget -nc -O $OUT/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz
wget -nc -O $OUT/probeMap-illuminaMethyl450_hg19_GPL16304_TCGAlegacy https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/probeMap%2FilluminaMethyl450_hg19_GPL16304_TCGAlegacy

# Copy number
wget -nc -O $OUT/TCGA.PANCAN.sampleMap-Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_data_by_genes.gz
wget -nc -O $OUT/TCGA.PANCAN.sampleMap-Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz

# Clinical data
wget -nc -O $OUT/Survival_SupplementalTable_S1_20171025_xena_sp https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp
wget -nc -O $OUT/TCGA_phenotype_denseDataOnlyDownload.tsv.gz https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz

# Immune types
wget -nc -O $OUT/Subtype_Immune_Model_Based.txt.gz https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Subtype_Immune_Model_Based.txt.gz
