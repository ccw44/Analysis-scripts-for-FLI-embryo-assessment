# TPM and FPKM from the featureCounts output of count.sh.
#
#   Rscript count2tpm.R        (run inside the count output directory)
#
# Reads  gene_counts_raw.txt / transcript_counts_raw.txt
# Writes Calculated_Gene_TPM.csv, Calculated_Gene_FPKM.csv,
#        Calculated_Transcript_TPM.csv, Calculated_Transcript_FPKM.csv
#
# Scope: TPM/FPKM are for expression display only. Differential expression uses
# the raw counts from Gene_Counts_Clean.csv, never these values.

# FeatureCounts  TPM and FPKM

counts_data <- read.table("gene_counts_raw.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)


gene_lengths <- counts_data$Length       
gene_symbols <- counts_data$gene_name    

bam_cols <- grep("\\.bam$", colnames(counts_data))
count_matrix <- counts_data[, bam_cols]

colnames(count_matrix) <- gsub(".STAR.Aligned.sortedByCoord.out.bam", "", colnames(count_matrix))
colnames(count_matrix) <- gsub(".Hisat_aln.sorted.bam", "", colnames(count_matrix))


rpk <- count_matrix / (gene_lengths / 1000)

scaling_factors <- colSums(rpk) / 1e6

tpm_matrix <- t(t(rpk) / scaling_factors)

tpm_out <- data.frame(Gene_ID = rownames(tpm_matrix), Gene_Name = gene_symbols, tpm_matrix)
write.csv(tpm_out, "Calculated_Gene_TPM.csv", row.names = FALSE)


total_reads_millions <- colSums(count_matrix) / 1e6

fpkm_matrix <- t(t(rpk) / total_reads_millions)

fpkm_out <- data.frame(Gene_ID = rownames(fpkm_matrix), Gene_Name = gene_symbols, fpkm_matrix)
write.csv(fpkm_out, "Calculated_Gene_FPKM.csv", row.names = FALSE)



counts_data <- read.table("transcript_counts_raw.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

tx_lengths <- counts_data$Length       
gene_ids <- counts_data$gene_id
gene_symbols <- counts_data$gene_name    

bam_cols <- grep("\\.bam$", colnames(counts_data))
count_matrix <- counts_data[, bam_cols]

colnames(count_matrix) <- gsub(".STAR.Aligned.sortedByCoord.out.bam", "", colnames(count_matrix))
colnames(count_matrix) <- gsub(".Hisat_aln.sorted.bam", "", colnames(count_matrix))

rpk <- count_matrix / (tx_lengths / 1000)
scaling_factors <- colSums(rpk) / 1e6
tpm_matrix <- t(t(rpk) / scaling_factors)

tpm_out <- data.frame(Transcript_ID = rownames(tpm_matrix), Gene_ID = gene_ids, Gene_Name = gene_symbols, tpm_matrix)
write.csv(tpm_out, "Calculated_Transcript_TPM.csv", row.names = FALSE)


total_reads_millions <- colSums(count_matrix) / 1e6
fpkm_matrix <- t(t(rpk) / total_reads_millions)

fpkm_out <- data.frame(Transcript_ID = rownames(fpkm_matrix), Gene_ID = gene_ids, Gene_Name = gene_symbols, fpkm_matrix)
write.csv(fpkm_out, "Calculated_Transcript_FPKM.csv", row.names = FALSE)
