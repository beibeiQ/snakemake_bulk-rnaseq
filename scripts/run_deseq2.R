suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(DESeq2)
  library(tidyverse)
})

option_list <- list(
  make_option("--salmon_dir", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--tx2gene", type = "character"),
  make_option("--design", type = "character"),
  make_option("--condition_col", type = "character"),
  make_option("--contrasts", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--pca", type = "character"),
  make_option("--norm", type = "character"),
  make_option("--reference_level", type = "character", default = "")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

meta <- read_tsv(opt$metadata, show_col_types = FALSE)
stopifnot(all(c("sample", opt$condition_col) %in% colnames(meta)))

quant_files <- file.path(opt$salmon_dir, meta$sample, "quant.sf")
names(quant_files) <- meta$sample
if (any(!file.exists(quant_files))) {
  stop("以下 quant.sf 不存在: ", paste(quant_files[!file.exists(quant_files)], collapse = ", "))
}

tx2gene <- read_tsv(opt$tx2gene, show_col_types = FALSE)
stopifnot(all(c("TXNAME", "GENEID") %in% colnames(tx2gene)))

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene[, c("TXNAME", "GENEID")])

dds <- DESeqDataSetFromTximport(txi = txi, colData = meta, design = as.formula(opt$design))

if (nzchar(opt$reference_level)) {
  ref <- strsplit(opt$reference_level, ",")[[1]]
  if (length(ref) == 2 && ref[1] %in% colnames(colData(dds))) {
    dds[[ref[1]]] <- relevel(as.factor(dds[[ref[1]]]), ref = ref[2])
  }
}

dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE) |>
  as.data.frame() |>
  rownames_to_column("gene_id")
write_tsv(norm_counts, opt$norm)

vsd <- vst(dds, blind = FALSE)
pdf(opt$pca)
plotPCA(vsd, intgroup = opt$condition_col)
dev.off()

contrast_defs <- strsplit(opt$contrasts, ",")[[1]]
for (contrast in contrast_defs) {
  parts <- strsplit(contrast, ":")[[1]]
  if (length(parts) != 3) {
    warning("跳过无效 contrast: ", contrast)
    next
  }
  res <- results(dds, contrast = c(parts[1], parts[2], parts[3])) |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    arrange(padj)

  out_file <- file.path(opt$outdir, paste0(gsub(":", "_", contrast), "_DEG.tsv"))
  write_tsv(res, out_file)
}
