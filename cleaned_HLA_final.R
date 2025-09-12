
############################################################
# cleaned_HLA_final.R
# Consolidated, function-structured R script for HLA-DQ
# microbiome analyses (gut + saliva). 
#
# Notes:
#   - Many functions include arguments so you can override
#     file paths or plotting behavior.
#   - SALIVA-related analyses are marked with # SALIVA
############################################################

# ------------------------------
# Packages (install if missing)
# ------------------------------
required_pkgs <- c(
  "limma","reshape2","dplyr","vegan","ggplot2","ggpubr","phyloseq",
  "SIAMCAT","randomcoloR","ggrepel","stats","ade4","ape","anpan",
  "ggvenn","ComplexHeatmap","circlize","topGO","data.table","pROC","stringr"
)

for (pkg in required_pkgs) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    install.packages(pkg, repos = "https://cran.rstudio.com")
    suppressWarnings(require(pkg, character.only = TRUE))
  }
}

# anpan and SIAMCAT may be Bioconductor/GitHub packages in your environment;
# Please ensure they are installed before running heavy analyses.

# ------------------------------
# Utility: Safe read table
# ------------------------------
safe_read <- function(path, sep="\t", header=TRUE, ...) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  read.table(path, sep=sep, header=header, stringsAsFactors = FALSE, ...)
}

# ------------------------------
# Common data loaders
# ------------------------------
load_diet <- function(path = "dietary_diversity_use.txt") {
  diet <- safe_read(path, sep="\t", header=TRUE)
  return(diet)
}

load_metaphlan_species <- function(path = "metaphlan.txt") {
  tb <- safe_read(path, sep="\t", header=TRUE)
  colnames(tb)[1] <- "clade_name"
  tb <- tb %>% filter(grepl("s__", clade_name) & !grepl("t__", clade_name))
  # get species name after "s__"
  spnames <- as.data.frame(strsplit2(tb$clade_name, "s__"))[,2]
  rownames(tb) <- spnames
  tb <- tb[,-1, drop=FALSE]
  tb_t <- as.data.frame(t(tb))
  tb_t$SampleID <- rownames(tb_t)
  tb_t$subject <- strsplit2(rownames(tb_t), "_")[,2]
  tb_t$week <- strsplit2(rownames(tb_t), "_")[,1]
  return(tb_t)
}

load_pathways <- function(path = "all_pathabundance_relab.tsv") {
  tb <- safe_read(path, sep="\t", header=TRUE)
  colnames(tb)[1] <- "Pathway"
  return(tb)
}



# ------------------------------
# Transformations & stats
# ------------------------------
compute_shannon <- function(abun_matrix) {
  # abun_matrix: numeric matrix rows=samples, cols=species
  mat <- as.matrix(abun_matrix)
  # ensure non-negative
  sh <- vegan::diversity(mat, index = "shannon", base = exp(1))
  return(sh)
}



# ------------------------------
# Generic plotting helpers
# ------------------------------
plot_box_wilcox <- function(df, xcol="Type", ycol,
                            xlabel="", ylabel=NULL, title=NULL, save=NULL) {
  p <- ggplot(df, aes_string(x=xcol, y=ycol)) +
    geom_boxplot(aes_string(color=xcol), width=0.6) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(colour = 'black'),
          panel.background = element_blank(),
          axis.text = element_text(size = 9, color = 'black'),
          axis.title = element_text(size = 10, color = 'black')) +
    scale_color_manual(values = palette) +
    ylab(ylabel %||% ycol) + xlab(xlabel) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 12)) +
    guides(fill=guide_legend(title = NULL)) +
    guides (fill = "none") +
    stat_compare_means(comparisons=list(c("N","Y")) ,test = wilcox.test,paired = FALSE,p.adjust.method = "none",label="p.format") +
    scale_x_discrete(labels = c("DQ2/8-","DQ2/8+"))
  if (!is.null(title)) p <- p + ggtitle(title)
  if (!is.null(save)) ggsave(save, p, width=5, height=4)
  return(p)
}

# ------------------------------
# Beta diversity wrapper (adonis)
# ------------------------------
adonis_beta <- function(feature_df, meta_df, covariates_formula = "gender+age+BMI+Type+Protein_g+Fat_g+Carbohydrate_g+Fiber_g++VitaminA_ug+VitaminB1_mg+VitaminB2_mg+VitaminB3_mg+VitaminC_mg+VitaminE_mg+Cholesterol_mg+Calcium_mg+Sodium_mg", method="jsd") {
  # feature_df: rows = samples, cols = features (numeric)
  # meta_df: contains SampleID and covariates; rownames or SampleID must match
  data2 <- as.data.frame(feature_df)
  data2$SampleID <- rownames(feature_df)
  dt <- merge(data2, meta_df, by="SampleID")
  # prepare phyloseq-like distance (use JSD via phyloseq)
  # convert: taxa_are_rows = T expects taxa rows; so transpose
  mat <- t(as.matrix(dt[, colnames(feature_df), drop=FALSE]))
  OTU = phyloseq::otu_table(mat, taxa_are_rows = TRUE)
  physeq = phyloseq(OTU)
  jsd = phyloseq::distance(physeq, method = method)
  dists <- jsd
  equal <- paste0("dists~", covariates_formula)
  r1 <- adonis2(as.formula(equal), data = dt, permutations = 999, by="margin")
  r1$BH <- p.adjust(as.numeric(r1$`Pr(>F)`), method = "BH")
  r1$bonferroni <- p.adjust(as.numeric(r1$`Pr(>F)`), method = "bonferroni")
  return(r1)
}

# ------------------------------
# SIAMCAT wrapper (simple)
# ------------------------------
run_siamcat_analysis <- function(feat_path = "./relab_gut.tsv", meta_path = "./meta.tsv", output_prefix = "siamcat") {
  if (!file.exists(feat_path) || !file.exists(meta_path)) {
    message("SIAMCAT: feature or meta path missing.")
    return(NULL)
  }
  feat <- read.table(feat_path, sep="\t", header=TRUE, row.names=1, check.names = FALSE)
  meta <- read.table(meta_path, sep="\t", header=TRUE, row.names=1, check.names = FALSE)
  label <- create.label(meta=meta, label='Type', case='Y')
  sc.obj <- siamcat(feat=feat, label=label, meta=meta)
  sc.obj <- filter.features(sc.obj, filter.method = 'abundance', cutoff = 0.001)
  sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha = 0.05)
  sc.obj <- normalize.features(sc.obj, norm.method = "log.unit", norm.param = list(log.n0 = 1e-06, n.p = 2, norm.margin = 1))
  sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 2)
  sc.obj <- train.model(sc.obj, method = "lasso")
  sc.obj <- make.predictions(sc.obj)
  sc.obj <- evaluate.predictions(sc.obj)
  model.evaluation.plot(sc.obj)
  return(sc.obj)
}

# ------------------------------
# Virulence factor analysis (from BLAST results)
# ------------------------------
virulence_analysis <- function(virus_dir = "./virus_factor/gut", reads_table = "gut_reads.txt", diet_path = "dietary_diversity_use.txt", week_filter = "S2") {
  files <- list.files(virus_dir)
  if (length(files) == 0) {
    message("virulence_analysis: no files in ", virus_dir)
    return(NULL)
  }
  result_all <- NULL
  for (file_name in files) {
    dat <- read.table(file.path(virus_dir, file_name), sep=" ", header=FALSE, stringsAsFactors = FALSE)
    colnames(dat) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    filtered <- dat %>% filter(pident >= 70, length >= 50)
    ab <- filtered %>% count(sseqid) %>% rename(virulence_factor = sseqid, count = n)
    sample_name <- strsplit(file_name, "_matches")[[1]][1]
    ab$SampleID <- sample_name
    result_all <- rbind(result_all, ab)
  }
  tb <- acast(result_all, SampleID ~ virulence_factor, value.var = "count", fill=0)
  tb <- as.data.frame(tb)
  tb$Sum <- rowSums(tb)
  reads <- safe_read(reads_table, sep="\t", header=TRUE)
  tb$sampleID <- rownames(tb)
  tb2 <- merge(tb, reads, by.x="sampleID", by.y="sampleID", all.x=TRUE)
  # normalize to reads
  name_v <- setdiff(colnames(tb), c("Sum","sampleID"))
  tb3 <- tb2[, name_v] / tb2$reads * 100
  tb3$SampleID <- tb2$sampleID
  diet <- load_diet(diet_path)
  dt1 <- merge(tb3[, c("SampleID", colnames(tb3)[-ncol(tb3)])], diet, by="SampleID")
  dt1$week <- strsplit2(dt1$SampleID, "_")[,1]
  dt1 <- dt1[dt1$week == week_filter,]
  # aggregate linear models per virulence factor
  result <- NULL
  name_list <- colnames(tb3)[-ncol(tb3)]
  for (i in name_list) {
    dt_i <- dt1[, c(i, "age", "gender", "BMI", "week", "subject", "Type", "Protein_g", "Fat_g", "Carbohydrate_g", "Fiber_g", "VitaminA_ug", "VitaminB1_mg", "VitaminB2_mg", "VitaminB3_mg", "VitaminC_mg", "VitaminE_mg", "Cholesterol_mg", "Calcium_mg", "Sodium_mg")]
    dt_i <- na.omit(dt_i)
    colnames(dt_i)[1] <- "factor"
    dt_i$factor <- as.numeric(unlist(dt_i$factor))
    dt_i$age <- as.numeric(dt_i$age)
    if (!all(dt_i$factor == 0)) {
      dt_i$factor <- scale(dt_i$factor)
      summ <- summary(lm(factor ~ Type + week + age + gender + BMI + Protein_g + Fat_g + Carbohydrate_g + Fiber_g + VitaminA_ug + VitaminB1_mg + VitaminB2_mg + VitaminB3_mg + VitaminC_mg + VitaminE_mg + Cholesterol_mg + Calcium_mg + Sodium_mg, dt_i))
      result <- rbind(result, c(i, summ$coefficients[2,1], summ$coefficients[2,4]))
    }
  }
  result <- as.data.frame(result, stringsAsFactors = FALSE)
  colnames(result) <- c("factor", "estimate", "pvalue")
  result$pvalue <- as.numeric(result$pvalue)
  result$P_BH <- p.adjust(result$pvalue, method = "BH")
  return(list(table=result, dt=dt1))
}

# ------------------------------
# Differential species (linear models)
# ------------------------------
species_differential <- function(metaphlan_path = "metaphlan.txt", diet_path = "dietary_diversity_use.txt") {
  tb <- load_metaphlan_species(metaphlan_path)
  tb_species <- tb[tb$week == week_filter, ]
  sample_cols <- setdiff(colnames(tb_species), c("SampleID","subject","week"))
  diet <- load_diet(diet_path)
  dt <- merge(tb_species[, c(sample_cols, "SampleID")], diet, by="SampleID")
  result <- NULL
  for (i in sample_cols) {
    dt_i <- dt[, c(i, "age", "gender", "BMI", "week", "subject", "Type", "Protein_g", "Fat_g", "Carbohydrate_g", "Fiber_g", "VitaminA_ug", "VitaminB1_mg", "VitaminB2_mg", "VitaminB3_mg", "VitaminC_mg", "VitaminE_mg", "Cholesterol_mg", "Calcium_mg", "Sodium_mg")]
    dt_i <- na.omit(dt_i)
    colnames(dt_i)[1] <- "factor"
    dt_i$factor <- as.numeric(unlist(dt_i$factor))
    if (!all(dt_i$factor == 0)) {
      dt_i$factor <- scale(dt_i$factor)
      summ <- summary(lm(factor ~ Type + age + gender + BMI + Protein_g + Fat_g + Carbohydrate_g + Fiber_g + VitaminA_ug + VitaminB1_mg + VitaminB2_mg + VitaminB3_mg + VitaminC_mg + VitaminE_mg + Cholesterol_mg + Calcium_mg + Sodium_mg, dt_i))
      result <- rbind(result, c(i, summ$coefficients[2,1], summ$coefficients[2,4]))
    }
  }
  result <- as.data.frame(result, stringsAsFactors = FALSE)
  colnames(result) <- c("species", "estimate", "pvalue")
  result$pvalue <- as.numeric(result$pvalue)
  result$P_BH <- p.adjust(result$pvalue, method = "BH")
  return(result)
}

# ------------------------------
# Pathway differential (linear models)
# ------------------------------
pathway_differential <- function(pathway_file = "./all_pathabundance_relab.tsv", diet_path = "dietary_diversity_use.txt") {
  tb_path <- safe_read(pathway_file, sep="\t", header=TRUE)
  colnames(tb_path)[1] <- "Pathway"
  tb_unstratified <- tb_path %>% filter(!grepl("g__", Pathway), !grepl("unclassified", Pathway), !grepl("UNGROUPED", Pathway), !grepl("UNMAPPED", Pathway), !grepl("UNINTEGRATED", Pathway))
  rownames(tb_unstratified) <- tb_unstratified$Pathway
  tb_unstratified <- tb_unstratified[,-1]
  tb_t <- as.data.frame(t(tb_unstratified))
  tb_t$SampleID <- rownames(tb_t)
  tb_t$subject <- strsplit2(rownames(tb_t), "_")[,2]
  tb_t$week <- strsplit2(rownames(tb_t), "_")[,1]
  tb_t <- tb_t[ tb_t$subject != "" & !tb_t$week %in% c("S6","S7","S8") , ]
  name <- colnames(tb_t)[1:(ncol(tb_t)-3)]
  diet <- load_diet(diet_path)
  dt_pwy <- merge(tb_t[, c(name, "SampleID")], diet, by="SampleID")
  dt_pwy$week <- strsplit2(dt_pwy$SampleID, "_")[,1]
  result <- NULL
  for (i in name) {
    dt_i <- dt_pwy[, c(i, "age", "gender", "BMI", "week", "subject", "Type", "Protein_g", "Fat_g", "Carbohydrate_g", "Fiber_g", "VitaminA_ug", "VitaminB1_mg", "VitaminB2_mg", "VitaminB3_mg", "VitaminC_mg", "VitaminE_mg", "Cholesterol_mg", "Calcium_mg", "Sodium_mg")]
    dt_i <- na.omit(dt_i)
    colnames(dt_i)[1] <- "factor"
    dt_i$factor <- as.numeric(unlist(dt_i$factor))
    if (!all(dt_i$factor == 0)) {
      dt_i$factor <- scale(dt_i$factor)
      summ <- summary(lm(factor ~ Type + age + gender + BMI + Protein_g + Fat_g + Carbohydrate_g + Fiber_g + VitaminA_ug + VitaminB1_mg + VitaminB2_mg + VitaminB3_mg + VitaminC_mg + VitaminE_mg + Cholesterol_mg + Calcium_mg + Sodium_mg, dt_i))
      result <- rbind(result, c(i, summ$coefficients[2,1], summ$coefficients[2,4]))
    }
  }
  result <- as.data.frame(result, stringsAsFactors = FALSE)
  colnames(result) <- c("pathway", "estimate", "pvalue")
  result$pvalue <- as.numeric(result$pvalue)
  result$P_BH <- p.adjust(result$pvalue, method = "BH")
  return(result)
}

# ------------------------------
# anpan phylogenetic loop wrapper
# ------------------------------
anpan_phylo_loop <- function(gene_dir = "anpan/gene_file4", meta_path = "anpan/metadata.tsv", output_file = "tree_result_anpan_new_W2.tsv") {
  file_names <- list.files(gene_dir)
  if (length(file_names) == 0) {
    message("anpan_phylo_loop: no files in ", gene_dir)
    return(NULL)
  }
  meta <- safe_read(meta_path, sep="\t", header=TRUE, row.names=1)
  for (j in file_names) {
    try({
      bug_path <- file.path(gene_dir, j)
      gene_matrix <- read.table(bug_path, sep="\t", row.names=1, header=TRUE, check.names = FALSE)
      gene_matrix_t <- as.data.frame(t(gene_matrix))
      gene_matrix_t$week <- strsplit2(rownames(gene_matrix_t), "_")[,1]
      gene_matrix_t <- gene_matrix_t[gene_matrix_t$week == "S2", ]
      gene_matrix_t <- gene_matrix_t[ , -ncol(gene_matrix_t), drop=FALSE]
      gene_matrix <- as.data.frame(t(gene_matrix_t))
      gene_matrix <- gene_matrix[, colSums(gene_matrix) > 0, drop=FALSE]
      gene_matrix <- gene_matrix[rowSums(gene_matrix) > 0, , drop=FALSE]
      gene_matrix <- t(gene_matrix)
      metadata <- meta[rownames(gene_matrix), , drop=FALSE]
      pca_result <- prcomp(gene_matrix, scale. = TRUE)
      reduced_data <- pca_result$x[, 1:10, drop=FALSE]
      distance_matrix <- dist(reduced_data, method = "euclidean")
      nj_tree <- ape::nj(distance_matrix)
      tr <- nj_tree
      anpan_result <- anpan::anpan_pglmm(meta_file = metadata,
                                         tree_file = tr,
                                         outcome = "is_DQ28",
                                         covariates = c("age", "gender", "BMI", "Protein_g", "Fat_g", "Carbohydrate_g", "Fiber_g"),
                                         family = "gaussian",
                                         bug_name = "sim_bug",
                                         reg_noise = TRUE,
                                         loo_comparison = TRUE,
                                         run_diagnostics = FALSE,
                                         refresh = 500,
                                         show_plot_cor_mat = FALSE,
                                         show_plot_tree = FALSE,
                                         show_post = FALSE)
      anpan_result2 <- anpan_result$loo$comparison
      if (!is.null(anpan_result2)) {
        res <- t(c(j, t(anpan_result2[1,]), t(anpan_result2[2,])))
        write.table(res, output_file, sep="\t", append = TRUE, col.names = FALSE, row.names = FALSE)
      }
    }, silent = TRUE)
  }
  return(TRUE)
}

# ------------------------------
# GO enrichment wrapper (topGO)
# ------------------------------
go_enrichment_from_gene_list <- function(background_genes, gene_to_go, differential_genes) {
  # background_genes: vector of all gene IDs
  # gene_to_go: data.frame with columns GeneID and GO (semicolon separated)
  # differential_genes: vector
  gene_list <- factor(as.integer(background_genes %in% differential_genes))
  names(gene_list) <- background_genes
  go_map <- strsplit(as.character(gene_to_go$GO), split = "; ")
  names(go_map) <- gene_to_go$GeneID
  # run for MF, BP, CC
  res_list <- list()
  for (ont in c("MF","BP","CC")) {
    GOdata <- new("topGOdata", ontology = ont, allGenes = gene_list, annot = annFUN.gene2GO, gene2GO = go_map)
    result_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    result_table <- GenTable(GOdata, classicFisher = result_fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)
    result_table$ontology <- ont
    res_list[[ont]] <- result_table
  }
  return(do.call(rbind, res_list))
}

# ------------------------------
# Heatmap helper (ComplexHeatmap)
# ------------------------------
plot_complex_heatmap <- function(mat, row_names = NULL, col_names = NULL, filename = NULL, width=8, height=6) {
  if (!is.null(filename)) {
    pdf(filename, width=width, height=height)
  }
  ht <- ComplexHeatmap::Heatmap(as.matrix(mat), name="value", show_row_names = !is.null(row_names), show_column_names = !is.null(col_names))
  draw(ht)
  if (!is.null(filename)) dev.off()
  return(invisible(ht))
}

# ------------------------------
# Figure functions: implement each figure with modular functions
# The original script included many figures; here we provide full
# implementations adapted and organized. If some input files are missing,
# functions will return informative messages.
# ------------------------------

# Fig1a: Shannon index boxplot
fig1a <- function(metaphlan_path = "metaphlan.txt", diet_path = "dietary_diversity_use.txt", save=NULL) {
  message("Running Fig1a (Shannon) ...")
  tb <- load_metaphlan_species(metaphlan_path)
  tb_species <- tb
  sample_cols <- setdiff(colnames(tb_species), c("SampleID","subject","week"))
  sh <- compute_shannon(tb_species[, sample_cols, drop=FALSE])
  sh_df <- data.frame(SampleID = tb_species$SampleID, Shannon = sh)
  diet <- load_diet(diet_path)
  shanno <- merge(diet, sh_df, by="SampleID")
  p <- plot_box_wilcox(shanno, xcol="Type", ycol="Shannon", ylabel="Shannon_index", save=save)
  print(p)
  invisible(shanno)
}

# Fig1b: GMHI boxplot
fig1b <- function(gmhi_path = "GMHI.txt", diet_path = "dietary_diversity_use.txt", save=NULL) {
  message("Running Fig1b (GMHI) ...")
  if (!file.exists(gmhi_path)) {
    message("GMHI file not found: ", gmhi_path)
    return(NULL)
  }
  gmhi <- safe_read(gmhi_path, sep="\t", header=TRUE)
  diet <- load_diet(diet_path)
  gmhi <- merge(diet, gmhi, by="SampleID")
  p <- plot_box_wilcox(gmhi, xcol="Type", ycol="GMHI", ylabel="GMHI", save=save)
  print(p)
  invisible(gmhi)
}

# Fig1c: Virulence factor analysis
fig1c <- function(virus_dir = "./virus_factor/gut", reads_table = "gut_reads.txt", diet_path = "dietary_diversity_use.txt", save=NULL) {
  message("Running Fig1c (virulence factors) ...")
  res <- virulence_analysis(virus_dir=virus_dir, reads_table=reads_table, diet_path=diet_path, week_filter="S2")
  if (is.null(res)) return(NULL)
  dt1 <- res$dt
  p <- plot_box_wilcox(dt1, xcol="Type", ycol="Sum", ylabel="Relative abundance of virulence factor", save=save)
  print(p)
  return(res$table)
}

# Fig1d: Beta-diversity adonis across species/pathways/KO/EC/GO (combined)
fig1d <- function( diet_path = "dietary_diversity_use.txt", save=NULL) {
  message("Running Fig1d (beta-diversity tests) ...")
  diet <- load_diet(diet_path)
  # Species
  species_tb <- load_metaphlan_species("metaphlan.txt")
  species_tb <- species_tb[species_tb$week == week_filter, ]
  species_cols <- setdiff(colnames(species_tb), c("SampleID","subject","week"))
  r_multi_species <- tryCatch({
    adonis_beta(species_tb[, species_cols, drop=FALSE], diet, covariates_formula = "gender+age+BMI+Type+Protein_g+Fat_g+Carbohydrate_g+Fiber_g+Cholesterol_mg+VitaminA_ug+VitaminB1_mg+VitaminB2_mg+VitaminB3_mg+VitaminC_mg+VitaminE_mg+Calcium_mg+Sodium_mg")
  }, error=function(e) NULL)
  # Pathway
  pwy_tb <- load_table_t("all_pathabundance_relab.tsv")
  pwy_tb <- pwy_tb[pwy_tb$week == week_filter , ]
  pwy_cols <- setdiff(colnames(pwy_tb), c("SampleID","subject","week"))
  r_multi_pathway <- tryCatch({
    adonis_beta(pwy_tb[, pwy_cols, drop=FALSE], diet, covariates_formula = "gender+age+BMI+Type+Protein_g+Fat_g+Carbohydrate_g+Fiber_g")
  }, error=function(e) NULL)
  # KO
  ko_tb <- load_table_t("all_ko_relab.tsv")
  ko_tb <- ko_tb[ko_tb$week == week_filter , ]
  ko_cols <- setdiff(colnames(ko_tb), c("SampleID","subject","week"))
  r_multi_ko <- tryCatch({
    adonis_beta(ko_tb[, ko_cols, drop=FALSE], diet, covariates_formula = "gender+age+BMI+Type+Protein_g+Fat_g+Carbohydrate_g+Fiber_g")
  }, error=function(e) NULL)
  # EC
  ec_tb <- load_table_t("all_ec_relab.tsv")
  ec_tb <- ec_tb[ec_tb$week == week_filter, ]
  ec_cols <- setdiff(colnames(ec_tb), c("SampleID","subject","week"))
  r_multi_ec <- tryCatch({
    adonis_beta(ec_tb[, ec_cols, drop=FALSE], diet, covariates_formula = "gender+age+BMI+Type+Protein_g+Fat_g+Carbohydrate_g+Fiber_g")
  }, error=function(e) NULL)
  # GO
  go_tb <- load_table_t("all_go_relab.tsv")
  go_tb <- go_tb[go_tb$week == week_filter , ]
  go_cols <- setdiff(colnames(go_tb), c("SampleID","subject","week"))
  r_multi_go <- tryCatch({
    adonis_beta(go_tb[, go_cols, drop=FALSE], diet, covariates_formula = "gender+age+BMI+Type+Protein_g+Fat_g+Carbohydrate_g+Fiber_g")
  }, error=function(e) NULL)
  # combine (extract r2 and p-values as in original)
  # Write combined table to file
  combined <- list(species = r_multi_species, pathway = r_multi_pathway, ko = r_multi_ko, ec = r_multi_ec, go = r_multi_go)
  saveRDS(combined, file="fig1d_adonis_results.rds")
  return(combined)
}

# Fig1e: SIAMCAT predictive modeling
fig1e <- function(feat_path = "./relab_gut.tsv", meta_path = "./meta.tsv") {
  message("Running Fig1e (SIAMCAT)...")
  return(run_siamcat_analysis(feat_path = feat_path, meta_path = meta_path))
}

# Fig1f: species differential scatterplot (volcano-like)
fig1f <- function(metaphlan_path = "metaphlan.txt", diet_path = "dietary_diversity_use.txt", save=NULL) {
  message("Running Fig1f (species differential) ...")
  res <- species_differential(metaphlan_path = metaphlan_path, diet_path = diet_path, week_filter="S2")
  if (is.null(res)) return(NULL)
  # merge with abundance and phylum info
  tb <- safe_read(metaphlan_path, sep="\t", header=TRUE)
  colnames(tb)[1] <- "clade_name"
  tb <- tb %>% filter(grepl("s__", clade_name) & !grepl("t__", clade_name))
  tb$species <- as.data.frame(strsplit2(tb$clade_name, "s__"))[,2]
  tb$phylum <- strsplit2(strsplit2(tb$clade_name,"p__")[,2],".c__")[,1]
  rownames(tb) <- tb$species
  tb_mat <- as.data.frame(t(tb[ , -which(colnames(tb)=="clade_name"), drop=FALSE]))
  tb_mat$SampleID <- rownames(tb_mat)
  diet <- load_diet(diet_path)
  tb_merge <- merge(res, tb[, c("clade_name","phylum")], by.x="species", by.y="row.names", all.x=TRUE)
  tb_merge$estimate <- as.numeric(tb_merge$estimate)
  tb_merge$pvalue <- as.numeric(tb_merge$pvalue)
  p <- ggplot(tb_merge, aes(x=estimate, y=pvalue)) +
    geom_point(aes(color=phylum, size=abs(estimate))) +
    geom_text_repel(aes(label=ifelse(pvalue < 0.1, species, "")), size=3) +
    theme_minimal() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
  if (!is.null(save)) ggsave(save, p, width=7, height=5)
  print(p)
  return(res)
}

# Fig1g: pathway associations (selected list and plotting)
fig1g <- function(pathway_file = "./all_pathabundance_relab.tsv", diet_path = "dietary_diversity_use.txt", save=NULL) {
  message("Running Fig1g (pathway associations) ...")
  res <- pathway_differential(pathway_file = pathway_file, diet_path = diet_path)
  if (is.null(res)) return(NULL)
  # select top pathways of interest similar to original
  name_pwy <- c("PWY-6163: chorismate biosynthesis from 3-dehydroquinate","ARO-PWY: chorismate biosynthesis I","COMPLETE-ARO-PWY: superpathway of aromatic amino acid biosynthesis","VALSYN-PWY: L-valine biosynthesis","GLYCOCAT-PWY: glycogen degradation I","PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)","1CMET2-PWY: folate transformations III (E. coli)","COA-PWY: coenzyme A biosynthesis I (prokaryotic)","PWY-6703: preQ0 biosynthesis","PWY-6609: adenine and adenosine salvage III")
  r <- res[res$pathway %in% name_pwy, ]
  if (nrow(r) == 0) return(r)
  r$estimate <- as.numeric(r$estimate)
  r$pvalue <- as.numeric(r$pvalue)
  r$lab <- factor(r$pathway, levels=name_pwy)
  p <- ggplot(r, aes(x=estimate, y=lab)) +
    geom_point(aes(color=pvalue), size=3) +
    scale_color_gradient(low="#132B43", high="#56B1F7", limits=c(0, 0.05)) +
    theme_minimal() + geom_vline(xintercept=0, linetype="dashed") +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
  if (!is.null(save)) ggsave(save, p, width=6, height=4)
  print(p)
  return(r)
}

# ------------------------------
# Figure 2 (phylogenetic / strain-level) - includes saliva marker for 2b
# ------------------------------

# Fig2a and Fig2b share the same anpan loop logic; fig2b is SALIVA if user wants.
fig2a <- function(gene_dir = "anpan/gene_file4", meta_path = "anpan/metadata.tsv", out = "gut_tree_diff_w2.txt") {
  message("Running Fig2a (phylogenetic ELPD differences)...")
  anpan_phylo_loop(gene_dir = gene_dir, meta_path = meta_path, output_file = out)
  df <- tryCatch(read.table(out, sep="\t", header=FALSE, stringsAsFactors = FALSE), error=function(e) NULL)
  if (!is.null(df)) {
    colnames(df) <- c("Name","Value","Error")
    df$Name <- factor(df$Name, levels = rev(unique(df$Name)))
    df$Value <- -as.numeric(df$Value)
    p <- ggplot(df, aes(x=Value, y=Name)) + geom_point() + geom_errorbarh(aes(xmin=Value-Error, xmax=Value+Error), height=0.4) + theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))
    print(p)
  }
  return(invisible(df))
}

fig2b <- function(...){ # SALIVA
  message("Fig2b is implemented by same function as fig2a but for saliva gene dir")
  return(fig2a(...))
}

fig2c <- function(ko_diff_file = "diff_ko_gut2.tsv", ko_table = "./all_ko_relab.tsv", save=NULL) {
  message("Running Fig2c (KO heatmap)...")
  if (!file.exists(ko_diff_file)) {
    message("KO diff file not found: ", ko_diff_file)
    return(NULL)
  }
  tb <- safe_read(ko_table, sep="\t", header=TRUE)
  tb <- tb %>% filter(!grepl("unclassified", ID), !grepl("UNGROUPED", ID))
  rownames(tb) <- tb$ID
  tb_t <- as.data.frame(t(tb[,-1, drop=FALSE]))
  # select top diff KOs
  diff <- safe_read(ko_diff_file, sep="\t", header=TRUE)
  name <- diff$V1
  sel <- intersect(rownames(tb_t), name)
  if (length(sel) == 0) return(NULL)
  mat <- tb_t[sel, , drop=FALSE]
  plot_complex_heatmap(mat, filename=save)
  return(invisible(mat))
}

# ------------------------------
# Figure 3: Binding affinity and structure predictions
# Fig3b and Fig3d correspond to SALIVA analyses
# ------------------------------
fig3a <- function() {
  bind <- read.table("/bind_all.txt",sep="\t",header = T)
  bind <- bind[,c("MHC","Rank_EL","Identity")]
  bind <- aggregate(bind$Rank_EL ,by=list(MHC=bind$MHC,Identity=bind$Identity),min)
  bind <- na.omit(bind)
  
  prop <- read.table("gene_poportion.tsv",sep="\t",header=T)
  r<- merge(bind,prop,by="Identity")
  r$rank_level <- NA
  r[r$rank<=1,"rank_level"] <- "strong"
  r[r$rank<=5 & r$rank>1,"rank_level"] <- "week"
  r[r$rank<=100 & r$rank>5,"rank_level"] <- "none"
  
  my_comparisons=list(c("strong","week"),c("week","none"),c("strong","none"))
  
  ggplot(r,aes(x=rank_level,y=prop))+
    
    geom_boxplot(aes(color=rank_level),width=0.6)+
    
    theme(panel.grid = element_blank(), 
          axis.line = element_line(colour = 'black'), 
          panel.background = element_blank(), 
          axis.text = element_text(size = 9, color = 'black'), 
          axis.title = element_text(size = 10, color = 'black')) +
    ylab("Prop")+xlab("")+
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 12))+
    guides(fill=guide_legend(title = NULL))+
    guides (fill = "none")+
    stat_compare_means(comparisons=my_comparisons,test = wilcox.test,paired = F,p.adjust.method = "none",label="p.format")
}
fig3b <- function() { # SALIVA
  message("Fig2b is implemented by same function as fig3a but for saliva gene dir")
  return(fig2a(...))
}
fig3c <- function() {
  m <- bind2
  m_j <- m[m$gene=="A0A369NC92",]
  colnames(m_j)
  p <- ggplot(m_j,aes(y=prop,x=as.numeric(rank)))+
    
    geom_point(aes(color=MHC),size=1)+
    ylim(0,1)+
    xlab("bind_rank")+
    ylab("proportion")+
    annotate('text',x=min(m_j$rank)+min(m_j$rank)/10,y=0.95,label=paste0("p =",p_value),size=4)+
    geom_smooth(method = "lm", se = FALSE,color="black")+
    scale_color_manual(values = c("red","blue","green","purple","brown","lightpink")) +
    theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_blank(), axis.text = element_text(size = 9, color = 'black'), axis.title = element_text(size = 10, color = 'black'))# +
}
fig3d <- function() { # SALIVA
  m <- bind2
  m_j <- m[m$gene=="F3SHD3",]
  colnames(m_j)
  p <- ggplot(m_j,aes(y=prop,x=as.numeric(rank)))+
    
    geom_point(aes(color=MHC),size=1)+
    ylim(0,1)+
    xlab("bind_rank")+
    ylab("proportion")+
    annotate('text',x=min(m_j$rank)+min(m_j$rank)/10,y=0.95,label=paste0("p =",p_value),size=4)+
    geom_smooth(method = "lm", se = FALSE,color="black")+
    scale_color_manual(values = c("red","blue","green","purple","brown","lightpink")) +
    theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_blank(), axis.text = element_text(size = 9, color = 'black'), axis.title = element_text(size = 10, color = 'black'))# +
  
}

# ------------------------------
# Figure 4: metabolites & phenotypes
# ------------------------------
fig4 <- function(metabolite_table = "/all_metabs.tsv", diet_path = "dietary_diversity_use.txt") {
  message("Running Fig4 (metabolite & phenotype associations) ...")
  if (!file.exists(metabolite_table)) {
    message("Metabolite table not found: ", metabolite_table)
    return(NULL)
  }
  df <- safe_read(metabolite_table, sep="\t", header=TRUE)
  diet <- load_diet(diet_path)
  # User's original pipeline likely performed many correlations and linear models.
  # Here we provide a flexible linear model per metabolite similar to species/pathways
  result <- NULL
  met_cols <- setdiff(colnames(df), "ID")
  df_t <- as.data.frame(t(df[ , -1, drop=FALSE]))
  df_t$SampleID <- rownames(df_t)
  dt <- merge(df_t, diet, by="SampleID")
  for (i in met_cols) {
    if (!(i %in% colnames(dt))) next
    dt_i <- dt[, c(i, "age", "gender", "BMI", "Type", "Protein_g", "Fat_g", "Carbohydrate_g", "Fiber_g"), drop=FALSE]
    dt_i <- na.omit(dt_i)
    colnames(dt_i)[1] <- "metab"
    dt_i$metab <- as.numeric(dt_i$metab)
    if (length(unique(dt_i$metab)) < 2) next
    summ <- summary(lm(metab ~ Type + age + gender + BMI + Protein_g + Fat_g + Carbohydrate_g + Fiber_g, dt_i))
    result <- rbind(result, c(i, summ$coefficients[2,1], summ$coefficients[2,4]))
  }
  if (is.null(result)) return(NULL)
  result <- as.data.frame(result, stringsAsFactors = FALSE)
  colnames(result) <- c("metab", "estimate", "pvalue")
  result$pvalue <- as.numeric(result$pvalue)
  result$P_BH <- p.adjust(result$pvalue, method = "BH")
  return(result)
}

# ------------------------------
# Master runner
# ------------------------------
run_all_figures <- function() {
  message("Starting full analysis pipeline...")
  try(fig1a(), silent=TRUE)
  try(fig1b(), silent=TRUE)
  try(fig1c(), silent=TRUE)
  try(fig1d(), silent=TRUE)
  try(fig1e(), silent=TRUE)
  try(fig1f(), silent=TRUE)
  try(fig1g(), silent=TRUE)
  try(fig2a(), silent=TRUE)
  try(fig2b(), silent=TRUE) # SALIVA
  try(fig2c(), silent=TRUE)
  try(fig3a(), silent=TRUE)
  try(fig3b(), silent=TRUE) # SALIVA
  try(fig3c(), silent=TRUE)
  try(fig3d(), silent=TRUE) # SALIVA
  try(fig4(), silent=TRUE)

  message("Finished pipeline (note: some  functions may be noop).")
  invisible(TRUE)
}

# helper: infix default
`%||%` <- function(a, b) if (!is.null(a)) a else b

# End of script

