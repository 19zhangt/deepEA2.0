#! /usr/bin/Rscript

# Setting environment variables
options(stringsAsFactors = F)
setwd('/storage_server/home/zhangt/freePEA/mlPEA_denovo/results')

# Loading Packages
source("part1_denovo_transcripts/functions.R")

genome_dir <- "/distorage_server/home/zhangt/genomes"
data_dir <- "/storage_server/home/zhangt/freePEA/data"
script_dir <- "/storage_server/home/zhangt/freePEA/scripts"
sp_all <- c("Ath", "Gar", "Ghi", "Pvu", "Gma", "Sbi", "Zma", "Ata", "Tdi", "Tae", "Osa", "Ppa")


########################################
## 1.1 Summary of the number of assembled transcripts 
########################################
summary_table <- list()
for(sp in sp_all){ summary_table[[sp]] <- fTransfragsCount(i = sp) }
summary_table_df <- do.call('rbind', summary_table)
write.table(summary_table_df, file = "part1_denovo_transcripts/results_table/transcripts_count.txt", 
            quote = F, row.names = F, col.names = F, sep = "\t")


########################################
## 1.2 Percentage
########################################
run_sp <- "Ath"
tree_info <- fNumFlow(df = summary_table[[run_sp]]); plot(tree_info)

## Percentage of different types of classification
tool_list <- c("TransABySS", "rnaSPAdes", "Trinity")
r_level <- c("known", "long", "short", "unmapped", "sc_chimera", "dc_chimera", "other")
radar_plots <- list()
for(tool in tool_list){
    hc_mapping_info <- fSeq2Genome(i = run_sp, fasta_file = sprintf("%s/%s/02_merge_assembly/rep/HC_%s.fasta",
                                                                 data_dir, run_sp, tool))
    lc_mapping_info <- fSeq2Genome(i = run_sp, fasta_file = sprintf("%s/%s/02_merge_assembly/rep/LC_%s.fasta",
                                                                    data_dir, run_sp, tool))
    hc_val <- table(hc_mapping_info)/sum(table(hc_mapping_info))
    lc_val <- table(lc_mapping_info)/sum(table(lc_mapping_info))
    hc_val <- hc_val[r_level]
    lc_val <- lc_val[r_level]
    t_df <- do.call('rbind', list(c("HC", hc_val), c("LC", lc_val)))
    class(t_df) <- "numeric"
    t_df <- as.data.frame(t_df)
    t_df[, 1] <- c("HC", "LC")
    # devtools::install_github("ricardo-bion/ggradar")
    cat(sprintf("%s_%s", run_sp, tool), "\n")
    radar_plots[[sprintf("%s_%s", run_sp, tool)]] <- fRadarPlot(t_df)
}

pdf("part1_denovo_transcripts/results_plot/part1_pipeflow_radar.pdf", width = 15, height = 6)
patchwork::wrap_plots(radar_plots, ncol = 3, guides = "collect")
dev.off()
names(radar_plots)

tool_count <- fToolCount(i = run_sp)
ggsave("part1_denovo_transcripts/results_plot/part1_tool_overlap.pdf", width = 7, height = 4)


########################################
## 1.3 ML learning
########################################
# BUSCO doesn't contain three species: "Gar", "Ghi", "Tdi" 
for(tmp_sp in sp_all){
  ml_result <- fMLpipe(i = tmp_sp)
  save(ml_result, file = sprintf("part1_denovo_transcripts/RData/%s_prediction.RData", tmp_sp))
  cat(tmp_sp, sum(table(ml_result$df$ml_res)), sum(table(ml_result$df$ml_res)[c("busco_names", "ml_remain")]), "\n")
}



mc_val <- table(mapping_info)/sum(table(mapping_info))
mc_val <- mc_val[r_level]
class(mc_df) <- "numeric"
mc_df <- as.data.frame(mc_df)
mc_df[, 1] <- c("Positiver")
fRadarPlot(mc_df)


hc_mapping_info <- mapping_info[ml_result_list[[run_sp]]$pu_remain]
lc_mapping_info <- mapping_info[ml_result_list[[run_sp]]$pu_remove]
hc_val <- table(hc_mapping_info)/sum(table(hc_mapping_info))
lc_val <- table(lc_mapping_info)/sum(table(lc_mapping_info))
hc_val <- hc_val[r_level]
lc_val <- lc_val[r_level]
t_df <- do.call('rbind', list(c("HC", hc_val), c("LC", lc_val)))
class(t_df) <- "numeric"
t_df <- as.data.frame(t_df)
t_df[, 1] <- c("HC", "LC")

add_radar_plots <- radar_plots
add_radar_plots[['PU']] <- fRadarPlot(t_df, f_col = c("#E04B87", "#BBB3BC"))


pdf("part1_denovo_transcripts/results_plot/part1_radar_plots.pdf", width = 10, height = 8)
patchwork::wrap_plots(add_radar_plots, ncol = 2, byrow = T, guides = "collect")
dev.off()

eval_list <- list()
for(sp in names(ml_result_list)){
  eval_list[[sp]] <- ml_result_list[[sp]]$ml
}

eval_val <-  do.call('rbind', eval_list)
write.table(eval_val, file = "part1_denovo_transcripts/results_table/evaluation_PUlearning.txt", 
            quote = F, row.names = T, col.names = T, sep = "\t")


AUC_val <- eval_val[, "AUC"]
AUC_df <- data.frame('sp'= names(AUC_val), 'value'= as.numeric(AUC_val))

AUC_df %>% 
  mutate(sp = factor(sp, levels = sp_all)) %>% 
  ggplot(aes(x=sp, y=value, fill=sp)) +
  geom_bar(stat = "identity", position = 'dodge', color=NA, size=1/4) +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c("#F5BC2D", "#E7916E", "#058085", "#A72E5B", "#FDE6BA",
                               "#E62787", "#DAACD1", "#208738", "#F5D01D", "#DCDCDC",
                               "#1490C3", "#854383")) +
  coord_cartesian(ylim = c(0.5, 1)) + theme_bw() + BaseTheme() + ylab("AUC") + xlab("")
ggsave("part1_denovo_transcripts/results_plot/part1_AUC.pdf", width = 6.5, height = 3.5)


## Cross Species
for (x in sp_all){
  for (y in sp_all){
    if ( x != y ){
      f_tmp <- sprintf("%s_%s", x, y)
      f_outfile <- sprintf("./part1_denovo_transcripts/ml_out/%s.txt", f_tmp)
      if (!file.exists(f_outfile)){
        print(f_tmp)
        cross_species_pred[[f_tmp]] <- fCrossSpeciesPred(x = c(x, y))
        write.table(cross_species_pred[[f_tmp]],
                    file = f_outfile, quote = F, sep = "\t", col.names = F)
      }
    }
  }
}

save(cross_species_pred, file = "part1_denovo_transcripts/RData/20220122_corss_prediction.RData")

prediction_mat <- matrix(nrow = length(sp_all), ncol = length(sp_all))
rownames(prediction_mat) <- colnames(prediction_mat) <- sp_all
for (x in sp_all){
  for (y in sp_all){
    if ( x != y ){
      f_tmp <- sprintf("%s_%s", x, y)
      prediction_mat[x, y] <- as.numeric(cross_species_pred[[f_tmp]]['AUC'])
    } else {
      prediction_mat[x, y] <- NA
    }
  }
}

pheatmap::pheatmap(prediction_mat, cluster_rows = F, cluster_cols = F)



# cl <- makeCluster(2); registerDoParallel(cl)
# out <- foreach(nn = seq_len(nrow(combind_vec)), .combine = list,.multicombine = TRUE) %dopar% {
#   x <- combind_vec[nn, ]
#   fCrossSpeciesPred(x = x)
# }
# stopImplicitCluster()



# ## Remove Chimera
# length(intersect(pred.prediction, busco_names))
# merge_transfrags <- c(pred.prediction, busco_names)
# 
# out_dir <- sprintf("%s/%s/03_evaluation/PU_filter", data_dir, sp)
# out_file <- paste0(out_dir, "/remain.fasta")
# out_chimera <- sprintf("%s/chimera.txt", out_dir)
# 
# if(!file.exists(out_chimera)){
#   system(command = paste0("mkdir -p ", out_dir))
#   seqinr::write.fasta(sequences = f_transfrags[merge_transfrags],
#                       names = merge_transfrags,
#                       file.out = out_file)
#   
#   command_faindex <- sprintf("samtools faidx %s", out_file)
#   tool_path <- "~/miniconda3/envs/denovota/bin"
#   command_bed <- sprintf("%s/bedtools makewindows -g %s.fai -w 60 > %s/window.bed", 
#                          tool_path, out_file, out_dir)
#   command_windows <- sprintf("%s/bedtools getfasta -bed %s/window.bed -fi %s -fo %s/window.fa", 
#                              tool_path, out_dir, out_file, out_dir)
#   command_mapindex <- sprintf("~/publicdir/software/hisat2-2.2.1/hisat2-build-s %s %s_index", 
#                               out_file, out_file)
#   hisat_command <- "~/publicdir/software/hisat2-2.2.1/hisat2 --norc --no-spliced-alignment --score-min L,0,-0 -k 20 -p 40 -f "
#   command_align <- sprintf("%s -x %s_index -U %s/window.fa | samtools view -F 4 -Sb - | samtools sort -o %s/windows.bam -", 
#                            hisat_command, out_file, out_dir, out_dir)
#   # a_R3_TRINITY_DN514_c0_g1_i2
#   command_bamindex <- sprintf("samtools index %s/windows.bam", out_dir)
#   system(command = command_faindex)
#   system(command = command_bed)
#   system(command = command_windows)
#   if(!file.exists(sprintf("%s_index.1.ht2", out_file))){
#     system(command = command_mapindex)
#   }
#   system(command = command_align)
#   system(command = command_bamindex)
#   
#   py_chimera_iden = sprintf("~/miniconda3/envs/metaPlants/bin/python %s/chimera.py %s/windows.bam  %s/chimera.txt", script_dir, out_dir, out_dir)
#   
#   system(command = py_chimera_iden)
# }
# 
# predict_chimera <- read.table(sprintf("%s/chimera.txt", out_dir))
# dim(predict_chimera)
# 
# fPieRatio(f_itable = table(seq_type[predict_chimera$V1]), main = "predict_chimera")
# length(names(seq_type)[seq_type=="dc_chimera"])
# length(intersect(names(seq_type)[seq_type=="dc_chimera"], predict_chimera$V1))
# 
# library(Biostrings)
# pairwiseAlignment(pattern = f_transfrags[["a_R3_TRINITY_DN962_c0_g1_i8"]], 
#                   subject = f_transfrags[["a_R3_TRINITY_DN962_c0_g1_i8"]],
#                   type = "local")
# 
# a <- read.table("/storage_server/home/zhangt/freePEA/data/Zma/03_evaluation/PU_filter/tmp/tmp.txt")
# fPieRatio(f_itable = table(seq_type[a$V1]), main = "predict_chimera")
# 
# 
# grid::grid.newpage()
# grid::grid.draw(VennDiagram::venn.diagram(filename = NULL, 
#                                           list("predict"=predict_chimera$V1,
#                                                "Chimeric2"=names(seq_type)[seq_type==""]),
#                                           disable.logging = TRUE))
# 
# cat(length(transfrags), length(pu_remain),
#     length(setdiff(pu_remain, predict_chimera$V1)))
# 
# 
# candidate_transfrags <- setdiff(pu_remain, predict_chimera$V1)
# fPieRatio(f_itable = table(contig.tree$class[candidate_transfrags]), 
#           main = "candidate_transfrags")
# fPieRatio(f_itable = table(compare_gff_mat[rownames(compare_gff_mat)%in%candidate_transfrags, 4]), 
#           main = "candidate_transfrags")
# 
# setdiff(remove_list$Chimeric2, predict_chimera$V1)[1:10]
# setdiff(predict_chimera$V1, c(remove_list$Chimeric1,remove_list$Chimeric2))[1:10]

# for unmapped sequences ####
# # Upset plot #### Ghi Gma Zma Ata Tdi Tae Ppa Gar Pvu Sbi Osa
# diff_sp <- list()
# for (y in c("Gar", "Ghi", "Pvu", "Gma", "Sbi", "Zma", "Osa", "Ppa")){
#   f_sp_ctgs <- read.table(paste0("data/Ath/02_merge_assembly/filter/species_alignment/", y, "_align.txt"), 
#                           sep = ",", stringsAsFactors = F, header = T)
#   f_sp_ctgs <- f_sp_ctgs[f_sp_ctgs$coverage>80&f_sp_ctgs$identity>80, ]
#   f_sp_ctgs$sp <- y
#   diff_sp[[y]] <- f_sp_ctgs
# }
# diff_sp_mat <- do.call('rbind', diff_sp)
# diff_sp_mat <- reshape2::dcast(diff_sp_mat, eli.qname~sp)
# rownames(diff_sp_mat) <- diff_sp_mat[,1]
# diff_sp_mat <- diff_sp_mat[,-1]
# diff_sp_mat[diff_sp_mat!=0] <- TRUE
# diff_sp_mat[diff_sp_mat==0] <- FALSE
# ComplexUpset::upset(diff_sp_mat, colnames(diff_sp_mat), 
#                     name='', width_ratio=0.1, min_size=40)


# diff_supp_num_list <- tree_class_list <- list()
# for(j in 1:4){
#   t_count_range <- vector()
#   for (t_num in which(cluster_count == j)) {
#     t_count_range <- c(t_count_range, cluster_start[t_num]:cluster_end[t_num])
#   }
#   # extract names
#   t_extract_tf <- apply(cluster_df[t_count_range, ], 1, function(x){strsplit(x[2], " >| ")[[1]][2]})
#   t_extract_tf <- gsub("\\.\\.\\.$", "", t_extract_tf)
#   t_extract_names <- t_extract_tf[t_extract_tf %in% names(transfrags)]
#   diff_supp_num_list[[j]] <- t_extract_names
#   # Decision tree
#   t_tree <- fBinaryTree(input_df = cigar_df, true_contigs = t_extract_names,
#                           remove_list = list("Chimeric1"=intersect(ct_df$name[ct_df$cn_num==1],t_extract_names),
#                                              "Chimeric2"=intersect(ct_df$name[ct_df$cn_num>1],t_extract_names),
#                                              "Unmapped"=intersect(unique(unmapped_transfrags),t_extract_names)),
#                           intersect_file = map_file$intersect_txt)
#   tree_class_list[[j]] <- t_tree$class
#   htmlwidgets::saveWidget(plot(t_tree$tree), paste0("Part_1/figures/", i, "software_support_", j, ".html"))
# }
# 
# fPieRatio(f_itable = table(tree_class_list[[1]]), main = "s1_class")
# fPieRatio(f_itable = table(tree_class_list[[2]]), main = "s2_class")
# fPieRatio(f_itable = table(tree_class_list[[3]]), main = "s3_class")
# fPieRatio(f_itable = table(tree_class_list[[4]]), main = "s4_class")

## 
# hc_class <- c("Intergenic", "Genic", "LC", "Unperfect", "Chimera")
# aligned_pool <- names(contig.tree$class)[contig.tree$class%in%hc_class]
# unaligned_pool <- names(contig.tree$class)[!contig.tree$class%in%hc_class]

# psol.feature.df <- feature.filter.df[c(i.busco.ctgs, supp.one.ctgs), ]
# psolRes.save <- PEA::PSOL(featureMatrix = psol.feature.df, 
#                           positives = i.busco.ctgs,
#                           unlabels = supp.one.ctgs,
#                           PSOLResDic = tmpPSOLdic,
#                           cpus = 40, negNum = 2000 )
# save.image(file = "1204.RData")
# load("1204.RData")
# psol.feature.tmp.df <- feature.filter.df
# extract feature ####
# trasrate.features <- c("prop_gc", "orf_length", "sCnuc", "sCcov", "sCord", "sCseg")
# detonate.features <- c("length", "effective_length", "FPKM", "contig_impact_score")
# 
# feature_list <- busco.list <- list()
# for(tool.name in TOOLNAMES){
#   for(rep.name in REPNAMES){
#     # import results (trasrate and detonate)
#     EVALUDIR <- paste0("data/", rep.name, "/02evaluation/")
#     i.trasrate <- read.table(paste0(EVALUDIR, "TransRate/", tool.name, "/", tool.name, "/contigs.csv"), 
#                              sep = ",", header = T, row.names = 1)
#     i.detonate <- read.table(paste0(EVALUDIR, "detonate/", tool.name,"/rsem_score.score.isoforms.results"),
#                              sep = "\t", header = T, row.names = 1)
#     i.detonate[, "contig_impact_score"] <- i.detonate[, "contig_impact_score"] + 
#       abs(min(i.detonate[, "contig_impact_score"])) + 1
#     i.featureMat <- cbind(i.detonate[, detonate.features],
#                           i.trasrate[rownames(i.detonate), trasrate.features])
#     rownames(i.featureMat) <- rownames(i.detonate)
#     i.featureMat[is.na(i.featureMat)] <- 0
#     # scaling some values
#     for(iii in c("length", "effective_length", "orf_length", "contig_impact_score", "FPKM")){
#       i.featureMat[, iii] <- log2(i.featureMat[,iii]+1)
#       colnames(i.featureMat)[colnames(i.featureMat)==iii] <- paste(iii, "_log2")
#     }
#     # rename
#     abbr.tool <- letters[which(TOOLNAMES == tool.name)]
#     abbr.rep <- strsplit(rep.name, "_")[[1]][2]
#     rownames(i.featureMat) <- paste0(abbr.tool, "_", abbr.rep, "_", rownames(i.featureMat))
#     i.featureMat$name <- rownames(i.featureMat)
#     feature_list[[paste0(tool.name, "_", rep.name)]] <- i.featureMat
#     cat(paste0(tool.name, "_", rep.name), "\n")
#     # import busco
#     i.busco <- read.delim2(paste0(EVALUDIR, "BUSCO/", tool.name, 
#                                   "/run_brassicales_odb10/full_table.tsv"), sep = "\t", header = F)
#     i.busco.complete <- i.busco[i.busco[,2]=="Complete", ]
#     i.busco.complete.names <- unique(apply(i.busco.complete, 1, function(x){
#       strsplit(x[3], split = ":")[[1]][1]}))
#     busco.list[[paste0(tool.name, "_", rep.name)]] <- paste0(abbr.tool, "_", abbr.rep, "_", i.busco.complete.names)
#   }
# }

# feature plot in different types ####
# tmp.class <- rep(NA, nrow(fetMat.filter))
# names(tmp.class) <- rownames(fetMat.filter)
# for(ii in 1:4) tmp.class[diff_supp_num_list[[ii]]] <- ii
# 
# fetMat.filter$type <- factor(tmp.class, levels = sort(unique(tmp.class)))
# fetMat.filter$name <- rownames(fetMat.filter)
# fetMat.filter.df <- reshape2::melt(fetMat.filter, id=c("type","name"))
# 
# pdf("features_in_different_tool_support.pdf", height = 8, width = 12.5)
# ggplot(fetMat.filter.df, aes(x=value, color=type)) +
#   geom_line(stat = "density", size=0.9, alpha=0.8) + facet_wrap(variable~., scales = "free", ncol = 4) +
#   theme_bw(base_size = 14) +
#   scale_color_manual(values = c("#8BC33E", "#00ACED", "#F5911F", "#EB1B24")) +
#   BaseTheme()
# dev.off()
# 
# tmp.class2 <- rep(NA, nrow(fetMat.filter))
# names(tmp.class2) <- rownames(fetMat.filter)
# for(ii in 1:4) {
#   tmp.class2[tree_class_list[[ii]]$align] <- "align"
#   tmp.class2[tree_class_list[[ii]]$unalign] <- "unalign"
# }
# fetMat.filter$type2 <- factor(tmp.class2, levels = sort(unique(tmp.class2)))
# fetMat.filter.df2 <- reshape2::melt(fetMat.filter, id=c("type","type2","name"))
# 
# pdf("features2_in_different_tool_support.pdf", height = 8, width = 12.5)
# ggplot(fetMat.filter.df2, aes(x=value, color=type2)) +
#   geom_line(stat = "density", size=0.9, alpha=0.8) + facet_wrap(variable~., scales = "free", ncol = 4) +
#   theme_bw(base_size = 14) + scale_color_manual(values = c("#EB1B24", "#8BC33E")) +
#   BaseTheme()
# dev.off()
# 
# pdf("features3_in_different_tool_support.pdf", height = 8, width = 12.5)
# ggplot(fetMat.filter.df2, aes(x=value, color=paste0(type, type2))) +
#   geom_line(stat = "density", size=0.9) + facet_wrap(variable~., scales = "free", ncol = 4) +
#   theme_bw(base_size = 14) + scale_color_manual(values = c("#8BC33E", "#45B753", "#00ACED", "#00C5EC",
#                                                            "#F5911F", "#B99600", "#EB1B24", "#E09A8E")) +
#   BaseTheme()
# dev.off()

############################## Application of ML (PSoL)  ########################################
# extract negative contigs ####
# tmpPSOLdic <- "RData/tmp_psol_61446/"
# system(command = paste0("mkdir ",tmpPSOLdic ))
# pos.names <- diff_supp_num_list[[4]]
# unlabel.names <- setdiff(rownames(fetMat.filter), pos.names)
# 
# psolRes.save <- PEA::PSOL(featureMatrix = fetMat.filter, positives = pos.names,
#                           unlabels = unlabel.names,
#                           PSOLResDic = tmpPSOLdic,
#                           cpus = 30 )
# 
# load("RData/tmp_psol_61446/PSOL_InitialNegativeSelection_Res.RData")


# # using BUSCO results to find negative
# pos.names <- busco.filter
# unlabel.names <- setdiff(diff_supp_num_list[[1]], pos.names)
# psol.fea.mat <- fetMat.filter[c(pos.names, unlabel.names), ]
# system(paste0("mkdir ./RData/tmp_psol_24478_5features/"))
# # busco.findPSoL
# 
# busco.6feat.findPSoL <- PEA::PSOL(featureMatrix = psol.fea.mat[, 1:6], positives = pos.names,
#                                   unlabels = unlabel.names, PSOLResDic = "./RData/tmp_psol_24478_6features/",
#                                   cpus = 40 , negNum = 2000)
# load("RData/tmp_psol_24478_6features/PSOL_InitialNegativeSelection_Res.RData")
# 
# FeaturePlot4PU(fea.mat.tmp = psol.fea.mat[, 1:6], pos.tmp = res$positives, 
#                neg.tmp = res$negatives, tree.list = tree_class_list)
# 
# # tmp.psolRes <- .PSOL_NegativeExpansion(featureMat = fetMat.filter, positives = pos.names, 
# #                                        unlabels = setdiff(rownames(fetMat.filter), c(pos.names, res$negatives)),negatives = res$negatives, TPR = 0.9999, PSOLResDic = tmpPSOLdic, cpus = 30)
# load("RData/tmp_psol_24478/PSOL_InitialNegativeSelection_Res.RData")
# 
# FeaturePlot4PU(fea.mat.tmp = fetMat.filter, pos.tmp = res$positives, 
#                neg.tmp = res$negatives, tree.list = tree_class_list)
