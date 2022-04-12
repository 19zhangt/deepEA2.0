#! /usr/bin/Rscript

# Setting environment variables
options(stringsAsFactors = F)
setwd('/storage_server/home/zhangt/freePEA/mlPEA_denovo/results')

# Loading Packages
library(dplyr)
library(ggplot2)
library(data.table)
library(future.apply)

data_dir <- "/storage_server/home/zhangt/freePEA/data"
genome_dir <- "/distorage_server/home/zhangt/genomes"
script_dir <- "/storage_server/home/zhangt/freePEA/scripts"
source("part2_peak_calling/functions.R")

########################################
## Peak calling
########################################
for(run_sp in c("Ath", "Osa", "Ppa")){
  ## call peak in three replicates and mergeing them
  peak_res <- sprintf("part2_peak_calling/RData/%s_peak.RData", run_sp)
  if(file.exists(peak_res)){
    load(peak_res)
  }else{
    system(command = "mkdir -p part2_peak_calling/RData")
    peaks <- fPEATransciptPeak(input_sp = run_sp)
    save(peaks, file = peak_res)
    # peaks2 <- fCallRelativePeak(input_sp = "Ath")
    # peaks2 <- as.data.frame(ath_peaks2)
  }
  peak_grangs <- eval(parse(text = sprintf("%s_peaks", tolower(run_sp))))
  load(sprintf("part1_denovo_transcripts/RData/%s_prediction.RData", run_sp))

  transfrags_peak <- unique(as.character(peak_grangs@seqnames))
  peak_num <- table(as.character(peak_grangs@seqnames))
  tmp_vec <- rep("non-m6A", nrow(ml_result$df))
  tmp_vec[rownames(ml_result$df)%in%transfrags_peak] <- "m6A"
  table(tmp_vec)
  overview_df <- ml_result$df
  overview_df$modification <- tmp_vec
  overview_df$peak_num <- table(peak_grangs@seqnames)[rownames(overview_df)]
  
  overview_df$gene <- ml_result$gffcompare[rownames(overview_df), ][, 2]
  overview_df$transcript <- ml_result$gffcompare[rownames(overview_df), ][, 3]
  
  length(which(overview_df$peak_num>=1))
  nrow(overview_df)
  
  length(intersect(which(overview_df$peak_num>=1),which(overview_df$ml_res=="busco_names")))
  
  
  
  ## coverage and identity of m6A transcripts
  plot(overview_df[overview_df$modification=="m6A", "genome_cov"],
       overview_df[overview_df$modification=="m6A", "genome_iden"], cex=0.5)
  
  ##
  fCoveragePlot(input_sp = run_sp, peak_names = "a_R3_TRINITY_DN1330_c0_g1_i1")
  ## salmon: expression values
  salmon_list <- list()
  for( rr in 1:3 ){
    salmon_list[[rr]] <- read.table(sprintf("%s/%s/03_evaluation/salmon/%s_R%s/quant.sf", 
                                            data_dir, run_sp, run_sp, rr),
                                    sep = "\t", header = T)
    salmon_list[[rr]]$repsample <- sprintf("R%s", rr)
  }
  salmon_df <- do.call('rbind', salmon_list)
  salmon_df <- reshape2::dcast(data = salmon_df, Name~repsample, value.var = "NumReads")
  rownames(salmon_df) <- salmon_df$Name
  salmon_df <- salmon_df[, -1]
  salmon_df$mean <- apply(salmon_df, 1, mean)
  salmon_df$cv <- apply(salmon_df, 1, function(x){sd(x)/mean(x)})
  salmon_df$outline <- apply(salmon_df, 1, function(x){max(x)!=0&min(x)==0})
  overview_df$expression <- salmon_df[rownames(overview_df), "mean"]
  overview_df$cv <- salmon_df[rownames(overview_df), "cv"]
  overview_df$outline <- salmon_df[rownames(overview_df), "outline"]
  ##
  filter_class <- c(overview_df$modification=="m6A"&!overview_df$outline)
  plot(overview_df[filter_class, "genome_cov"],
       overview_df[filter_class, "genome_iden"], cex=0.5)
  
  ## genome-based and denovo-based
  
  ## machine learning
  # diffbind_count <- fDiffBindCount(input_sp = run_sp)
  # hic_peak <- diffbind_count$Chr[which(diffbind_count$ratio > median(diffbind_count$ratio) &
  #                                        diffbind_count$diff > 5)]
  hic_peak <- fBUSCO(i = run_sp)
  peaks_df <- peak_grangs
  input_sp <- run_sp
}


fweakRM <- function(input_sp, peaks_df, hic_peak) {
  bed_path <- "~/miniconda3/envs/metaPlants/bin/bedtools"
  fasta_file <- sprintf("%s/%s/02_merge_assembly/tool/tools_merge.fasta", data_dir, input_sp)
  peaks_df <- as.data.frame(peaks_df)
  ## 
  pos_peak <- peaks_df[peaks_df$seqnames%in%hic_peak, 1:3]
  weakrm <- sprintf("part2_peak_calling/weakRM/%s", input_sp)
  system(sprintf("mkdir -p %s", weakrm))
  # pose
  write.table(pos_peak, file = sprintf("%s/pospre.bed", weakrm), row.names = F,
              quote = F, col.names = F, sep = "\t")
  system(command = sprintf("%s getfasta -fi %s -bed %s/pospre.bed -fo %s/pospre.fa",
                           bed_path, fasta_file, weakrm, weakrm))
  
  peakseq <- seqinr::read.fasta(sprintf("%s/pospre.fa", weakrm), 
                                seqtype = "DNA", as.string = T, forceDNAtolower = F)

  peakmotif1 <- sapply(peakseq, function(x){grepl('[AGag][AGag][Aa][Cc][ACTact]', x)})
  peakmotif2 <- sapply(peakseq, function(x){grepl('[Tt][AGag][Tt][Aa][CTct]', x)})
  table(peakmotif1&peakmotif2)
  
  filtered_hic <- names(peakseq)[peakmotif1&peakmotif2]
  filtered_hic <- lapply(filtered_hic, function(x){
    strsplit(x, split = ":|-")[[1]]
  })
  filtered_hic <- do.call('rbind', filtered_hic)
  write.table(filtered_hic, file = sprintf("%s/pos.bed", weakrm), row.names = F,
              quote = F, col.names = F, sep = "\t")
  
  
  system(command = sprintf("%s getfasta -fi %s -bed %s/pos.bed -fo %s/pos.fa",
                           bed_path, fasta_file, weakrm, weakrm))
  
  
  transfrags <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = T, forceDNAtolower = F)
  seq_length <- sapply(transfrags, nchar) 
  write.table(data.frame("name"=names(seq_length), "len"=seq_length), 
              file = sprintf("%s/transcript.size", weakrm), quote = F, row.names = F, col.names = F, sep = "\t")
  
  system(command = sprintf("%s shuffle -i %s/pos.bed -g %s/transcript.size -chrom -excl > %s/neg.bed",
                           bed_path, weakrm, weakrm, weakrm))  
  system(command = sprintf("%s getfasta -fi %s -bed %s/neg.bed -fo %s/neg.fa",
                           bed_path, fasta_file, weakrm, weakrm))
  
  
  python_dir <- "~/miniconda3/envs/denovota/bin/python"
  weakrm_path <- "/storage_server/home/zhangt/freePEA/scripts/WeakRM/Scripts"
  
  system(command = sprintf("mkdir -p %s/train", weakrm))
  system(command = sprintf("%s part2_peak_calling/split_data.py %s", python_dir, weakrm))
  
  system(command = sprintf("%s %s/token2npy.py --input_dir='%s/train/' --output_dir='%s/train/processed/'",
                           python_dir, weakrm_path, weakrm, weakrm))
  system(command = sprintf("%s %s/main.py --training=True --input_dir='%s/train/processed/' --epoch 50 --cp_dir='%s/train/processed/' --saving=True >%s_train.txt",
                           python_dir, weakrm_path, weakrm, weakrm, input_sp))
  
  # abundance and ratio
  # Test AUC:  0.8505881861300905
  # Test PRC:  0.8035540655607584
  
  # busco
  # Test AUC:  0.8633101761796254
  # Test PRC:  0.8200994168287946
  
  # Test AUC:  0.9074809160305344
  # Test PRC:  0.8518838429731139
}

final_peaks <- as.data.frame(peak_grangs)
busco_names <- fBUSCO(i = run_sp)
fasta_path <- sprintf("%s/%s/02_merge_assembly/tool/tools_merge.fasta", data_dir, run_sp)

prediction <- function(x){
  prediction_peak <- final_peaks[!final_peaks[,1]%in%busco_names, 1:3]
  
  write.table(prediction_peak, file = "part2_peak_calling/weakRM/Ath/prediction.bed", row.names = F,
              quote = F, col.names = F, sep = "\t")
  system(command = sprintf("%s -fi %s -bed %s -fo %s",
                           "~/miniconda3/envs/metaPlants/bin/bedtools getfasta", fasta_path, 
                           "./part2_peak_calling/weakRM/Ath/prediction.bed", 
                           "./part2_peak_calling/weakRM/Ath/prediction.fa"))
  system("mkdir -p part2_peak_calling/weakRM/Ath/prediction")
  system(command = "~/miniconda3/envs/denovota/bin/python part2_peak_calling/split_prediction_data.py")
  
  preddata_dir <- "part2_peak_calling/weakRM/Ath/prediction/"
  system(command = sprintf("%s %s/token2npy.py --input_dir='%s' --output_dir='%s/processed/'",
                           python_dir, weakrm_path, preddata_dir, preddata_dir))
  
  system(command = sprintf("%s %s/main.py --training=False --input_dir='%s/prediction/processed/' --cp_dir='%s/processed/'",
                           python_dir, weakrm_path, mldata_dir, mldata_dir))
  
  mldata_dir <- "part2_peak_calling/weakRM/Ath"
  prediction_score <- read.table(file = sprintf("%s/prediction/processed/prediction_score.txt", mldata_dir))
  plot(density(prediction_score$V1))
  table(prediction_score$V1>0.5)
  remove_peak_names <- prediction_peak[prediction_score$V1 < 0.5, 1]
}
