# functions

library(ggplot2)
library(dplyr)
library(snowfall)
library(randomForest)
library(pROC)
library(foreach)
library(doParallel)
library(data.tree)
library(htmlwidgets)

genome_dir <- "/distorage_server/home/zhangt/genomes"
data_dir <- "/storage_server/home/zhangt/freePEA/data"
script_dir <- "/storage_server/home/zhangt/freePEA/scripts"


fToolCount <- function(i){
  ## software support and UpSet plot ####
  ## clster file
  ff_output <- list()
  fasta_file <- sprintf("%s/%s/02_merge_assembly/tool/tools_merge.fasta", data_dir, i)
  cluster_df <- read.delim2(paste0(fasta_file, ".clstr"), sep = "\t", row.names = NULL, header = F)
  line_number <- grep("Cluster", cluster_df[,1])
  cluster_start <- line_number + 1
  cluster_end <- c(line_number[2:length(line_number)]-1, nrow(cluster_df))
  cluster_count <- sapply(1:length(cluster_start), function(x){ #56349
    t_range <- cluster_start[x]:cluster_end[x]
    t_replicates <- do.call('rbind', strsplit(cluster_df[t_range, 2], ">|_"))[,2]
    t_replicates <- unique(t_replicates)
    length(t_replicates)
  })
  cluster_soft <- sapply(1:length(cluster_start), function(x){ #56349
    t_range <- cluster_start[x]:cluster_end[x]
    t_replicates <- do.call('rbind', strsplit(cluster_df[t_range, 2], ">|_"))[,2]
    unique(t_replicates)
  })
  cluster_soft_df <- data.frame("cluster"=rep(paste0("Cluster", 1:length(cluster_start)), 
                                              cluster_count), 
                                "names"=unlist(cluster_soft))
  cluster_soft_df <- reshape2::dcast(cluster_soft_df, cluster~names)
  rownames(cluster_soft_df) <- cluster_soft_df[,1]
  cluster_soft_df <- cluster_soft_df[, -1]
  colnames(cluster_soft_df) <- c("Trinity", "rnaSPAdes", "TransABySS")
  cluster_soft_df[!is.na(cluster_soft_df)] <- TRUE
  cluster_soft_df[is.na(cluster_soft_df)] <- FALSE
  ff_output[['Upset']] <- ComplexUpset::upset(cluster_soft_df, colnames(cluster_soft_df),
                                              name='', width_ratio=0.1, min_size=10)
  
  ff_output[['Count']] <- table(cluster_count)
  ff_output[['data']] <- cluster_soft_df
  ff_output
}


fTransfragsCount <- function(i){
  ## Summary the number of transcripts
  ## species
  tool_list <- c("Trinity", "rnaSPAdes", "TransABySS")
  data_dir <- "/storage_server/home/zhangt/freePEA/data"
  tmp_species <- rep(i, each = 3)
  tmp_tool_list <- rep(tool_list, 1)
  df_transfrag <- matrix(nrow = length(tmp_tool_list), ncol = 11)
  df_transfrag[,1] <- tmp_species
  df_transfrag[,2] <- tmp_tool_list
  
  t_index1 <- df_transfrag[,2]=="Trinity"
  for(j in 1:3){
    t_trinity <- system(command = sprintf('grep ">" %s/%s/01_assembly/01_R%s_assembly/Trinity.fasta | wc -l',
                                          data_dir, i, j), intern = T)
    df_transfrag[t_index1, j+2] <- t_trinity
  }
  t_index2 <- df_transfrag[,2]=="rnaSPAdes"
  for(j in 1:3){
    t_trinity <- system(command = sprintf('grep ">" %s/%s/01_assembly/01_R%s_assembly/rnaSPAdes/transcripts.fasta | wc -l',
                                          data_dir, i, j), intern = T)
    df_transfrag[t_index2, j+2] <- t_trinity
  }
  t_index3 <- df_transfrag[,2]=="TransABySS"
  for(j in 1:3){
    t_trinity <- system(command = sprintf('grep ">" %s/%s/01_assembly/01_R%s_assembly/TransABySS/transabyss-final.fa | wc -l',
                                          data_dir, i, j), intern = T)
    df_transfrag[t_index3, j+2] <- t_trinity
  }
  for(j in tool_list){
    t_tool <- system(command = sprintf('grep ">" %s/%s/02_merge_assembly/rep/%s.fasta | wc -l',
                                       data_dir, i, j), intern = T)
    t_index4 <- df_transfrag[,2]==j
    df_transfrag[t_index4, 6] <- t_tool
  }
  for(j in tool_list){
    t_hc <- system(command = sprintf('grep ">" %s/%s/02_merge_assembly/rep/HC_%s.fasta | wc -l',
                                     data_dir, i, j), intern = T)
    t_index5 <- df_transfrag[,2]==j
    df_transfrag[t_index5, 7] <- t_hc
  }
  t_merge <- system(command = sprintf('grep ">" %s/%s/02_merge_assembly/tool/tools_merge.fasta | wc -l',
                                      data_dir, i), intern = T)
  df_transfrag[1, 8] <- t_merge
  count_tool <- fToolCount(i)
  df_transfrag[1, 9:11] <- count_tool$Count
  df_transfrag
}


fNumFlow <- function(df){
  ## results from fTransfragsCount
  # final_tree <- data.tree::Node$new(sprintf("Final: %s", df[1,8]))
  # positive_tree <- final_tree$AddChild("Positive: 55")
  # pooled_tree <- positive_tree$AddChild("Pooled: 55")
  pooled_tree <- data.tree::Node$new(sprintf("Pooled: %s", df[1,8]))
  TransABySS_tree <- pooled_tree$AddChild(sprintf("1_replicate supported:\n %s", df[3,7]))
  TransABySS_merge <- TransABySS_tree$AddChild(sprintf("1_cd-hit: %s", df[3,6]))
  TransABySS_merge$AddChild(sprintf("1_TransABySS\nR1: %s\nR2: %s\nR3: %s", df[3,3], df[3,4], df[3,5]))
  rnaSPAdes_tree <- pooled_tree$AddChild(sprintf("2_replicate supported:\n %s", df[2,7]))
  rnaSPAdes_merge <- rnaSPAdes_tree$AddChild(sprintf("2_cd-hit: %s", df[2,6]))
  rnaSPAdes_merge$AddChild(sprintf("2_rnaSPAdes\nR1: %s\nR2: %s\nR3: %s", df[2,3], df[2,4], df[2,5]))
  Trinity_tree <- pooled_tree$AddChild(sprintf("3_replicate supported:\n %s", df[1,7]))
  Trinity_merge <- Trinity_tree$AddChild(sprintf("3_cd-hit: %s", df[1,6]))
  Trinity_merge$AddChild(sprintf("3_Trinity\nR1: %s\nR2: %s\nR3: %s", df[1,3], df[1,4], df[1,5]))
  data.tree::SetGraphStyle(pooled_tree, rankdir = "BT")
  data.tree::SetEdgeStyle(pooled_tree, arrowhead = "NA", color = "grey35", penwidth = 2)
  data.tree::SetNodeStyle(pooled_tree, style = "filled, rounded", shape = "box", fillcolor = "#A984A8", 
                          fontcolor = "black", fontname = "helvetica", tooltip = data.tree::GetDefaultTooltip)
  data.tree::SetNodeStyle(TransABySS_tree, inherit = T, fillcolor = "#BDD1C5", penwidth = "1.5px")
  data.tree::SetNodeStyle(rnaSPAdes_tree, inherit = T, fillcolor = "#558ED5", penwidth = "1.5px")
  data.tree::SetNodeStyle(Trinity_tree, inherit = T, fillcolor = "#E7D86E", penwidth = "1.5px")
  pooled_tree
}

# '#558ED5',"#A984A8"
fMapGenome <- function(fasta_file, genome_file, genome_anno){
  # Mapping transfrags to genome
  tool_path <- "~/miniconda3/envs/denovota/bin"
  bam_file <- paste0(fasta_file, ".bam")
  chimera_file <- paste0(fasta_file, ".chimera.bam")
  bed_file <- paste0(fasta_file, ".bed")
  bed_intersect <- paste0(fasta_file, ".intersect.bed")
  minimap_para <- "-I 18G -ax splice:hq -uf -G 10000 --MD --secondary=no -t 32" 
  index_command <- sprintf("%s/minimap2 -I 18G -x splice -t 32 -d %s.mmi %s", 
                           tool_path, genome_file, genome_file)
  map_command <- sprintf("%s/minimap2 %s %s.mmi %s | samtools sort -O BAM - >%s",
                         tool_path, minimap_para, genome_file, fasta_file, bam_file)
  chimera_command <- sprintf("%s/samtools view -bS -f 2048 %s >%s", 
                             tool_path, bam_file, chimera_file)
  bam_index <- sprintf("%s/samtools index %s", tool_path, bam_file)
  bamtobed_command <- sprintf("%s/bedtools bamtobed -bed12 -i %s > %s", 
                              tool_path, bam_file, bed_file)
  intersect_command <- sprintf("awk '$3~/gene/' %s | %s/bedtools intersect -wo -a %s -b - > %s", 
                               genome_anno, tool_path, bed_file, bed_intersect)
  # running mapping
  if (!file.exists(bam_file)) {
    system(command = index_command)
    system(command = map_command)
    system(command = bam_index)
    system(command = chimera_command)
    system(command = bamtobed_command)
  }
  system(command = intersect_command)
  list("raw_bam" = bam_file, 
       "chimera_bam" = chimera_file, 
       "intersect_txt" = bed_intersect)

# "~/miniconda3/envs/denovota/bin/gmap -D ", dirname(index.file),
# " -d index --max-intronlength-middle=2000 --split-large-introns ",
# " -t 20 -f psl -n 0 --nofails -p 3 ", fasta.file, " > ", psl.file)
}


fGffCompare <- function(input_fasta, input_gtf, script_dir){
  ## gffcompare
  ## fasta fasta.bam gtf_file
  track_file <- sprintf("%s_gffcompare.tracking", input_fasta)
  if (!file.exists(track_file)){
    tool_path <- "~/miniconda3/envs/denovota/bin"
    gffc1 <- sprintf("%s/bedtools bamtobed -bed12 -i %s.bam > %s.bed ", 
                     tool_path, input_fasta, input_fasta)
    gffc2 <- sprintf('%s/bedToGenePred %s.bed %s.genepred', 
                     tool_path, input_fasta, input_fasta)
    gffc3 <- sprintf('%s/genePredToGtf "file" -utr %s.genepred %s.gtf', 
                     tool_path, input_fasta, input_fasta)
    gffcompare <- sprintf("%s/gffcompare/gffcompare", script_dir)
    gffc4 <- sprintf("%s -r %s -o %s_gffcompare %s.gtf", 
                     gffcompare, input_gtf, input_fasta, input_fasta)
    system(command = gffc1)
    system(command = gffc2)
    system(command = gffc3)
    system(command = gffc4)
  }
  
  compare_gff_df1 <- read.table(sprintf("%s_gffcompare.%s.gtf.tmap",input_fasta, basename(input_fasta)),
                                sep = "\t", header =T)
  compare_gff_df <- read.table(track_file, sep = "\t")
  compare_gff_df <- compare_gff_df[compare_gff_df$V3!="-", ]
  compare_gff_mat <- as.data.frame(t(apply(compare_gff_df, 1, function(x){
    f_gene_names <- strsplit(x[3], "\\|")[[1]]
    f_transfrags_names <- strsplit(x[5], "\\|")[[1]][1]
    f_transfrags_names <- gsub("q1:", "", f_transfrags_names)
    f_type <- x[4]
    return(c(f_transfrags_names, f_gene_names[1], f_gene_names[2], f_type))
  })))
  rownames(compare_gff_mat) <- compare_gff_mat[,1]
  compare_gff_mat
}


fBinaryTree <- function(input_df, true_contigs, remove_list, intersect_file){
  ## map information
  # html_file <- paste0("part1_denovo_transcripts/results_plot/Part_1/plot.html")
  # htmlwidgets::saveWidget(plot(mapping_info$tree), file = html_file)
  input_df <- input_df[intersect(true_contigs, rownames(input_df)), ]
  filter_index <- input_df$coverage>=80&input_df$identity>=80
  table(filter_index)
  remain_df <- input_df[filter_index,]; f_map <- unique(rownames(remain_df))
  filter_df <- input_df[!filter_index,]
  # 
  f_intersect <- read.table(intersect_file, sep="\t")
  f_i_ncol <- ncol(f_intersect)
  f_intersect$coverage1 <- round(f_intersect[, f_i_ncol]/(f_intersect[,3]-f_intersect[,2])*100, digits = 2)
  f_intersect$coverage2 <- round(f_intersect[, f_i_ncol]/(f_intersect[,f_i_ncol-5]-f_intersect[,f_i_ncol-6]+1)*100, digits = 2)
  f_intersect_filter <- f_intersect %>% filter(coverage1 >= 50 | coverage2 >= 50 )
  # 
  f_genic <- unique(f_intersect_filter$V4)
  # decision tree
  f_HC_df <- subset(remain_df, coverage>=98&identity>=98); f_HC <- unique(rownames(f_HC_df))
  f_LC_df <- subset(remain_df, coverage<98|identity<98); f_LC <- unique(rownames(f_LC_df))
  f_LC1 <- unique(rownames(subset(remain_df, coverage>=80&coverage<90|identity>=80&identity<90)))
  f_LC2 <- unique(rownames(subset(remain_df, coverage>=90&coverage<95|identity>=90&identity<95)))
  f_LC3 <- unique(rownames(subset(remain_df, coverage>=95&coverage<98|identity>=95&identity<98)))
  f_perfect <- unique(rownames(subset(f_HC_df, coverage==100&identity==100)))
  # 
  f_unperfect_df <- subset(f_HC_df, coverage!=100|identity!=100)
  f_unperfect <- unique(rownames(f_unperfect_df))
  f_incom_iden <- unique(rownames(subset(f_unperfect_df, coverage==100&identity<100)))
  f_incom_cov <- unique(rownames(subset(f_unperfect_df, coverage<100&identity==100)))
  f_incom_both <- unique(rownames(subset(f_unperfect_df, coverage<100&identity<100)))
  # The genome distribution of perfect alignment
  f_com_genic <- unique(intersect(f_perfect, f_genic))
  f_com_inter <- unique(setdiff(f_perfect, f_genic))
  # The genome distribution of unperfect alignment
  f_incom_igenic <- unique(intersect(f_incom_iden, f_genic))
  f_incom_iinter <- unique(setdiff(f_incom_iden, f_genic))
  f_incom_cgenic <- unique(intersect(f_incom_cov, f_genic))
  f_incom_cinter <- unique(setdiff(f_incom_cov, f_genic))
  f_incom_bgenic <- unique(intersect(f_incom_both, f_genic))
  f_incom_binter <- unique(setdiff(f_incom_both, f_genic))
  # Output gene class
  f_class <- rep("raw", length(true_contigs)); names(f_class) <- true_contigs
  f_class[f_map] <- "Remain"; f_class[setdiff(true_contigs, f_map)] <- "Remove"
  f_class[f_HC] <- "HC"; f_class[f_LC] <- "LC"
  f_class[f_perfect] <- "Perfect"; f_class[f_unperfect] <- "Unperfect"
  f_class[f_com_genic] <- "P_Genic"; f_class[f_com_inter] <- "P_Intergenic"
  f_class[remove_list$Chimeric1] <- "Chimera1"
  f_class[remove_list$Chimeric2] <- "Chimera2"
  f_class[remove_list$Unmapped] <- "Unmapped"
  # Construct tree --------------
  f_tree <- data.tree::Node$new(paste0("The assembled\ntranscript fragments\n", length(true_contigs)))
  f_t_map <- f_tree$AddChild(paste0("Mapped\n", nrow(input_df)))
  f_t_chimera <- f_tree$AddChild(paste0("Chimeric\n", length(c(remove_list$Chimeric1, remove_list$Chimeric2))))
  f_t_chimera$AddChild(paste0("Same\nchromosome\n", length(remove_list$Chimeric1)))
  f_t_chimera$AddChild(paste0("Different\nchromosomes\n", length(remove_list$Chimeric2)))
  f_t_unmap <- f_tree$AddChild(paste0("Unmapped\n", length(remove_list$Unmapped)))
  f_t_80map <- f_t_map$AddChild(paste0("Remain\nCov>=0.8&\nIden>=0.8\n", length(f_map)))
  f_t_remove <- f_t_map$AddChild(paste0("Remove\n", length(unique(rownames(filter_df)))))
  f_t_remove$AddChild(paste0("Genic\n", length(unique(intersect(rownames(filter_df), f_genic)))))
  f_t_remove$AddChild(paste0("Intergenic\n", length(unique(setdiff(rownames(filter_df), f_genic)))))
  f_t_LC <- f_t_80map$AddChild(paste0("Low-confidence\n", length(f_LC)))
  f_t_HC <- f_t_80map$AddChild(paste0("High-confidence\nCov>=0.98&\nIden>=0.98\n", 
                                                 length(f_HC)))
  f_t_LC1 <- f_t_LC$AddChild(paste0("Cov∈[0.8,0.9)&\nIden∈[0.8.9)\n", 
                                    length(f_LC1)))
  f_t_LC1$AddChild(paste0("Genic\n", length(unique(intersect(f_LC1, f_genic)))))
  f_t_LC1$AddChild(paste0("Intergenic\n", length(unique(setdiff(f_LC1, f_genic)))))
  f_t_LC2 <- f_t_LC$AddChild(paste0("Cov∈[0.9,0.95)&\nIden∈[0.9,0.95)\n", 
                                          length(f_LC2)))
  f_t_LC2$AddChild(paste0("Genic\n", length(unique(intersect(f_LC2, f_genic)))))
  f_t_LC2$AddChild(paste0("Intergenic\n", length(unique(setdiff(f_LC2, f_genic)))))
  f_t_LC3 <- f_t_LC$AddChild(paste0("Cov∈[0.95,0.98)&\nIden∈[0.95,0.98)\\n", 
                                          length(f_LC3)))
  f_t_LC3$AddChild(paste0("Genic\n", length(unique(intersect(f_LC3, f_genic)))))
  f_t_LC3$AddChild(paste0("Intergenic\n", length(unique(setdiff(f_LC3, f_genic)))))
  f_t_com <- f_t_HC$AddChild(paste0("Perfect\nCov=1 &den=1\\n", length(f_perfect)))
  f_t_com$AddChild(paste0("Genic\n", length(f_com_genic)))
  f_t_com$AddChild(paste0("Intergenic\n", length(f_com_inter)))
  # unperfect
  f_t_incom <- f_t_HC$AddChild(paste0("Unperfect\n", length(f_unperfect)))
  f_t_incom_iden <- f_t_incom$AddChild(paste0("Cov=1&\nIden∈[0.98,1)\n", length(f_incom_iden)))
  f_t_incom_iden$AddChild(paste0("Genic\n", length(f_incom_igenic)))
  f_t_incom_iden$AddChild(paste0("Intergenic\n", length(f_incom_iinter)))
  f_t_incom_cov <- f_t_incom$AddChild(paste0("Cov∈[0.98,1)&\nIden=1\n", length(f_incom_cov)))
  f_t_incom_cov$AddChild(paste0("Genic\n", length(f_incom_cgenic)))
  f_t_incom_cov$AddChild(paste0("Intergenic\n", length(f_incom_cinter)))
  f_t_incom_both <- f_t_incom$AddChild(paste0("Cov∈[0.98,1)&\nIden∈[0.98,1)\n", length(f_incom_both)))
  f_t_incom_both$AddChild(paste0("Genic\n", length(f_incom_bgenic)))
  f_t_incom_both$AddChild(paste0("Intergenic\n", length(f_incom_binter)))
  
  # Tree Style
  data.tree::SetGraphStyle(f_tree, rankdir = "TB")
  data.tree::SetEdgeStyle(f_tree, arrowhead = "vee", color = "grey35", penwidth = 2)
  data.tree::SetNodeStyle(f_tree, style = "filled, rounded", shape = "box", fillcolor = "LightBlue", 
               fontcolor = "black", fontname = "helvetica", tooltip = data.tree::GetDefaultTooltip)
  data.tree::SetNodeStyle(f_t_80map, inherit = T, fillcolor = "GreenYellow", penwidth = "1.5px")
  
  return(list("tree" = f_tree, "class"= f_class))
}


fPieRatio <- function(f_itable, main){
  pie_data <- as.data.frame(f_itable, stringsAsFactors = F)
  pie_data <- pie_data[order(pie_data[,2], decreasing = T), ]
  colnames(pie_data)[1] <- "Var1"
  rownames(pie_data) <- pie_data[,1]
  tmp_level <- intersect(c("known", "long", "long_m", "long_k", "short", "short_j/n", "short_c/e", "unmapped", 
                           "sc_chimera", "dc_chimera", "other"), unique(pie_data$Var1))
  pie_data <- pie_data[tmp_level, ]
  pie_data$Var1 <- factor(pie_data$Var1, levels = tmp_level)
  pie_data$per <- round(pie_data$Freq/sum(pie_data$Freq), digits = 4)
  pie_data$label <- paste0(pie_data$per*100, "%")

  print(ggplot(data = pie_data,
         mapping = aes(x = "Count", y = per, fill = Var1)) +
    # scale_fill_manual(values = c('#BDD1C5','#4D6372','#E7D86E','#558ED5',"#A984A8"), guide = FALSE)
    geom_bar(stat = 'identity', position = 'fill', width = 0.5) +
    labs(x = '', y = '', title = '') +
    coord_polar("y", direction = -1) +
    theme_void(base_size = 12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    scale_fill_manual(values = rev(soilpalettes::soil_palette("redox2", 
                                                              length(unique(pie_data$Var1))))) +
    ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) +
    geom_text(aes(x=1.3, y = sum(per) - cumsum(per) + per/2, label=label)) + 
    guides(fill = guide_legend(title = NULL))
  )
}


fSeq2Genome <- function(i, fasta_file){
  # sequence to genome
  ## fMapGenome
  genome_dir <- "/distorage_server/home/zhangt/genomes"
  data_dir <- "/storage_server/home/zhangt/freePEA/data"
  script_dir <- "/storage_server/home/zhangt/freePEA/scripts"
  
  genome_file <- sprintf("%s/%s/Genome/%s.fa", genome_dir, i, i)
  genome_anno <- sprintf("%s/%s/Annotation/%s.exons.gff3", genome_dir, i, i)
  gtf_file <- sprintf("%s/%s/Annotation/%s.exons.gtf", genome_dir, i, i)
  # fasta_file <- sprintf("%s/tools_merge.fasta", merge_dir)
  # bufasta_file <- sprintf("%s/bu_tools_merge.fasta", merge_dir)
  # system(command = sprintf('cp %s %s && sed -r "s/A+$//" %s > %s', 
  #                          fasta_file, bufasta_file, bufasta_file, fasta_file))
  transfrags <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = T, forceDNAtolower = F)
  # map contigs to genome using minimap2
  map_file <- fMapGenome(fasta_file = fasta_file, genome_file = genome_file, 
                         genome_anno = genome_anno)
  # chimera from map_file
  chimera_bam_file <- Rsamtools::BamFile(map_file$chimera_bam)
  chimera_aln <- Rsamtools::scanBam(chimera_bam_file)[[1]]
  chimera_transfrags <- unique(chimera_aln$qname)
  bam_file <- Rsamtools::BamFile(map_file$raw_bam)
  total_aln <- Rsamtools::scanBam(bam_file)[[1]]
  # unmapped transfrags
  unmapped_transfrags <- total_aln$qname[which(is.na(total_aln$cigar))]
  # Two chimera types: same chromosome and different chromosomes
  ct_df <- data.frame("name"=total_aln$qname[total_aln$qname%in%chimera_transfrags],
                      "chr"=total_aln$rname[total_aln$qname%in%chimera_transfrags])
  ct_df <- ct_df %>% group_by(name) %>% summarise(cn_num=length(unique(chr)))
  # remove unmapped and chimera
  mapped_transfrags <- total_aln$qname[!total_aln$qname%in%c(chimera_transfrags, 
                                                             unmapped_transfrags)]
  
  # # Decision Tree
  genome_cov_ide <- sprintf("%s_bam_iden_cov.csv", fasta_file)
  genome_bam_file <- sprintf("%s.bam", fasta_file)
  python_dir <- "~/miniconda3/envs/denovota/bin/python"
  # if (!file.exists(genome_cov_ide)){
    system(command = sprintf("%s %s/bam2idencov.py %s %s",
                             python_dir, script_dir,
                             genome_bam_file, genome_cov_ide))
  # }
  # python output of identity and coverage of contigs map to genome
  cigar_df <- read.table(genome_cov_ide, sep = ",", header = T)
  
  mapped_cigar_df <- cigar_df[cigar_df[,1]%in%mapped_transfrags, ]
  rownames(mapped_cigar_df) <- mapped_cigar_df[,1]
  
  # 
  remove_list <- list("Chimeric1"=ct_df$name[ct_df$cn_num==1], 
                      "Chimeric2"=ct_df$name[ct_df$cn_num>1],
                      "Unmapped"=unique(unmapped_transfrags))
  
  cat("Total: ", length(names(transfrags)),
      "Unmapped: ", length(unique(unmapped_transfrags)),
      "Chimeric1: ", length(unique(ct_df$name[ct_df$cn_num==1])),
      "Chimeric2: ", length(unique(ct_df$name[ct_df$cn_num>1])),
      "\n")
  # inter_file <- read.table(map_file$intersect_txt, sep = "\t")
  # hybrid_transfrags <- names(table(inter_file$V4))[table(inter_file$V4)>=2]
  # contig.tree <- fBinaryTree(input_df = mapped_cigar_df, 
  #                            true_contigs = names(transfrags),
  #                            remove_list = remove_list,
  #                            intersect_file = map_file$intersect_txt)
  # contig.tree
  compare_out <- fGffCompare(input_fasta = fasta_file, input_gtf = gtf_file, 
                             script_dir =script_dir)
  ## "=": known, "k", "m": long, "c","n","j","e":short  "other"
  # Transcript classification
  transcript_out <- rep("other", length(transfrags))
  names(transcript_out) <- names(transfrags)
  transcript_out[unique(unmapped_transfrags)] <- "unmapped"
  transcript_out[unique(ct_df$name[ct_df$cn_num==1])] <- "sc_chimera"
  transcript_out[unique(ct_df$name[ct_df$cn_num>1])] <- "dc_chimera"
  
  transcript_out[names(transfrags)%in%compare_out$V1[compare_out$V4%in%c("c","e")]] <- "short_c/e"
  transcript_out[names(transfrags)%in%compare_out$V1[compare_out$V4%in%c("n","j")]] <- "short_j/n"
  transcript_out[names(transfrags)%in%compare_out$V1[compare_out$V4%in%c("k")]] <- "long_k"
  transcript_out[names(transfrags)%in%compare_out$V1[compare_out$V4%in%c("m")]] <- "long_m"
  transcript_out[names(transfrags)%in%compare_out$V1[compare_out$V4%in%c("=")]] <- "known"
  
  # transcript_out[hybrid_transfrags] <- "hybrid_transfrags"
  output_list <- list()
  output_list[['gff_table']] <- compare_out
  output_list[['out']] <- transcript_out
  output_list[['cov_iden']] <- mapped_cigar_df
  output_list
}


fRadarPlot <- function(df, f_col){
  library(ggradar)
  print(ggradar(df, 
                axis.label.size = 5,
                group.line.width = 0.8,group.point.size = 2,
                grid.label.size = 4,
                legend.text.size =8, 
                values.radar = c("0%", "35%", "70%"),
                grid.mid = 0.35,
                grid.max = 0.7,
                legend.position = "right",
                background.circle.colour = 'white',
                gridline.min.colour = 'grey80',
                gridline.mid.colour = 'grey80',
                gridline.max.colour = 'grey80',
                gridline.min.linetype = 'solid',
                gridline.mid.linetype = 'solid',
                gridline.max.linetype = 'solid',
                axis.line.colour = "grey80",
                plot.extent.x.sf = 1.2) + 
          scale_color_manual(values = alpha(f_col, 0.7)))
}


plotROC <- function(cvRes) {
  
  cvListPredictions <- list()
  cvListLabels <- list()
  AUCVec <- rep(0, length(cvRes) )
  for( i in 1:length(cvRes) ) {
    curCV <- cvRes[[i]]
    cvListPredictions[[i]] <- c( curCV$positives.test.score, curCV$negatives.test.score )
    cvListLabels[[i]] <- c( rep(1, length(curCV$positives.test.score)), rep(0, length(curCV$negatives.test.score) ) )
    AUCVec[i] <- curCV$test.AUC
  }
  mAUC <- format( mean(AUCVec), digits= 3)
  
  #if( !require(ROCR) ) {
  #   install.packages("ROCR")
  #   library(ROCR)
  #}
  pred <- ROCR::prediction(cvListPredictions, cvListLabels)
  perf <- ROCR::performance(pred,"tpr","fpr")
  
  
  par(mar=c(5,6,4,2))   
  plot(perf, col= "gray", lty=3, main = paste( "AUC = ", mAUC, sep = ""), cex.lab = 2.5, cex.axis = 2, cex.main = 3, mgp = c(4,1.8,0) )
  plot(perf, col = "black",  lwd= 3, avg="vertical", spread.estimate="none", add=TRUE)  
  
}


evalPrediction <- function(threshold,
                           posScore = NULL, negScore = NULL,
                           beta = 1,
                           TP, TN, FP, FN){
  
  if(is.null(posScore) & is.null(negScore)){
    Sn <- TP/(TP+FN)
    Sp <- TN/(TN+FP)
    Pr <- TP/(TP+FP)
    Acc <- (TP+TN)/(TP+TN+FP+FN)
    
    Fscore <- ((1+beta^2)*Pr*Sn)/(beta^2*Pr+Sn)
    a <- as.numeric(TP+FP)
    b <- as.numeric(TP+FN)
    c <- as.numeric(TN+FP)
    d <- as.numeric(TN+FN)
    MCC <- (TP*TN-FP*FN)/sqrt( a*b*c*d )
    res <- c(round(c(Sn, Sp, Pr, Acc, Fscore, MCC), digits = 2), TP, FP, TN, FN)
    names(res) <- c("Sn", "Sp", "Pr", "Acc", "F-score", "MCC", "TP", "FP", "TN", "FN")
  }else{
    TP <- as.numeric(length(which(posScore >= threshold)))
    FN <- as.numeric(length(which(posScore < threshold)))
    TN <- as.numeric(length(which(negScore < threshold)))
    FP <- as.numeric(length(which(negScore >= threshold)))
    
    Sn <- TP/(TP+FN)
    Sp <- TN/(TN+FP)
    Pr <- TP/(TP+FP)
    Acc <- (TP+TN)/(TP+TN+FP+FN)
    
    Fscore <- ((1+beta^2)*Pr*Sn)/(beta^2*Pr+Sn)
    MCC <- (TP*TN-FP*FN)/sqrt(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
    AUC <- pROC::roc( c(rep(1, length(posScore)), rep(0, length(negScore))),
                      c(posScore, negScore) )$auc[1]
    res <- c(Sn, Sp, Pr, Acc, Fscore, MCC, AUC)
    names(res) <- c("Sn", "Sp", "Pr", "Acc", "F-score", "MCC", "AUC")
  }
  res
}


.winCount <- function(x, y, z){
  data.frame(
    "Position" = subset(z, !duplicated(x)),
    "window" = subset(x, !duplicated(x)),
    "curWave" = as.numeric(by(y, factor(x, levels = unique(x)), mean))
  )
}


.intervalCov <- function(dzFile){
  library(dplyr)
  baseCov <- data.table::fread(file = dzFile, sep = "\t", header = F, stringsAsFactors = F)
  # baseCov <- data.frame(baseCov, stringsAsFactors = F)
  winDow <- rep(NA, nrow(baseCov))
  baseCov <- cbind(baseCov, winDow)
  colnames(baseCov) <- c("Chr", "Position", "readsNumber", "window")
  baseCov$Position <- baseCov$Position + 1
  baseCov$window <- ceiling(baseCov$Position/25)
  
  baseCovOutput <- baseCov %>% 
    group_by(Chr) %>% 
    do(.winCount(x=.$window, y=.$readsNumber, z=.$Position))
  
  baseCovOutput <- split(baseCovOutput, factor(baseCovOutput$Chr, levels = unique(baseCovOutput$Chr)))
  baseCovOutput
}


fBUSCO <- function(i){
  # busco result ####
  busco_class <- list("Ath"="brassicales", "Gar"="eudicots", "Ghi"="eudicots", 
                      "Pvu"="fabales", "Gma"="fabales", "Sbi"="poales", "Zma"="poales", 
                      "Ata"="poales", "Tdi"="poales", "Tae"="poales", "Osa"="poales",
                      "Ppa"="embryophyta")
  data_dir <- "/storage_server/home/zhangt/freePEA/data"
  busco_ctg <- read.delim2(sprintf("%s/%s/03_evaluation/BUSCO/busco/run_%s_odb10/full_table.tsv", 
                                   data_dir, i, busco_class[[i]]), sep = "\t", header = F)
  busco_ctg <- busco_ctg[busco_ctg[,2]%in%c("Complete", "Duplicated"), ]
  busco_ctg_names <- unique(apply(busco_ctg, 1, function(x){strsplit(x[3], split = ":")[[1]][1]}))
  busco_ctg_names
}


fTransforDF <- function(x){
  f_data <- x^2+1
  f_data <- sqrt(f_data)
  f_data <- x + f_data
  log(f_data)
}


fROCplot <- function(pred){
  f_pred <- lapply(pred, function(x){
    perf <- ROCR::performance(x,"tpr","fpr")
    dat= data.frame(fpr_pu = perf@x.values[[1]],
                    tpr_pu = perf@y.values[[1]])
    dat
  })
  f_pred_mat <- do.call('rbind', f_pred)
  f_pred_mat$type <- rep(names(pred), sapply(pred, function(x){length(ROCR::performance(x,"tpr","fpr")@x.values[[1]])}))
  
  f_auc <- sapply(pred, function(x){
    auc_ROCR <- ROCR::performance(x, measure = "auc")
    auc_ROCR@y.values[[1]]
  })
  cat(f_auc, "\n")
  ggplot(f_pred_mat, aes(.data$fpr_pu, .data$tpr_pu, color=type))+
    geom_path(size=1, linetype="solid")+
    # geom_path(aes(.data$fpr,.data$tpr),linetype="solid",color="black",size=0.5)+
    geom_abline(slope = 1, linetype="dashed",color="darkblue",size=0.8)+
    theme_bw(base_size = 14)+
    xlab("false-positive rate")+
    ylab("true-positive rate")+
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))+
    annotate("text",label=paste("AUC: ", paste(paste(names(pred), ": ", round(f_auc,digits = 2)), collapse = "\n"), sep="\n"), x=0.75, y=0.2, size=3)+
    coord_fixed() + 
    BaseTheme() +
    scale_color_manual(values = soilpalettes::soil_palette("podzol", length(pred))) +
    guides(fill = guide_legend(title = ''))
}


fMLpipe <- function(i){
  ## types
  fasta_path <- sprintf("%s/%s/02_merge_assembly/tool/tools_merge.fasta", data_dir, i)
  f_transfrags <- seqinr::read.fasta(fasta_path, seqtype = "DNA", as.string = T,
                                     forceDNAtolower = F)
  seq_type <- fSeq2Genome(i = i, fasta_file = fasta_path)
  true_positive <- names(seq_type$out)[seq_type$out%in%c("known", "long_k", "long_m")]
  true_negative <- setdiff(names(f_transfrags), true_positive)
  
  # combind_vec <- vector()
  # for(rep in 1:3){
  #   for(type in c("input", "ip")){
  #     combind_vec <- c(combind_vec,
  #                      sprintf("%s/%s/03_evaluation/coverage/%s_R%s.%s.txt",
  #                              data_dir, i, i, rep, type))
  #   }
  # }
  # library(foreach)
  # library(doParallel)
  # cl <- makeCluster(6)
  # registerDoParallel(cl)
  # return_list <- foreach(input = combind_vec, .combine = list, .multicombine = TRUE)  %dopar%   {
  #   readcov <- .intervalCov(dzFile = input)
  #   f_name <- gsub(".txt", "", basename(input))
  #   feature_list <- list()
  #   feature_list[["meancov"]] <- sapply(readcov, function(x){mean(x$curWave)})[names(f_transfrags)]
  #   feature_list[["sdcov"]] <- sapply(readcov, function(x){sd(x$curWave)})[names(f_transfrags)]
  #   feature_list[["cvcov"]] <- sapply(readcov, function(x){sd(x$curWave)/mean(x$curWave)*100})[names(f_transfrags)]
  #   feature_list[["zero"]] <- sapply(readcov, function(x){sum(x$curWave<5)/length(x$curWave)*100})[names(f_transfrags)]
  #   names(feature_list[["meancov"]]) <- names(feature_list[["sdcov"]]) <- 
  #     names(feature_list[["cvcov"]]) <- names(feature_list[["zero"]]) <- names(f_transfrags)
  #   # 8 bins
  #   contigs.covbin <- lapply(readcov, function(zz){
  #     zz.len <- ceiling(nchar(f_transfrags[as.character(zz[1,1])])/25)
  #     zz.tmp.vec <- rep(0, zz.len)
  #     zz.tmp.vec[zz$window] <- zz$curWave
  #     tmpvalue <- aggregate(zz.tmp.vec, by=list(cut(1:length(zz.tmp.vec),breaks = 8)), mean)$x
  #     round(tmpvalue, digits = 2)
  #   })
  #   contigs.covbin.mat <- do.call('rbind', contigs.covbin)
  #   contigs.covbin.df <- as.data.frame(contigs.covbin.mat)
  #   contigs.covbin.df <- contigs.covbin.df[names(f_transfrags), ]
  #   rownames(contigs.covbin.df) <- names(f_transfrags)
  #   colnames(contigs.covbin.df) <- as.character(1:8)
  #   feature_list[["bincov"]] <- t(contigs.covbin.df)
  #   # merge features
  #   feature.mat <- t(do.call('rbind', feature_list))
  #   feature.df <- as.data.frame(feature.mat)
  #   feature.df[is.na(feature.df)] <- 0
  #   colnames(feature.df) <- paste0(f_name, "_", colnames(feature.df))
  #   # feature.df$ID <- rownames(feature.df)
  #   # feature.df$type <- strsplit(f_name, split = "_|\\.")[[1]][3]
  #   feature.df
  # }
  # stopImplicitCluster()
  # return_df <- do.call('cbind', return_list)
  # # read_cov_df <- return_df %>% group_by(ID, type) %>% summarise_all("mean")
  # 
  # salmon_data_list <- list()
  # for(rep in 1:3){
  #   feature_data <- read.table(sprintf("%s/%s/03_evaluation/salmon/%s_R%s/quant.sf", 
  #                                      data_dir, i, i, rep),header = T, row.names = 1)
  #   salmon_data_list[[rep]] <- feature_data[names(f_transfrags), c(3,4)]
  #   rownames(salmon_data_list[[rep]]) <- names(f_transfrags)
  #   colnames(salmon_data_list[[rep]]) <- paste0("R", rep, "_", colnames(salmon_data_list[[rep]]))
  # }
  # salmon_data_df <- do.call('cbind', salmon_data_list)
  
  # feature_other <- read.table(sprintf("%s/%s/03_evaluation/TransRate/R1/tools_merge/contigs.csv", data_dir, i),
  #                             header = T, sep = ",", row.names = 1)[names(f_transfrags),1:3]
  # rownames(feature_other) <- names(f_transfrags)
  # # 72 6 3
  # feature.merge.df <- cbind(return_df, salmon_data_df, feature_other) #transrate_data_df
  # 
  # system(command = "mkdir -p part1_denovo_transcripts/RData")
  # save(feature.merge.df, file = sprintf("part1_denovo_transcripts/RData/%s_feauteMat.RData", i))
  # ## save result
  # feature_df_list[[i]] <- feature.merge.df
  
  load(sprintf("part1_denovo_transcripts/RData/%s_feauteMat.RData", i))
  ## busco
  busco_names <- fBUSCO(i = i)

  # # fFeaturePlot(fea.mat = fTransforDF(feature.merge.df), pos.tmp = true_positive,
  # #              neg.tmp = true_negative)
  # # ggsave(filename = "tmp.pdf", width = 5, height = 40)
  
  system(command = "mkdir -p part1_denovo_transcripts/PU_dir")
  pu_python_input <- sprintf("part1_denovo_transcripts/PU_dir/%s_PU_input.txt", i)
  pu_python_output <- sprintf("part1_denovo_transcripts/PU_dir/%s_PU_output.txt", i)
  
  if(!file.exists(pu_python_output)){
    pu_input_df <- feature.merge.df
    all(busco_names%in%rownames(pu_input_df))
    rm(feature.merge.df)
    pu_input_df[is.na(pu_input_df)] <- 0
    pu_input_df <- as.data.frame(scale(fTransforDF(pu_input_df)))
    tmp.type <- rep(0, nrow(pu_input_df))
    tmp.type[rownames(pu_input_df) %in% busco_names] <- 1
    # table(tmp.type)
    # true.tmp <- rep(0, nrow(pu_input_df))
    # true.tmp[rownames(pu_input_df) %in% true_positive] <- 1
    # table(true.tmp)
    pu_input_df$authentic <- tmp.type
    # pu_input_df$final <- true.tmp
    # cat(sum(pu_input_df[,"authentic"]==1),
    #     sum(pu_input_df[,"authentic"]!=1&pu_input_df[,"final"]==1),
    #     sum(pu_input_df[,"authentic"]!=1&pu_input_df[,"final"]!=1), "\n")

    # reod.cal <- c(which(pu_input_df[,"authentic"]==1),
    #               which(pu_input_df[,"authentic"]!=1&pu_input_df[,"final"]==1),
    #               which(pu_input_df[,"authentic"]!=1&pu_input_df[,"final"]!=1))
    # pu_input_df <- pu_input_df[reod.cal, 1:(ncol(pu_input_df)-1)]

    write.table(pu_input_df, file = pu_python_input, quote = F, sep = "\t")

    PU_command <- sprintf("~/miniconda3/envs/metaPlants/bin/python %s/PUlearning.py %s %s",
                          script_dir, pu_python_input, pu_python_output)

    system(command = PU_command)
  }
  
  ## recall
  pu_result <- read.table(pu_python_output, sep = ",", row.names = 1)
  # plot(density(pu_result[,1]))
  
  Vthresold <- 0.1
  pu_remove <- unique(rownames(pu_result)[pu_result[,1] <= Vthresold])
  pu_remain <- setdiff(rownames(pu_result), pu_remove)
  
  # all(pu_remove%in%names(mapping_info)); all(pu_remain%in%names(mapping_info))
  # fPieRatio(f_itable = table(mapping_info[pu_remain]), main = "Positive")
  # fPieRatio(f_itable = table(mapping_info[pu_remove]), main = "Negative")
  # re_prediction <- fPrediction(positive = busco_names, negative = pu_remove, 
  #             unlabeled = pu_remain, featuresmat = feature.merge.df)
  # pu_remove <- c(pu_remove, setdiff(pu_remain, re_prediction))
  # pu_remain <- re_prediction
  
  TP <- length(intersect(true_positive, pu_remain))
  FN <- length(intersect(true_positive, pu_remove))
  FP <- length(intersect(true_negative, pu_remain))
  TN <- length(intersect(true_negative, pu_remove))

  # evalPrediction(TP=TP, TN=TN, FP=FP, FN=FN)
  
  ## save result
  ml_result <- evalPrediction(TP=TP, TN=TN, FP=FP, FN=FN)
  
  pu_type1 <- intersect(rownames(pu_result), true_positive)
  pu_type2 <- setdiff(rownames(pu_result), true_positive)
  cvListPredictions <- c(pu_result[pu_type1, 1], pu_result[pu_type2, 1])
  cvListLabels <- c( rep(1, length(pu_type1)), rep(0, length(pu_type2) ) )
  PU_pred <- ROCR::prediction(cvListPredictions, cvListLabels)
  auc_ROCR <- ROCR::performance(PU_pred, measure = "auc")
  ml_result <- c(ml_result, as.character(round(auc_ROCR@y.values[[1]], digits = 3)))
  names(ml_result)[length(ml_result)] <- "AUC"
  
  # TransRate
  transrate_data_list <- list()
  for(rep in 1:3){
    feature_data <- read.table(sprintf("%s/%s/03_evaluation/TransRate/R%s/tools_merge/contigs.csv",
                                       data_dir, i, rep),header = T, sep = ",", row.names = 1)
    feature_col <- c("score", "sCnuc", "sCcov", "sCord", "sCseg")
    transrate_data_list[[rep]] <- feature_data[names(f_transfrags), feature_col]
    rownames(transrate_data_list[[rep]]) <- names(f_transfrags)
    colnames(transrate_data_list[[rep]]) <- paste0("R", rep, "_", colnames(transrate_data_list[[rep]]))
  }
  transrate_data_df <- do.call('cbind', transrate_data_list)
  
  tr_type1 <- intersect(names(f_transfrags), true_positive)
  tr_type2 <- setdiff(names(f_transfrags), true_positive)
  tr_cvListPredictions <- c(transrate_data_df[tr_type1, "R1_score"], transrate_data_df[tr_type2, "R1_score"])
  tr_cvListLabels <- c( rep(1, length(tr_type1)), rep(0, length(tr_type2) ) )
  tr_pred <- ROCR::prediction(tr_cvListPredictions, tr_cvListLabels)

  # RCOC plot
  fROCplot(pred = list("TransRate"=tr_pred, "PU-learning"=PU_pred))
  ggsave(sprintf("part1_denovo_transcripts/results_plot/%s_PUlearning.pdf", i), 
         width = 5, height = 5)

  tmp_vec <- rep("ml_remove", length(names(f_transfrags)))
  tmp_vec[names(f_transfrags)%in%pu_remain] <- "ml_remain"
  tmp_vec[names(f_transfrags)%in%busco_names] <- "busco_names"
  
  out_df <- data.frame("transfrags"=names(f_transfrags), 
                       "length"=sapply(f_transfrags, nchar),
                       "genome_type" = seq_type$out[names(f_transfrags)],
                       "genome_cov" = seq_type$cov_iden[names(f_transfrags), "coverage"],
                       "genome_iden"= seq_type$cov_iden[names(f_transfrags), "identity"],
                       "ml_res" = tmp_vec,
                       stringsAsFactors = F)
  
  return(list("ml"=ml_result, "df"=out_df, 
              "gffcompare"=seq_type$gff_table, "cov_iden"=seq_type$cov_iden))
}


fCrossSpeciesPred <- function(x) {
  # species 1
  sp1_fasta_path <- sprintf("%s/%s/02_merge_assembly/tool/tools_merge.fasta", data_dir, x[1])
  sp1_transfrags <- seqinr::read.fasta(sp1_fasta_path, seqtype = "DNA", as.string = T,
                                     forceDNAtolower = F)
  sp1_seq_type <- fSeq2Genome(i = x[1], fasta_file = sp1_fasta_path)
  sp1_true_positive <- names(sp1_seq_type)[sp1_seq_type%in%c("known", "long")]
  sp1_true_negative <- setdiff(names(sp1_transfrags), sp1_true_positive)
  load(sprintf("part1_denovo_transcripts/RData/%s_feauteMat.RData", x[1]))
  sp1_feature_mat <- feature.merge.df
  rm(feature.merge.df)
  sp1_model.denovo <- cross_validation(method="randomForest",
                                   featureMat = as.data.frame(scale(fTransforDF(sp1_feature_mat))),
                                   positives = sp1_true_positive, 
                                   negatives = sp1_true_negative,
                                   cross = 5, cpus = 5)
  maxAUC_Classifer <- .find_ClassifierWithMaxAUC( sp1_model.denovo )
  # species 2
  sp2_fasta_path <- sprintf("%s/%s/02_merge_assembly/tool/tools_merge.fasta", data_dir, x[2])
  sp2_transfrags <- seqinr::read.fasta(sp2_fasta_path, seqtype = "DNA", as.string = T,
                                       forceDNAtolower = F)
  sp2_seq_type <- fSeq2Genome(i = x[2], fasta_file = sp2_fasta_path)
  sp2_true_positive <- names(sp2_seq_type)[sp2_seq_type%in%c("known", "long")]
  sp2_true_negative <- setdiff(names(sp2_transfrags), sp2_true_positive)
  load(sprintf("part1_denovo_transcripts/RData/%s_feauteMat.RData", x[2]))
  sp2_feature_mat <- feature.merge.df
  rm(feature.merge.df)
  colnames(sp2_feature_mat) <- colnames(sp1_feature_mat)

  pred.prediction.score <- .predictor( method = "randomForest", 
                                       classifier = maxAUC_Classifer$classifier,
                                       featureMat = as.data.frame(scale(fTransforDF(sp2_feature_mat))))
  plot(density(pred.prediction.score))
  Fthresold <- 0.4
  pred.positive <- names(pred.prediction.score)[pred.prediction.score > Fthresold]
  pred.negative <- setdiff(names(pred.prediction.score), pred.positive)
  
  TP <- length(intersect(sp2_true_positive, pred.positive))
  FN <- length(intersect(sp2_true_positive, pred.negative))
  FP <- length(intersect(sp2_true_negative, pred.positive))
  TN <- length(intersect(sp2_true_negative, pred.negative))
  
  eval_vec <- evalPrediction(TP=TP, TN=TN, FP=FP, FN=FN)
  
  name_type1 <- intersect(names(pred.prediction.score), sp2_true_positive)
  name_type2 <- setdiff(names(pred.prediction.score), sp2_true_positive)
  cvListPredictions <- c(pred.prediction.score[name_type1], pred.prediction.score[name_type2])
  cvListLabels <- c( rep(1, length(name_type1)), rep(0, length(name_type2) ) )
  f_pred <- ROCR::prediction(cvListPredictions, cvListLabels)
  f_ROCR <- ROCR::performance(f_pred, measure = "auc")
  eval_vec <- c(eval_vec, as.character(round(f_ROCR@y.values[[1]], digits = 3)))
  names(eval_vec)[length(eval_vec)] <- "AUC"
  eval_vec
}


cross_validation <- function( seed = 1, method = c("randomForest", "svm"), 
                              featureMat, positives, negatives, cross = 5, 
                              cpus = 1, ... ){
  
  call <- match.call()
  
  if( length(method) > 1){
    method <- method[1]
  } 
  
  #sample index for cv
  posSample_cv <- .cvSampleIndex(length(positives), cross = cross, seed = seed)
  negSample_cv <- .cvSampleIndex(length(negatives), cross = cross, seed = seed)
  
  cvRes <- list()
  if( cpus > 1 ) {
    #require(snowfall)
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport("classifier", namespace = "PEA")
    sfExport(".predictor", namespace = "PEA")
    sfExport(".one_cross_validation", namespace = "PEA")
    sfLibrary( "pROC", character.only = TRUE)
    sfLibrary( "e1071", character.only = TRUE)
    sfLibrary( "randomForest", character.only = TRUE )
    
    cvRes <- sfApply( matrix(1:cross, ncol = 1), 1,  .one_cross_validation, method = method, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv, ...)
    sfStop()
  }else {
    for( j in 1:cross ) {
      cvRes[[j]] <- .one_cross_validation( cv = j, method = method, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv, ... )
    }
  }
  cvRes
}


.find_ClassifierWithMaxAUC <- function( cvRes ) {
  
  classifier <- NA
  maxAUC <- 0
  for( i in 1:length(cvRes) ) {
    res <- cvRes[[i]]
    if( res$test.AUC > maxAUC) {
      maxAUC <- res$test.AUC
      classifier <- res$classifier
    }
  }#end for i
  
  return( list(maxAUC = maxAUC, classifier = classifier))
}

.predictor <- function( method = c("randomForest", "svm"), classifier, featureMat ) {
  
  if(length(method) > 1){
    method <- method[1]
  }
  
  if( method == "randomForest") {
    res <- predict(classifier, data.frame(featureMat), type= "vote" )[,"1"]
  }else {
    res <- predict( classifier, data.frame(featureMat), type = "raw") 
  }
  names(res) <- rownames(featureMat)
  res
}


.one_cross_validation <- function( cv, method, featureMat, positives, negatives, posSample_cv, negSample_cv, balanced = TRUE, ratio = 10, ... ) {
  call <- match.call()
  j <- cv
  
  #for train samples
  train_genes_p <- positives[ (posSample_cv$sample_matrix[,j][1:posSample_cv$train_lens[j]] ) ]
  test_genes_p <- positives[ (posSample_cv$sample_matrix[,j][-c(1:posSample_cv$train_lens[j])]) ]
  
  #trained negatives randomly selected, and tested on all negatives
  train_genes_n <- negatives[(negSample_cv$sample_matrix[,j][1:negSample_cv$train_lens[j]] ) ]
  test_genes_n <- negatives[ (negSample_cv$sample_matrix[,j][-c(1:negSample_cv$train_lens[j])]) ]
  
  #select part of train_genes_n
  if( balanced == TRUE ) {
    if( length(train_genes_n) > ratio*length(train_genes_p) ) {
      train_genes_n <- train_genes_n[sample(1:length(train_genes_n), replace=FALSE)[1:(ratio*length(train_genes_p))]]
    }
  }
  
  
  
  obj <- classifier( method = method, featureMat = featureMat, positiveSamples = train_genes_p, negativeSamples = train_genes_n, ... )
  bestmodel <- obj
  
  positives.train.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[train_genes_p,])
  negatives.train.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[train_genes_n,])
  positives.test.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[test_genes_p,])
  negatives.test.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[test_genes_n,])
  
  
  
  train.AUC <- roc( c(rep(1, length(train_genes_p)), rep(0, length(train_genes_n))), 
                    c(positives.train.score, negatives.train.score) )$auc[1]
  test.AUC <- roc( c(rep(1, length(test_genes_p)), rep(0, length(test_genes_n))), 
                   c(positives.test.score, negatives.test.score) )$auc[1]
  
  res <- ( list( positves.train = train_genes_p, negatives.train = train_genes_n, 
                 positives.test = test_genes_p, negatives.test = test_genes_n, 
                 ml = method, classifier = bestmodel, 
                 positives.train.score = positives.train.score,
                 negatives.train.score = negatives.train.score,
                 positives.test.score = positives.test.score,
                 negatives.test.score = negatives.test.score,
                 train.AUC = train.AUC,
                 test.AUC = test.AUC) )
  
  res
}

.cvSampleIndex <- function( len, cross = 5, seed = 1 ) {
  
  cv <- cross
  sample_matrix <- matrix(0, nrow = len, ncol = cv)
  colnames(sample_matrix) <- paste("cv", c(1:cv), sep = "" )
  
  #random samples 
  set.seed(seed)
  index <- sample(1:len, len, replace = FALSE )
  step = floor( len/cv )
  
  start <- NULL
  end <- NULL
  train_lens <- rep(0, cv)
  for( i in c(1:cv) ) {
    start <- step*(i-1) + 1
    end <- start + step - 1
    if( i == cv ) 
      end <- len
    
    train <- index[-c(start:end)]
    test <- index[start:end]
    train_lens[i] <- length(train)
    
    sample_matrix[,i] <- c(train, test)
  }#end for i
  
  return( list( train_lens = train_lens, sample_matrix = sample_matrix))
}


classifier <- function( method = c("randomForest", "svm"), featureMat, positiveSamples, 
                        negativeSamples, ...) {
  call <- match.call()
  
  if( length(method) > 1){
    method <- method[1]
  } 
  
  
  if( is.null(rownames(featureMat) ) )
    stop("Error: no row names (i.e., sample IDs) were assigned for featureMat." )
  if( is.null(colnames(featureMat) ) )
    stop("Error: no colnames were defined for featureMat." )
  
  positiveSamples <- intersect( rownames(featureMat), positiveSamples )
  negativeSamples <- intersect( rownames(featureMat), negativeSamples )
  posLen <- length(positiveSamples)
  negLen <- length(negativeSamples)
  if( posLen == 0 )
    stop("Error: no positive samples included in featureMat." )
  if( negLen == 0 )
    stop("Error: no negative samples were included in featureMat." )
  
  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(positiveSamples, negativeSamples), ] )
  tmpData <- cbind( fmat, label )
  colnames(tmpData) <- c(colnames(fmat), "Class")
  if( method == "randomForest" ) {
    obj <- randomForest(x = fmat, y = factor(label), ... )
  }else{
    obj <- svm(x = fmat, y = factor(label), ... )
  }
  obj
}


fPrediction <- function(positive, negative, unlabeled, featuresmat){
  train.pos <- positive
  train.neg <- negative
  
  cat(length(train.pos), length(train.neg), "\n")
  
  train.feature.mat <- featuresmat[c(train.pos, train.neg), ]
  train.feature.mat <- scale(fTransforDF(train.feature.mat))
  library(snowfall)
  library(randomForest)
  library(pROC)
  model.denovo <- cross_validation(method="randomForest",
                                   featureMat = train.feature.mat,
                                   positives = train.pos, negatives = train.neg,
                                   cross = 5, cpus = 5)
  
  # Ath_model <- model.denovo
  # model.denovo <- PEA::classifier(method = "randomForest", featureMat = train.feature.mat,
  #                            positiveSamples = train.pos, negativeSamples = train.neg)
  # plotROC(model.denovo)
  # posScore <- vector()
  # for (ii in 1:length(model.denovo)) {
  #   posScore <- c(posScore, model.denovo[[ii]]$positives.test.score)
  # }
  # negScore <- vector()
  # for (ii in 1:length(model.denovo)) {
  #   negScore <- c(negScore, model.denovo[[ii]]$negatives.test.score)
  # }
  # Fscore <- vector()
  # for (ii in seq(0, 1, 0.05)) {
  #   measures <- suppressMessages(evalPrediction(threshold = ii, posScore = posScore, negScore = negScore))
  #   cat(ii, measures[5], "\n")
  # }
  
  Fthresold <- 0.5
  maxAUC_Classifer <- .find_ClassifierWithMaxAUC( model.denovo )
  # important.features <- importance(maxAUC_Classifer$classifier)
  # important.features[order(important.features[,1]), ]
  
  pred.feature.mat <- feature.merge.df[unlabeled, ]
  pred.feature.mat <- scale(fTransforDF(pred.feature.mat))
  pred.prediction.score <- .predictor( method = "randomForest", 
                                       classifier = maxAUC_Classifer$classifier,
                                       featureMat = pred.feature.mat)
  
  plot(density(pred.prediction.score))
  
  pred.prediction <- names(pred.prediction.score)[pred.prediction.score > Fthresold]
  cat(length(pred.prediction.score), length(pred.prediction))
  
  return(pred.prediction)
}


ExtractExpression <- function(tool.name){
  exp.mat.list <- list()
  for(rep.name in c("Ath_R1", "Ath_R2", "Ath_R3")){
    tmp.detonate.tpm <- read.table(paste0("data/", rep.name, "/02evaluation/detonate/", 
                                          tool.name, "/rsem_score.score.isoforms.results"),
                                   sep = "\t", header = T)
    rep.name.abbr <- strsplit(rep.name, "_")[[1]][2]
    tmp.detonate.tpm[,1] <- paste0(rep.name.abbr, "_", tmp.detonate.tpm[,1])
    exp.mat.list[[rep.name]] <- tmp.detonate.tpm
  }
  exp.mat <- do.call('rbind', exp.mat.list)
  rownames(exp.mat) <- exp.mat[,1]
  exp.mat
}

GmapRatio <- function(input.names, tree.list){
  map.names <- c(tree.list[[1]]$align, tree.list[[2]]$align, # align  HC_contigs
                 tree.list[[3]]$align, tree.list[[4]]$align)
  length(intersect(input.names, map.names))/length(input.names)*100
}


BaseTheme <- function(){
  theme(panel.border = element_rect(colour = "black"),
        axis.line = element_line(colour = "black",size=0.5),
        axis.text = element_text(colour="black"),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        axis.ticks = element_line(colour = "black",size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
}

fFeaturePlot <- function(fea.mat, pos.tmp, neg.tmp){
  fea.mat.tmp <- fea.mat
  type.tmp <- rep("Unlabel", nrow(fea.mat.tmp))
  type.tmp[rownames(fea.mat.tmp)%in%pos.tmp] <- "Positive"
  type.tmp[rownames(fea.mat.tmp)%in%neg.tmp] <- "Negative"
  fea.mat.tmp$type <- type.tmp
  fea.mat.tmp$type <- factor(fea.mat.tmp$type, 
                             levels = intersect(c("Positive", "Negative", "Unlabel"), unique(type.tmp)))
  fea.mat.tmp$name <- rownames(fea.mat.tmp)
  fea.df.tmp <- reshape2::melt(fea.mat.tmp, id=c("type","name"))
  table(fea.mat.tmp$type)
  library(ggplot2)
  ggplot(fea.df.tmp, aes(x=value, color=type)) +
    geom_line(stat = "density", size=0.9, alpha=0.8) + 
    facet_wrap(variable~., scales = "free", ncol = 5) +
    theme_bw(base_size = 14) +
    scale_color_manual(values = c("#8BC33E", "#00ACED", "#F5911F", "#EB1B24")[1:length(unique(type.tmp))]) +
    BaseTheme()
}


fPCA <- function(matrix){
  data2.pca <- prcomp(matrix, center=F, scale=F)
  summary(data2.pca)
  plot(data2.pca$x[,1], data2.pca$x[,2], col=soilpalettes::soil_palette("redox2",7)[factor(contig.tree$class)], main="PCA")
}

# library(data.table)
# library(dplyr)
# 


.getPvalue <- function(inputVec, mappedInput, mappedRIP){
  testMat <- matrix(c(as.numeric(inputVec[8]),
                      mappedInput,
                      as.numeric(inputVec[9]),
                      mappedRIP), nrow = 2, ncol = 2)
  p.value <- fisher.test(x = testMat)$p
  ratio <- log2(((as.numeric(inputVec[8]) + 1)*mappedRIP)/((as.numeric(inputVec[9]) + 1)*mappedInput))
  res <- c(ratio, p.value)
  res
}

