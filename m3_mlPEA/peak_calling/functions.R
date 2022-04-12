

fCallRelativePeak <- function(input_sp){
  ## use median as background
  library(GenomicRanges)
  library(future)
  library(future.apply)
  library(progress)
  library(rlang)
  
  ratio <- -1
  level <- 0.05
  readsCount <- 10
  concatenate = 4
  
  # call peak
  .winCount <- function(x, y, z){
    data.frame(
      "Position" = subset(z, !duplicated(x)),
      "window" = subset(x, !duplicated(x)),
      "curWave" = as.numeric(by(y, factor(x, levels = unique(x)), mean))
    )
  }
  
  .intervalCov <- function(dzFile){
    baseCov <- fread(file = dzFile, sep = "\t", header = F, stringsAsFactors = F)
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
  
  # sequences and its lengths of transfrags
  fasta_file <- sprintf("%s/%s/02_merge_assembly/tool/tools_merge.fasta", data_dir, input_sp)
  transfrags <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = T, forceDNAtolower = F)
  seq_length <- sapply(transfrags, nchar)
  
  rangelist <- list()
  
  for(rep_num in 1:3){
    cat(sprintf("================== Running: rep %s  =================", rep_num))
    mappedInput <- as.numeric(system(command= sprintf('head -n 1 %s/%s/03_evaluation/hisat2/%s_R%s_mapping.info | grep -oP "[0-9]+"',
                                                      data_dir, input_sp, input_sp, rep_num), intern = T))
    mappedRIP <- as.numeric(system(command= sprintf('head -n 1 %s/%s/03_evaluation/hisat2_IP/%s_R%s_mapping.info | grep -oP "[0-9]+"',
                                                    data_dir, input_sp, input_sp, rep_num), intern = T))
    inputfile <- sprintf("%s/%s/03_evaluation/coverage/%s_R%s.input.txt",
                         data_dir, input_sp, input_sp, rep_num)
    ipfile <- sprintf("%s/%s/03_evaluation/coverage/%s_R%s.ip.txt",
                      data_dir, input_sp, input_sp, rep_num)
    
    input <- .intervalCov(dzFile = inputfile)
    RIP <- .intervalCov(dzFile = ipfile)
    
    if(!all(names(input) == names(RIP))){
      cat("Note: The chromosomes in the input and RIP are not consistent!\n",
          "the interactions will be used!")
      interNames <- intersect(names(input), names(RIP))
      input <- input[interNames]
      RIP <- RIP[interNames]
    }
    
    resList <- vector("list", length = length(RIP))
    names(resList) <- names(RIP)
    
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
    
    cl <- makeClusterPSOCK(20) #availableCores()
    plan(cluster, workers = cl)
    system.time(resList <- future_lapply(names(input), function(x){
      curMat <- merge(input[[x]], RIP[[x]], by = 'window', all=TRUE)
      curMat$curWave.x[is.na(curMat$curWave.x)] <- 0
      curMat$curWave.y[is.na(curMat$curWave.y)] <- 0
      curMat[, 'windowave.input'] <- round(curMat$curWave.x)
      curMat[, 'windowave.RIP'] <- round(curMat$curWave.y)
      curMat <- curMat[curMat$windowave.RIP != 0, ]

      curMat <- curMat[curMat$windowave.RIP > max(curMat$windowave.RIP)/20, ]
      
      median_input <- median(curMat$windowave.input)
      median_ip <- median(curMat$windowave.RIP)
      
      curPvalue <- rep(NA, nrow(curMat))
      curRatio <- rep(NA, nrow(curMat))
      curFDR <- rep(NA, nrow(curMat))
      curMat <- cbind(curMat, curPvalue, curRatio, curFDR)
      # cat(dim(curMat), "\n")
      pvalue <- t(apply(curMat, 1, .getPvalue, mappedInput = median_input, mappedRIP = median_ip))
      curMat$curPvalue <- as.numeric(pvalue[,2])
      curMat$curRatio <- as.numeric(pvalue[,1])
      curMat$curFDR <- p.adjust(curMat$curPvalue, "fdr")
      curMat
    }))
    # future:::ClusterRegistry("stop")
    
    parallel::stopCluster(cl)
    
    .findContinuous <- function(inputVec){
      Breaks <- c(0, which(diff(inputVec) != 1), length(inputVec))
      res <- lapply(seq(length(Breaks) - 1),
                    function(i) inputVec[(Breaks[i] + 1):Breaks[i+1]])
      res
    }
    
    resMat <- do.call(what = rbind, args = resList)
    resMat <- subset(resMat, resMat$curRatio < ratio & resMat$curFDR < level &
                       resMat$curWave.y >= readsCount)
    
    tt <- by(resMat$window, factor(x = resMat$Chr.y, levels = unique(resMat$Chr.y)), .findContinuous)
    resPeaks <- NULL
    
    pb <- progress_bar$new(total = length(tt), clear = FALSE)
    
    for(i in 1:length(tt)){
      pb$tick()
      curList <- tt[[i]]
      curLen <- unlist(lapply(curList, length))
      curList <- curList[which(curLen >= concatenate)]
      if(length(curList) == 0){
        next
      }
      Start <- unlist(lapply(curList, function(x) (x[1]-1)*25+1))
      End <- unlist(lapply(curList, function(x) (x[length(x)]*25)))
      curMat <- subset(resMat, resMat$Chr.y == names(tt)[i])
      curFDR <- lapply(curList, function(x)  curMat$curFDR[match(x, curMat$window)])
      meanFDR <- unlist(lapply(curFDR, mean))
      maxFDR <- unlist(lapply(curFDR, max))
      minFDR <- unlist(lapply(curFDR, min))
      curRatio <- lapply(curList, function(x)  curMat$curRatio[match(x, curMat$window)])
      meanRatio <- unlist(lapply(curRatio, mean))
      maxRatio <- unlist(lapply(curRatio, max))
      minRatio <- unlist(lapply(curRatio, min))
      windowNumber <- unlist(lapply(curList, length))
      curPeaks <- cbind(names(tt)[i], Start, End, windowNumber,
                        meanFDR, maxFDR, minFDR,
                        meanRatio, maxRatio, minRatio)
      resPeaks <- rbind(resPeaks, curPeaks)
    }
    colnames(resPeaks) <- c("Chromosome", "Start(1-based)", "End", "Bin number",
                            "Mean FDR", "Max FDR", "Minimum FDR",
                            "Mean Ratio", "Max Ratio", "Minimum Ratio")
    
    for(i in 1:nrow(resPeaks)){
      if(as.numeric(resPeaks[i, 3]) > seq_length[resPeaks[i, 1]]){
        resPeaks[i, 3] <- seq_length[resPeaks[i, 1]]
      }
    }
    cat(nrow(resPeaks), "\n")
    ## peaks
    rangelist[[sprintf("rep%s", rep_num)]] <- GRanges(seqnames = resPeaks[,1], 
                                    ranges = IRanges(start = as.numeric(resPeaks[,2]), 
                                                     end = as.numeric(resPeaks[,3])))
  }
  
  peak_grangeslist <- GRangesList(rangelist)
  peak_coverage <- coverage(peak_grangeslist)
  covered_ranges <- IRanges::slice(peak_coverage, lower=2, rangesOnly = T)
  rangelist[["merged_peaks"]] <- GRanges(covered_ranges)

  gc()
  rangelist
}


fPEATransciptPeak <- function(input_sp){
  ## use total reads as background
  library(GenomicRanges)
  library(future)
  library(future.apply)
  library(progress)
  library(rlang)
  
  ratio <- -1
  level <- 0.05
  readsCount <- 10
  concatenate = 4
  
  # call peak
  .winCount <- function(x, y, z){
    data.frame(
      "Position" = subset(z, !duplicated(x)),
      "window" = subset(x, !duplicated(x)),
      "curWave" = as.numeric(by(y, factor(x, levels = unique(x)), mean))
    )
  }
  
  .intervalCov <- function(dzFile){
    baseCov <- fread(file = dzFile, sep = "\t", header = F, stringsAsFactors = F)
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
  
  
  # sequences and its lengths of transfrags
  fasta_file <- sprintf("%s/%s/02_merge_assembly/tool/tools_merge.fasta", data_dir, input_sp)
  transfrags <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = T, forceDNAtolower = F)
  seq_length <- sapply(transfrags, nchar)
  
  rangelist <- list()
  library(data.table)
  for(rep_num in 1:3){
    cat(sprintf("================== Running: rep %s  =================", rep_num))
    mappedInput <- as.numeric(system(command= sprintf('head -n 1 %s/%s/03_evaluation/hisat2/%s_R%s_mapping.info | grep -oP "[0-9]+"',
                                                      data_dir, input_sp, input_sp, rep_num), intern = T))
    mappedRIP <- as.numeric(system(command= sprintf('head -n 1 %s/%s/03_evaluation/hisat2_IP/%s_R%s_mapping.info | grep -oP "[0-9]+"',
                                                    data_dir, input_sp, input_sp, rep_num), intern = T))
    inputfile <- sprintf("%s/%s/03_evaluation/coverage/%s_R%s.input.txt",
                         data_dir, input_sp, input_sp, rep_num)
    ipfile <- sprintf("%s/%s/03_evaluation/coverage/%s_R%s.ip.txt",
                      data_dir, input_sp, input_sp, rep_num)
    
    input <- .intervalCov(dzFile = inputfile)
    RIP <- .intervalCov(dzFile = ipfile)
    
    if(!all(names(input) == names(RIP))){
      cat("Note: The chromosomes in the input and RIP are not consistent!\n",
          "the interactions will be used!")
      interNames <- intersect(names(input), names(RIP))
      input <- input[interNames]
      RIP <- RIP[interNames]
    }
    
    resList <- vector("list", length = length(RIP))
    names(resList) <- names(RIP)
    
    cl <- makeClusterPSOCK(10) #availableCores()
    plan(cluster, workers = cl)
    system.time(resList <- future_lapply(names(input), function(x){
      curMat <- merge(input[[x]], RIP[[x]], by = 'window', all=TRUE)
      curMat$curWave.x[is.na(curMat$curWave.x)] <- 0
      curMat$curWave.y[is.na(curMat$curWave.y)] <- 0
      curMat[, 'windowave.input'] <- round(curMat$curWave.x)
      curMat[, 'windowave.RIP'] <- round(curMat$curWave.y)
      curMat <- curMat[curMat$windowave.RIP != 0, ]
      
      curPvalue <- rep(NA, nrow(curMat))
      curRatio <- rep(NA, nrow(curMat))
      curFDR <- rep(NA, nrow(curMat))
      curMat <- cbind(curMat, curPvalue, curRatio, curFDR)
      # cat(dim(curMat), "\n")
      pvalue <- t(apply(curMat, 1, .getPvalue, mappedInput = mappedInput, mappedRIP = mappedRIP))
      curMat$curPvalue <- as.numeric(pvalue[,2])
      curMat$curRatio <- as.numeric(pvalue[,1])
      curMat$curFDR <- p.adjust(curMat$curPvalue, "fdr")
      curMat
    }))
    # future:::ClusterRegistry("stop")
    
    parallel::stopCluster(cl)
    
    .findContinuous <- function(inputVec){
      Breaks <- c(0, which(diff(inputVec) != 1), length(inputVec))
      res <- lapply(seq(length(Breaks) - 1),
                    function(i) inputVec[(Breaks[i] + 1):Breaks[i+1]])
      res
    }
    
    resMat <- do.call(what = rbind, args = resList)
    resMat <- subset(resMat, resMat$curRatio < ratio & resMat$curFDR < level &
                       resMat$curWave.y >= readsCount)
    
    tt <- by(resMat$window, factor(x = resMat$Chr.y, levels = unique(resMat$Chr.y)), .findContinuous)
    resPeaks <- NULL
    
    pb <- progress_bar$new(total = length(tt), clear = FALSE)
    
    for(i in 1:length(tt)){
      pb$tick()
      curList <- tt[[i]]
      curLen <- unlist(lapply(curList, length))
      curList <- curList[which(curLen >= concatenate)]
      if(length(curList) == 0){
        next
      }
      Start <- unlist(lapply(curList, function(x) (x[1]-1)*25+1))
      End <- unlist(lapply(curList, function(x) (x[length(x)]*25)))
      curMat <- subset(resMat, resMat$Chr.y == names(tt)[i])
      curFDR <- lapply(curList, function(x)  curMat$curFDR[match(x, curMat$window)])
      meanFDR <- unlist(lapply(curFDR, mean))
      maxFDR <- unlist(lapply(curFDR, max))
      minFDR <- unlist(lapply(curFDR, min))
      curRatio <- lapply(curList, function(x)  curMat$curRatio[match(x, curMat$window)])
      meanRatio <- unlist(lapply(curRatio, mean))
      maxRatio <- unlist(lapply(curRatio, max))
      minRatio <- unlist(lapply(curRatio, min))
      windowNumber <- unlist(lapply(curList, length))
      curPeaks <- cbind(names(tt)[i], Start, End, windowNumber,
                        meanFDR, maxFDR, minFDR,
                        meanRatio, maxRatio, minRatio)
      resPeaks <- rbind(resPeaks, curPeaks)
    }
    colnames(resPeaks) <- c("Chromosome", "Start(1-based)", "End", "Bin number",
                            "Mean FDR", "Max FDR", "Minimum FDR",
                            "Mean Ratio", "Max Ratio", "Minimum Ratio")
    
    for(i in 1:nrow(resPeaks)){
      if(as.numeric(resPeaks[i, 3]) > seq_length[resPeaks[i, 1]]){
        resPeaks[i, 3] <- seq_length[resPeaks[i, 1]]
      }
    }
    cat(nrow(resPeaks), "\n")
    ## peaks
    rangelist[[sprintf("rep%s", rep_num)]] <- GRanges(seqnames = resPeaks[,1], 
                                    ranges = IRanges(start = as.numeric(resPeaks[,2]), 
                                                     end = as.numeric(resPeaks[,3])))
  }
  
  peak_grangeslist <- GRangesList(rangelist)
  peak_coverage <- coverage(peak_grangeslist)
  covered_ranges <- IRanges::slice(peak_coverage, lower=2, rangesOnly = T)
  rangelist[["merged_peaks"]] <- GRanges(covered_ranges)
  
  gc()
  rangelist
}


fDiffBindCount <- function(input_sp){
  # read count using diffbind
  ipbam_path <- sprintf("/storage_server/home/zhangt/freePEA/data/%s/03_evaluation/hisat2_IP/%s", input_sp, input_sp)
  inputbam_path <- sprintf("/storage_server/home/zhangt/freePEA/data/%s/03_evaluation/hisat2/%s", input_sp, input_sp)
  system("mkdir -p part2_peak_calling/Peak")
  bed_path <- sprintf("part2_peak_calling/Peak/%s_peak.bed", input_sp)
  out_path <- sprintf("part2_peak_calling/Peak/%s_peak_count.txt", input_sp)
  
  if (!file.exists(out_path)) {
    diffbind_conf <- sprintf("part2_peak_calling/Peak/%s_diffBind_conf.csv", input_sp) 
    
    writeLines(sprintf("SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
Sample1,Leaf,Leaf1,IP1,1,%s_R1.bam,Input1,%s_R1.bam,%s,bed
Sample2,Leaf,Leaf2,IP2,2,%s_R2.bam,Input2,%s_R2.bam,%s,bed
Sample3,Leaf,Leaf3,IP3,3,%s_R3.bam,Input3,%s_R3.bam,%s,bed",
                       ipbam_path, inputbam_path, bed_path,
                       ipbam_path, inputbam_path, bed_path,
                       ipbam_path, inputbam_path, bed_path),
               con = diffbind_conf )
    
    suppressMessages(library(DiffBind))
    suppressMessages(library(rtracklayer))
    
    peakMerge <- import(bed_path)
    tamoxifen <- dba(sampleSheet = diffbind_conf)
    txfCount <- dba.count(tamoxifen, peaks = peakMerge, summits=FALSE, minOverlap = 0)
    
    peak_count_df <- data.frame("Chr"= txfCount$peaks[[1]][,1], "Start"= txfCount$peaks[[1]][,2]-1, "End" = txfCount$peaks[[1]][,3],
                                "IP1"= round(as.numeric(txfCount$peaks[[1]][,5]), 2), "IP2"= round(as.numeric(txfCount$peaks[[2]][,5]), 2), 
                                "IP3"= round(as.numeric(txfCount$peaks[[3]][,5]), 2), "Input1"= round(as.numeric(txfCount$peaks[[1]][,7]), 2), 
                                "Input2"= round(as.numeric(txfCount$peaks[[2]][,7]), 2), "Input3"= round(as.numeric(txfCount$peaks[[3]][,7]), 2)
    )
    
    peak_count_df$ratio <- (peak_count_df$IP1+peak_count_df$IP2+peak_count_df$IP3+1)/(peak_count_df$Input1+peak_count_df$Input2+peak_count_df$Input3+1)
    
    peak_count_df$diff <- (peak_count_df$IP1+peak_count_df$IP2+peak_count_df$IP3)-(peak_count_df$Input1+peak_count_df$Input2+peak_count_df$Input3)
    
    write.table(peak_count_df, file = out_path, quote = F, row.names = F, sep = "\t")
    
  } else {
    peak_count_df <- read.table(out_path, header = T, sep = "\t")
  }
  
  peak_count_df
}


fCoveragePlot <- function(input_sp, peak_names){
  # coverage plot in input and ip samples
  rep_num <- 1
  
  inputfile <- sprintf("%s/%s/03_evaluation/coverage/%s_R%s.input.txt",
                       data_dir, input_sp, input_sp, rep_num)
  ipfile <- sprintf("%s/%s/03_evaluation/coverage/%s_R%s.ip.txt",
                    data_dir, input_sp, input_sp, rep_num)
  mappedInput <- as.numeric(system(command= sprintf('head -n 1 %s/%s/03_evaluation/hisat2/%s_R%s_mapping.info | grep -oP "[0-9]+"',
                                                    data_dir, input_sp, input_sp, rep_num), intern = T))
  mappedRIP <- as.numeric(system(command= sprintf('head -n 1 %s/%s/03_evaluation/hisat2_IP/%s_R%s_mapping.info | grep -oP "[0-9]+"',
                                                  data_dir, input_sp, input_sp, rep_num), intern = T))
  tmpPeakName <- peak_names
  
  resPeaks <- read.table(sprintf("part2_peak_calling/Peak/%s_peak.bed", input_sp), sep = "\t")
  
  peak_loc <- resPeaks[resPeaks[,1] == tmpPeakName, ]
  
  system("mkdir -p tmp_cov")
  system(command = paste0("grep ", tmpPeakName, " ", inputfile," > tmp_cov/tmp_input_", tmpPeakName, ".txt"))
  system(command = paste0("grep ", tmpPeakName," ", ipfile, " > tmp_cov/tmp_ip_", tmpPeakName, ".txt"))
  input_count <- read.table(paste0("tmp_cov/tmp_input_", tmpPeakName, ".txt"), sep = "\t")
  ip_count <- read.table(paste0("tmp_cov/tmp_ip_", tmpPeakName, ".txt"), sep = "\t")
  input_count$V3 <- input_count$V3/mappedInput*10^6
  ip_count$V3 <- ip_count$V3/mappedRIP*10^6
  
  cov.data <- data.frame(genome_location=c(ip_count$V2, input_count$V2), 
                         value=c(ip_count$V3, input_count$V3),
                         Group = factor( rep(c("IP", "Input"), 
                                             c(length(ip_count$V2), 
                                               length(input_count$V2)) ), 
                                         levels = c("IP", "Input") ))
  
  fasta_file <- sprintf("%s/%s/02_merge_assembly/tool/tools_merge.fasta", data_dir, input_sp)
  transfrags <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = T, forceDNAtolower = F)
  tmp_seq <- transfrags[[peak_names]]
  motif_site <- gregexpr('[AGag][AGag][Aa][Cc][ACTact]', tmp_seq)[[1]]
  motif_site <- as.data.frame(motif_site)
  
  # peak_loc2 <- ath_peaks22[ath_peaks22$seqnames == peak_names, ]
  
  tmp_max <- max(cov.data$value)
  
  ggplot(data = cov.data)+
    geom_line(aes(x = genome_location, y = value, colour = Group))+
    # geom_ribbon(aes(ymax = value,ymin=0,fill=Group), alpha = 0.4) +
    labs(y = "Normalized coverage", x = "" ) +
    annotate("rect", xmin = as.numeric(peak_loc[2]), xmax = as.numeric(peak_loc[3]), 
             ymin = -0.05*tmp_max, ymax = -0.03*tmp_max, 
             alpha = 0.99, colour = NA) +
    # annotate("rect", xmin = as.numeric(peak_loc2[2]), xmax = as.numeric(peak_loc2[3]), 
    #          ymin = -0.02*tmp_max, ymax = 0, 
    #          alpha = 0.99, colour = NA) +
    annotate("text", x = (as.numeric(peak_loc[2]) + as.numeric(peak_loc[3]))/2, 
             y= -0.01*tmp_max, label = "Peak", size = 3) +
    geom_segment(data = motif_site, aes(x = motif_site, y = rep(-0.05*tmp_max, length(motif_site)), 
                     xend = motif_site, yend = rep(-0.03*tmp_max, length(motif_site))), color = 'red') +
    scale_x_continuous(breaks = round(seq(min(cov.data$genome_location), 
                                          max(cov.data$genome_location), 
                                          by = ((max(cov.data$genome_location)-min(cov.data$genome_location))/5) )),
                       expand = c(0,0))+
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black", size = 0.5),
          axis.title = element_text(color = "black", size = 12),
          axis.text = element_text(colour = "black", size = 12) ) + 
    scale_y_continuous(expand = c(0, 0)) + labs(title = peak_names) + 
    scale_color_manual(values = c("#70AD45", "#4594D1"))
}

## need ml_result_list
fCountSummary <- function(input_sp){
  ## Summary
  load(paste0("part2_peak_calling/RData/", input_sp, "_peak.RData"))
  peaks_df <- as.data.frame(eval(parse(text = paste0(tolower(input_sp), "_peaks"))))
  peaks_df$seqnames <- as.character(peaks_df$seqnames)
  pu_remian <- ml_result_list[[input_sp]]$pu_remain
  peaks_df <- peaks_df[peaks_df$seqnames%in%pu_remian, ]
  sum_val <- sprintf("Transfrags: %s, with peak: %s, Number of peaks: %s",
                     length(pu_remian), length(unique(peaks_df$seqnames)), nrow(peaks_df))
  cat(sum_val, "\n")
  
  mapping_info <- fSeq2Genome(i = input_sp, 
                              fasta_file = sprintf("%s/%s/02_merge_assembly/tool/tools_merge.fasta", 
                                                   data_dir, input_sp))
  
  plot1 <- fPieRatio(f_itable = table(mapping_info$out[pu_remian]), main = paste0(input_sp, "_all"))
  plot2 <- fPieRatio(f_itable = table(mapping_info$out[unique(peaks_df$seqnames)]), main = paste0(input_sp, "_with_peak"))
  
  transfrags_withintron <- read.table("../../data/Ath/02_merge_assembly/tool/tools_merge.fasta_gffcompare_with_intron.txt")
  table(mapping_info$out[pu_remian])
  table(mapping_info$out[pu_remian])/length(pu_remian)*100
  
  tmp_names <- unique(peaks_df$seqnames)
  for(tt in c("known","long_m", "long_k", "short_j/n", "short_c/e")){
    cat(tt, "\n")
    cat(sum(mapping_info$out[tmp_names]%in%tt), "\n")
    cat(sum(names(mapping_info$out[tmp_names])[mapping_info$out[tmp_names]%in%tt]%in%transfrags_withintron$V1), "\n")
    cat(sum(!names(mapping_info$out[tmp_names])[mapping_info$out[tmp_names]%in%tt]%in%transfrags_withintron$V1), "\n")
    # 
    # setdiff(names(mapping_info$out[pu_remian])[mapping_info$out[pu_remian]%in%tt], 
    #           transfrags_withintron$V1)[1:2]
  }
  
  patchwork::wrap_plots(plot1 + plot2, guides = "collect")
  
  output_list <- list()
  output_list[['gff_table']] <- mapping_info
  output_list[['transfrags']] <- pu_remian
  output_list[['fragswithpeak']] <- unique(peaks_df$seqnames)
  output_list
  
  ## 
  multi_isoform_genes <- table(mapping_info$gff_table[,2])
  multi_isoform_genes <- names(multi_isoform_genes)[multi_isoform_genes>1]
  filter_mapinfo <- mapping_info$gff_table[!mapping_info$gff_table[,2]%in%multi_isoform_genes, ]
  single_isoform_trans <- rownames(filter_mapinfo)
}
