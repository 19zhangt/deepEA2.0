#! /usr/bin/Rscript

# Setting environment variables
options(stringsAsFactors = F)
setwd('/storage_server/home/zhangt/freePEA/mlPEA_denovo/results')

library(seqinr)


Input_df <- read.table("part0/ppa_kmer51/Input_union")
IP_df <- read.table("part0/ppa_kmer51/IP_union")

# range(Input_df$V2)
# range(IP_df$V2)

Input_totalread <- 99270040/4
IP_totalread <- 116905048/4


tmp_value <- (IP_df$V2/IP_totalread)/(Input_df$V2/Input_totalread)
percen_value <- vector()
for(i in 0:2){
  if(i != 2){
    percen_value <- c(percen_value, sum(tmp_value > i & tmp_value < i+1))
  } else {
    percen_value <- c(percen_value, sum(tmp_value > i))
  }
  cat(i, "\n")
}

barplot(percen_value/10000, las=1)

source("part1_denovo_transcripts/functions.R")


fPieRatio(f_itable = percen_value, main = "Enrichment_ratio")

pie_data <-  data.frame("name"= c("0-1", "1-2", ">2"), "Freq"=percen_value)
rownames(pie_data) <- c("0-1", "1-2", ">2")
pie_data <- pie_data[order(pie_data[,2], decreasing = T), ]
colnames(pie_data)[1] <- "Var1"
rownames(pie_data) <- pie_data[,1]
pie_data$Var1 <- factor(pie_data$Var1, levels = rownames(pie_data))
pie_data$per <- round(pie_data$Freq/sum(pie_data$Freq), digits = 4)
pie_data$label <- paste0(pie_data$per*100, "%")

p <- ggplot(data = pie_data,
             mapping = aes(x = "Count", y = per, fill = Var1)) +
        # scale_fill_manual(values = c('#BDD1C5','#4D6372','#E7D86E','#558ED5',"#A984A8"), guide = FALSE)
        geom_bar(stat = 'identity', position = 'fill', width = 0.5) +
        labs(x = '', y = '', title = '') +
        coord_polar("y", direction = -1) +
        theme_void(base_size = 12) +
        theme(axis.text = element_blank(), axis.ticks = element_blank()) +
        scale_fill_manual(values = rev(soilpalettes::soil_palette("redox2", 
                                                                  length(unique(pie_data$Var1))))) +
        ggtitle("Fold Enrichment") + theme(plot.title = element_text(hjust = 0.5)) +
        geom_text(aes(x=1.3, y = sum(per) - cumsum(per) + per/2, label=label)) + 
        guides(fill = guide_legend(title = NULL))


# enrichment
enrich_index  <- (IP_df$V2/IP_totalread)/(Input_df$V2/Input_totalread) > 2
table(enrich_index)

remove_df <- Input_df[!enrich_index, ]
extract_remove <- sample(1:nrow(remove_df), 10000)
remove_gc <- sapply(remove_df$V1[extract_remove], function(x){GC(s2c(x))*100})
pdf("part0/kmer1.pdf", width = 7, height = 5)
p

plot(density(remove_gc), main="", xlab="GC content", col="grey", lwd=2, las=1)

union_df <- cbind(Input_df[enrich_index, ], IP_df[enrich_index, 2])
extract_union <- sample(1:nrow(union_df), 10000)
union_gc <- sapply(union_df$V1[extract_union], function(x){GC(s2c(x))*100})
lines(density(union_gc), lwd=2)


plot(density(log2(union_df$V2)), main="", xlab="k-mer count", col="grey", lwd=2, las=1)
lines(density(log2(union_df$`IP_df[enrich_index, 2]`)), lwd=2)

dev.off()

# nrow(union_df)
# ## fisher test
# pvalue <- apply(union_df, 1, function(x){
#   fisher.test(matrix(c(as.numeric(x[2]), Input_totalread, 
#                        as.numeric(x[3]), IP_totalread), nrow = 2, byrow = T))$p.value
# })

## 0 
# GC_content <- apply(Input_df, 1, function(x){seqinr::GC(seqinr::s2c(x[1]))})


## 1
# kmer_index <- IP_df$V2/Input_df$V2 > 2
# table(kmer_index)
# 
# plot(log2(IP_df$V2[kmer_index]), log2(IP_df$V2/Input_df$V2)[kmer_index])
# 
# plot(density(IP_df$V2/Input_df$V2), xlim=c(0,10))
# 
# ## 2
# index2 <- IP_df$V2 > 100
# table(index2)
# 
# ## 3
# index3 <- order(IP_df$V2, decreasing = T)[1:100000]
# index3_gc <- apply(Input_df[index3, ], 1, function(x){seqinr::GC(seqinr::s2c(x[1]))*100})
# 
# index3 <- index3[index3_gc<30]
# xvalue3 <- log2(Input_df$V2[index3])
# yvalue3 <- log2(IP_df$V2[index3])
# 
# plot(xvalue3, yvalue3, cex = 0.5, xlim=c(0, max(c(xvalue3, yvalue3))), 
#      ylim=c(0, max(c(xvalue3, yvalue3))))
# 
# 
# ## 4
# index4 <- order(Input_df$V2, decreasing = T)[1:100000]
# xvalue4 <- log2(Input_df$V2[index4])
# yvalue4 <- log2(IP_df$V2[index4])
# 
# plot(xvalue4, yvalue4, cex = 0.5, xlim=c(0, max(c(xvalue4, yvalue4))), 
#      ylim=c(0, max(c(xvalue4, yvalue4))))
