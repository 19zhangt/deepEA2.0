#!/usr/bin/env Rscript

# loading runtime

species <- c("Ath", "Pvu", "Gar", "Sbi", "Ata", "Osa", "Ppa", "Gma", "Ghi", "Zma", "Tdi", "Tae")

library(dplyr)
library(ggplot2)
data_dir <- c("/storage_server/home/zhangt/freePEA/data")

## Trinity
runtime_list <- list()
for (sp in species){
  for(rep in c("R1", "R2", "R3")){
    runtime <- read.csv(sprintf("%s/%s/01_assembly/01_%s_assembly/runtime.log", 
                                data_dir, sp, rep), header = F)
    runtime_list[[paste0(sp, rep)]] <- runtime
  }
}
runtime_mat <- do.call('rbind', runtime_list)
runtime_mat$V4 <- round(as.numeric(runtime_mat$V4/3600), digits = 2)

## rnaSPAdes and TransABySS
rt_load <- read.csv(sprintf("%s/../time.log", data_dir), header = F)
rt_load$time <- apply(rt_load, 1, function(x){
  x[4] <- gsub(" ", "", x[4])
  if(grepl(pattern = "\\|", x[4])){
    tmp_time <- do.call('rbind', strsplit(strsplit(x[4], split = "\\|")[[1]], split = ":"))
    tmp_time <- apply(tmp_time,2,function(x){sum(as.numeric(x))})
    tmp_time <- tmp_time[1] + tmp_time[2]/60
  }else{
    tmp_time <- as.numeric(strsplit(x[4], split = "h|m|s")[[1]])
    tmp_time <- tmp_time[1]+tmp_time[2]/60
  }
  round(tmp_time, digits = 2)
})
rt_load <- rt_load[, -4]
colnames(rt_load) <- colnames(runtime_mat)

tool_list <- c("TransABySS", "rnaSPAdes", "Trinity")
runtime_merge <- rbind(runtime_mat, rt_load)

runtime_merge <- runtime_merge %>% group_by(V1, V2) %>% 
  summarise("mean"=mean(V4), 
            "sd"=sd(V4), 
            "se"=sd(V4)/sqrt(length(V4))) %>% 
  as.data.frame()

runtime_merge <- runtime_merge[runtime_merge$V1%in%tool_list, ]
runtime_merge<- runtime_merge[runtime_merge$V2%in%species, ]
runtime_merge$V1 <- factor(runtime_merge$V1, levels = tool_list)
runtime_merge$V2 <- factor(runtime_merge$V2, levels = intersect(species,unique(runtime_merge$V2)))


ggplot(data = runtime_merge, aes(x=V2, y=mean, fill=V1))+
  geom_bar(stat='identity', position="dodge") +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, position = position_dodge(0.9)) +
  theme_bw(base_size = 12) +
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + ylab("Runtime (h)") + 
  scale_fill_manual(values = rev(soilpalettes::soil_palette("redox2", length(unique(runtime_merge$V1)))))

system(command = "mkdir -p part1_denovo_transcripts/results_plot")
ggsave("part1_denovo_transcripts/results_plot/runtime_summary.pdf", width = 8, height = 4)


## Number of transcripts 
df_transfrag <- matrix(nrow = length(tmp_tool_list), ncol = 5)
tmp_tool_list <- rep(tool_list, 12)
tmp_species <- rep(species, each = 3)
df_transfrag[,1] <- tmp_species
df_transfrag[,2] <- tmp_tool_list

for(i in unique(tmp_species)){
  t_index1 <- df_transfrag[,1]==i&df_transfrag[,2]=="Trinity"
  for(j in 1:3){
    t_trinity <- system(command = sprintf('grep ">" ../../data/%s/01_assembly/01_R%s_assembly/Trinity.fasta | wc -l',
                                          i, j), intern = T)
    df_transfrag[t_index1, j+2] <- t_trinity
  }
  t_index2 <- df_transfrag[,1]==i&df_transfrag[,2]=="rnaSPAdes"
  for(j in 1:3){
    t_trinity <- system(command = sprintf('grep ">" ../../data/%s/01_assembly/01_R%s_assembly/rnaSPAdes/transcripts.fasta | wc -l',
                                          i, j), intern = T)
    df_transfrag[t_index2, j+2] <- t_trinity
  }
  t_index3 <- df_transfrag[,1]==i&df_transfrag[,2]=="TransABySS"
  for(j in 1:3){
    t_trinity <- system(command = sprintf('grep ">" ../../data/%s/01_assembly/01_R%s_assembly/TransABySS/transabyss-final.fa | wc -l',
                                          i, j), intern = T)
    df_transfrag[t_index3, j+2] <- t_trinity
  }
}

colnames(df_transfrag) <- c("species", "tool", "R1", "R2", "R3")
df_transfrag <- as.data.frame(df_transfrag, stringsAsFactors = F)
df_transfrag <- reshape2::melt(df_transfrag, id.vars = c("species", "tool"))
df_transfrag$value <- as.numeric(df_transfrag$value)

df_transfrag <- df_transfrag %>% group_by(species, tool) %>% 
  summarise("mean"=mean(value), 
            "sd"=sd(value), 
            "se"=sd(value)/sqrt(length(value))) %>% 
  as.data.frame()

df_transfrag$tool <- factor(df_transfrag$tool, levels = tool_list)
df_transfrag$species <- factor(df_transfrag$species, 
                               levels = intersect(species,unique(df_transfrag$species)))

ggplot(data = df_transfrag, aes(x=species, y=mean, fill=tool))+
  geom_bar(stat='identity', position="dodge") +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, position = position_dodge(0.9)) +
  theme_bw(base_size = 12) +
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + ylab("Number of transfrags") + 
  scale_fill_manual(values = rev(soilpalettes::soil_palette("redox2", length(unique(df_transfrag$tool)))))

ggsave("part1_denovo_transcripts/results_plot/numberoftransfrags.pdf", width = 8, height = 4)
