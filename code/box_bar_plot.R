library(ggplot2)

#box plot
p > ggplot(qpcrgene, aes(x = chr7_55808805_A_G, y = PSPHP1, color = chr7_55808805_A_G))+
geom_boxplot(alpha = 1, outlier.size = 0, size = 1, width = 0.5, fill = "transparent")+
geom_jitter(position = position_jitter(0.17), size = 1.5, alpha = 0.7)+
labs(x = "Genotype", y = "TPM value of PSPHP1", color = "rs34458430")+
theme(panel.background = element_blank(),
axis.text = element_text(hjust = 0.5, size = 18),
axis.title = element_text(hjust = 0.5, size = 18),
legend.title = element_text(size = 18, ),
legend.position = c(.05, .98),
legend.justification = c("left", "top"),
panel.border = element_rect(color = "black",fill = NA),
plot.title = element_text(hjust = 0.5, size = 26),
plot.subtitle = element_text(hjust = 0.5, size = 18))

#bar plot
#A.csv 、B.csv is the feature of the eQTL result

> a <- read.csv(file = "A.csv",header = TRUE, sep = ",")
> ce =ddply(a,"group", transform, percent_count = count /sum(count) * 100)
#order
> ce$type <- factor(ce$type, levels = c('eGene_pLI>0.9','eGene_pLI≤0.9','eGene_pLI=NA','eSNP_Gene_pLI>0.9',"eSNP_Gene_pLI≤0.9","eSNP_Gene_pLI=NA","Neither","GTEx_ColonSigmoid Only","GTEx_ColonTransverse Only","Both","eSNP_in_eGene","eSNP_near_eGene"))
#plot
p <- ggplot(ce,aes(x = group, y = percent_count, fill = type)) +
    geom_bar(stat ="identity", colour ="black")+
    scale_fill_brewer(palette = "Set3",direction = -1)+
    theme(panel.background = element_blank(),
          axis.text = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_text(size = 12),
          panel.border = element_rect(color = "black",fill = NA),
          plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 12))+
    labs(title = "Features of eQTLs", x="Items" , y = "Percentage", fill = "Types", size = 12)
#save
p


#another
> b <- read.csv(file = "B.csv",header = TRUE, sep = ",")
> ce =ddply(b,"group", transform, percent_count = count /sum(count) * 100)
> ce$type <- factor(ce$type, levels = c("CpG_island","NonCpG_island","ncRNA","exonic","intronic","UpstreamDownstream","splicing","UTR","intergenic","Active_promoter","Weak/poised_promoter","Strong_enhancer","Weak/poised_enhancer","Insulator","Polycomb_repressed","Transcriptional_transition","Transcriptional_elongation","Heterochromatin_LowSignal","Weak_transcribed","Others"))
p <- ggplot(ce,aes(x = group, y = percent_count, fill = type)) +
    geom_bar(stat ="identity", colour ="black")+
    theme(panel.background = element_blank(),
          axis.text = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_text(size = 12),
          panel.border = element_rect(color = "black",fill = NA),
          plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 12))+
    labs(title = "Features of eQTLs", x="Items" , y = "Percentage", fill = "Types", size = 12)
> p+ scale_fill_manual(values=c("#1f78b4","#a6cee3","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5"))

#plot for eQTL of the 38 controls
> c <- read.csv(file = "C.csv",header = TRUE, sep = ",")
> ce =ddply(c,"group", transform, percent_count = count /sum(count) * 100)
#order
> ce$type <- factor(ce$type, levels = c('eGene_pLI>0.9','eGene_pLI≤0.9','eGene_pLI=NA','eSNP_Gene_pLI>0.9',"eSNP_Gene_pLI≤0.9","eSNP_Gene_pLI=NA","Neither","GTEx_ColonSigmoid Only","GTEx_ColonTransverse Only","Both","eSNP_in_eGene","eSNP_near_eGene"))
#plot
p <- ggplot(ce,aes(x = group, y = percent_count, fill = type)) +
    geom_bar(stat ="identity", colour ="black")+
    scale_fill_brewer(palette = "Set3",direction = -1)+
    theme(panel.background = element_blank(),
          axis.text = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_text(size = 12),
          panel.border = element_rect(color = "black",fill = NA),
          plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 12))+
    labs(title = "Features of eQTLs", x="Items" , y = "Percentage", fill = "Types", size = 12)
#save
p


#the part2 of the eQTL of the 38 controls
d <- read.csv(file = "D.csv",header = TRUE, sep = ",")

ce =ddply(d,"group", transform, percent_count = count /sum(count) * 100)

ce$type <- factor(ce$type, levels = c("CpG_island","NonCpG_island","ncRNA","exonic","intronic","Upstream","UTR","intergenic","Active_promoter","Weak/poised_promoter","Weak_enhancer","Insulator","Polycomb_repressed","Transcriptional_transition","Transcriptional_elongation","Heterochromatin_LowSignal","Weak_transcribed","Others"))

p <- ggplot(ce,aes(x = group, y = percent_count, fill = type)) +
    geom_bar(stat ="identity", colour ="black")+
    theme(panel.background = element_blank(),
          axis.text = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_text(size = 12),
          panel.border = element_rect(color = "black",fill = NA),
          plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 12))+
    labs(title = "Features of eQTLs", x="Items" , y = "Percentage", fill = "Types", size = 12)

p+ scale_fill_manual(values=c("#1f78b4","#a6cee3","#b2df8a","#33a02c","#fb9a99","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5"))

#
