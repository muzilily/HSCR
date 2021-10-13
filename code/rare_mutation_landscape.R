library('GenVisR')
library(reshape2)

#no phenolyzer gene (222genes)
#read data
mutation_data=read.table(file="sam150.snp.indel.Gene_Scan.filter.HGMD.MGI.xls",header=TRUE, sep="\t")
#sample, genename, trv_type, chrom, pos, aachange
mutation=mutation_data[,c(1,2,3,4,5,6,17)]
#sample order
#add sample's clinical info
clinical=read.table(file="sam33.clinical.txt", header=TRUE, sep="\t")
clinicaldata <- melt(clinical, id.vars = c("sample"))
sampleorder = as.vector(clinical$sample)
#gene order
gene=read.table(file="gene16order",header=TRUE,sep="\t")
geneorder = as.vector(gene$genename)

pdf("gene222.snp.indel.glist.HGMD.landscape.addGender.pdf",width = 14,height = 9)
waterfall(mutation, fileType="MGI",clinDat = clinicaldata, clinLegCol = 2,
          mainXlabel=TRUE,
          sampOrder = sampleorder,
          geneOrder = geneorder, mainDropMut = TRUE)
dev.off()
