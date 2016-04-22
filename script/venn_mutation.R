library(VennDiagram)

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
title <- args[3]
fish_pval <- args[4]
ebcall_pval <- args[5]
realign_pval <- args[6]
variant_t <- args[7]
variant_n <- args[8]
genomon_total <- args[9]
firehose_total <- args[10]

x <- read.table(infile)
firehose_tmp<-x[,1]
genomon_tmp<-x[,2]

firehose = firehose_tmp[firehose_tmp > 0]
genomon = genomon_tmp[genomon_tmp > 0]

newtitle <- paste(title, "\nP-value(fisher) >=", fish_pval, "P-value(EBCall) >=", ebcall_pval, "P-value(fhsher_realignment) >=", realign_pval, "\nvariantPairNum_tumor >=", variant_t, "variantPairNum_normal <=", variant_n, sep = " ")
# venn.diagram( list(Firehose=firehose,Genomon=genomon), filename=outfile, ext.length = 0, ext.dist = -0.15, cex = 1.5, cat.cex = 1.5, main = newtitle, fill=c("dodgerblue", "goldenrod1"))
# venn.diagram( list(Firehose=firehose,Genomon=genomon), filename=outfile, ext.length = 0, ext.dist = -0.15, cex = 1.5, cat.cex = 1.5,  height = 4000, width = 4000, main = newtitle, fill=c("dodgerblue", "goldenrod1"))
if (firehose_total == 0) {
    venn.diagram( list(Genomon=genomon), filename=outfile, ext.length = 0, ext.dist = -0.15, cex = 1.5, cat.cex = 1.5,  height = 4000, width = 4000, main = newtitle, fill=c("goldenrod1"))
} else if (genomon_total == 0) {
    venn.diagram( list(Firehose=firehose), filename=outfile, ext.length = 0, ext.dist = -0.15, cex = 1.5, cat.cex = 1.5,  height = 4000, width = 4000, main = newtitle, fill=c("dodgerblue"))
} else {
    venn.diagram( list(Firehose=firehose,Genomon=genomon), filename=outfile, ext.length = 0, ext.dist = -0.15, cex = 1.5, cat.cex = 1.5,  height = 4000, width = 4000, main = newtitle, fill=c("dodgerblue", "goldenrod1"))
}

