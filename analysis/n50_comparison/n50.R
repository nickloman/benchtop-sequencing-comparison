#!/usr/bin/Rscript

library(reshape)
library(ggplot2)

assemblies<-read.table("n50_data.txt", sep="\t")
a<-melt(assemblies, id = c("V1", "V2"), measured = c( "V3"))
a$V1 <- factor(a$V1, levels = c("454 GS Junior (1)", "454 GS Junior (2)", "454 GS Junior (1+2)", "Ion Torrent PGM (1)", "Ion Torrent PGM (2)", "Ion Torrent PGM (1+2)", "MiSeq (contigs)", "MiSeq (scaffolds)"))
pdf("assembly_comparison_n50.pdf")
ggplot(a, aes(x = V1, y = value, shape = factor(a$V2))) + geom_point(aes(colour = factor(a$V2)), size=4) +   geom_point(colour="grey90", size = 1.5) + scale_colour_discrete("Assembler") + scale_shape_discrete("Assembler") + scale_x_discrete("Sequencing Run") + scale_y_continuous("Assembly N50 (bases)", limits=c(0, max(a$value))) + opts(axis.text.x=theme_text(angle=45, hjust=1.0, vjust=1.0))
dev.off()

