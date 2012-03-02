#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
library(ggplot2)
a<-read.table(args[1], header=TRUE, sep="\t")
pdf(args[2])
boxplot(qual ~ base, data=a, outline = FALSE, ylim=c(0,40))
dev.off()

