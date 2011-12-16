#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
indels<-read.table(args[1], sep="\t", header=TRUE)
mapped <- indels[indels$mapped == 1,]
sum(indels$insertions + indels$deletions) / sum(indels$rlen) * 100
sum(indels$insertions + indels$deletions) / nrow(indels)
