#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
indels<-read.table(args[1], sep="\t", header=TRUE)
mapped <- indels[indels$mapped == 1,]
# sum(indels$insertions + indels$deletions) / sum(indels$rlen) * 100
# sum(indels$insertions + indels$deletions) / nrow(indels)
indels_per_cent <- sum(mapped$insertions + mapped$deletions) / sum(mapped$rlen) * 100
indels_per_read <- sum(mapped$insertions + mapped$deletions) / nrow(mapped)
print(c(args[1], indels_per_cent, indels_per_read))
