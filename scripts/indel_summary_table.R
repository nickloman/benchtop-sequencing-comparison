#!/usr/bin/Rscript

#usage: <input from read_bam.py> <output file>

library(xtable)
args<-commandArgs(TRUE)
indels<-read.table(args[1], sep="\t", header=TRUE)
mapped<-indels[indels$mapped == 1,]
tab <- cbind( "insertions" = tapply(mapped$insertions, mapped$sample, sum),
        "deletions" = tapply(mapped$deletions, mapped$sample, sum),
        "bases" = tapply(mapped$rlen, mapped$sample, sum),
        "reads" = tapply(mapped$rid, mapped$sample, length )  )
fram<-as.data.frame(tab)
summary<-cbind(fram,
          "indels_per_100" = (fram$insertions + fram$deletions) / fram$bases * 100,  "indels_per_read" = (fram$insertions + fram$deletions) / fram$reads )
#summary2<-cbind(summary, "total_reads" = tapply(indels$rlen, indels$sample, length))
#summary3<-cbind(summary2, "% mapped reads" = summary2$total_reads / summary2$reads)
print(xtable(summary[,c(1,2,5,6)]), file = args[2])

