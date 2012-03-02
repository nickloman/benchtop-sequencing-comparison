#!/usr/bin/env Rscript

args<- commandArgs(TRUE)
library(reshape)
a<-read.table(args[1], header=TRUE)
b<-cast(a, pos ~ variable, sum)

total_mismatches<-b$subst_mismatches + b$indel_mismatches
calls<-total_mismatches + b$matches + b$unaligned

- 10 * log10(total_mismatches / calls)

-10 * log10(b$subst_mismatches / calls)

100 - (b$subst_mismatches / (b$matches+b$subst_mismatches+b$unaligned) * 100)
#c<- cast(a, qual ~ variable, mean)
#c<-cbind(c, "subst_score" = -10 * log10(c$subst_mismatches / (c$matches+c$subst_mismatches+c$unaligned)))
#c
