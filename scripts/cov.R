library(GenomeGraphs)
library(GenomicFeatures)
library(ShortRead)
library(HilbertVis)
aln<-readGappedAlignments("bwa/ILMN1_L5_280_filtered.sorted.bam")
ionaln<-readGappedAlignments("bwa/BHAM5_1.5.sorted.bam")
ion2aln<-readGappedAlignments("bwa/BHAM6_1.5.sorted.bam")
gsaln<-readGappedAlignments("bwa/HPA2.sorted.bam")
gs2aln<-readGappedAlignments("bwa/HPA4.sorted.bam")

cov<-coverage(aln)
ioncov<-coverage(ionaln)
ion2cov<-coverage(ion2aln)
gscov<-coverage(gsaln)
gs2cov<-coverage(gs2aln)

covnum<-as.numeric(cov$scaffold00001)
ioncovnum<-as.numeric(ioncov$scaffold00001)
ion2covnum<-as.numeric(ion2cov$scaffold00001)
gscovnum<-as.numeric(gscov$scaffold00001)
gs2covnum<-as.numeric(gs2cov$scaffold00001)

l<-list(makeGenomeAxis(T, F, T),
   "MiSeq (strain 280)" =makeBaseTrack(
       base = seq(0, length(covnum), length.out=10000),
       value=shrinkVector( as.vector(covnum), 10000 ),
       dp = DisplayPars(main = "Illumina MiSeq (strain 280)", ylim = c(0:100), pch = 20, cex=0.2, lwd = 0.2, color = 1)
   ),
   "Ion Torrent 1"=makeBaseTrack(
       base = seq(0, length(ioncovnum), length.out=10000),
       value=shrinkVector( as.vector(ioncovnum), 10000 ),
       dp = DisplayPars(main = "Ion Torrent PGM Run 1", ylim = c(0:100), pch = 20, cex=0.2, lwd = 0.2, color = 2 )
   ),
   "Ion Torrent 2"=makeBaseTrack(
       base = seq(0, length(ion2covnum), length.out=10000),
       value=shrinkVector( as.vector(ion2covnum), 10000 ),
       dp = DisplayPars(main = "Ion Torrent PGM Run 2", ylim = c(0:100), pch = 20, cex=0.2, lwd = 0.2, color = 3 )
   ),
   "454 Junior 1"=makeBaseTrack(
       base = seq(0, length(gscovnum), length.out=10000),
       value=shrinkVector( as.vector(gscovnum), 10000 ),
       dp = DisplayPars(main = "454 GS Junior Run 1", lim = c(0:100), pch = 20, cex=0.2, lwd = 0.2, color = 4 )
   ),
   "454 Junior 2"=makeBaseTrack(
       base = seq(0, length(gs2covnum), length.out=10000),
       value=shrinkVector( as.vector(gs2covnum), 10000 ),
       dp = DisplayPars(main = "454 GS Junior Run 2", ylim = c(0:100), pch = 20, cex=0.2, lwd = 0.2, color = 4 )
   ),
   makeGenomeAxis(T, F, T)
)

args<-commandArgs(TRUE)
pdf(args)
gdPlot(l, DisplayPars(cex.lab = 0.5))
dev.off()

