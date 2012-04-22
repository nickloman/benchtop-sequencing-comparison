#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)
library(reshape)
library(qrqc)

observed_quality_plot <- function(fn, title, lm) {
	a<-read.table(fn, header=TRUE)
	accuracy <- expected_accuracy(a)
	b<-cast(a, pos ~ variable, sum)
	total_mismatches<-b$subst_mismatches + b$indel_mismatches
	calls<-total_mismatches + b$matches
	qual<- 100 - (100 * (total_mismatches / calls))
	b<-cbind(b, qual)
	b<-b[b$matches > 1000,]
	accuracy<-accuracy[accuracy$matches > 1000,]
	print("plotting")
	# geom_smooth
	p<-ggplot(b, aes(x=pos, y=qual)) +
		geom_line(aes(linetype="Measured")) +
		geom_line(data = accuracy, aes(x=pos, y=eaccuracy, linetype="Expected")) +
		scale_x_continuous("Base Position", limits=c(0,lm)) +
		scale_y_continuous("Accuracy (%)", limits=c(95,100)) +
		scale_linetype_discrete("") +
		opts(title = title) + opts(legend.position="top")
	print("doneplot")
	p
}

aligned_plot <- function(fn, title, ylim) {
	a<-read.table(fn, header=TRUE)
	b<-cast(a, pos  ~ variable, sum)
	b<-cbind(b, "aligned" = 100 * (1 - b$unaligned / (b$indel_mismatches + b$matches + b$subst_mismatches + b$unaligned)))
	ggplot(b, aes(x=pos, y=aligned)) + geom_line() + scale_y_continuous("Reads aligned (%)", limits=c(0,100)) + scale_x_continuous("Base position", limits=c(0,ylim)) + opts(title = title)
}

expected_accuracy <- function(a) {
	b<-cast(a, pos + qual ~ variable, sum)
	b<-cbind(b, "totals" = b$matches + b$indel_mismatches + b$subst_mismatches)
	b<-cbind(b, "expected" = b$totals * 10 ^ ( -0.1 * b$qual ))
	c<-melt(b, id=c("pos", "qual"))
	d<-cast(c, pos ~ variable, sum)
	cbind(d, "eaccuracy" = 100 - (100 * (d$expected / d$totals)))
}

formatBack <- function(x) sub("[0]+$", "", format(100 - ((10^(-x/10)*100) ), nsmall = 3))

read_fastqc <- function(fn) {
	data=readLines(fn)
	ln<-grep("Per base sequence quality", data)
	lnend<-grep(">>END_MODULE", data)
	read.table(textConnection(sub('#', '', data[14:lnend[2]-1])), header=TRUE, sep="\t")
}


#../fastqc/MiSeq-1_280/MiSeq-1_280_fastqc/fastqc_data.txt
#../fastqc/IonTorrent-both_fastqc/fastqc_data.txt
#../fastqc/454Junior-both_fastqc/fastqc_data.txt

a<-read_fastqc("../fastqc/454Junior-both_fastqc/fastqc_data.txt")
b<-read_fastqc("../fastqc/IonTorrent-both_fastqc/fastqc_data.txt")
c<-read_fastqc("../fastqc/MiSeq-1_280/MiSeq-1_280_fastqc/fastqc_data.txt")

oa<-observed_quality_plot("454.accuracy.txt", "454 GS Junior 1+2", 600)
ob<-observed_quality_plot("iontor.accuracy.txt", "Ion Torrent PGM 1+2", 200)
od<-observed_quality_plot("../../tmap_alignments/iontor.tmap.map1.accuracy.txt", "Ion Torrent PGM - tmap1", 200)
oc<-observed_quality_plot("miseq.accuracy.txt", "MiSeq (strain 280)", 200)

aa<-aligned_plot("454.accuracy.txt", "454 GS Junior 1+2", 600)
ab<-aligned_plot("iontor.accuracy.txt", "Ion Torrent PGM 1+2", 200)
ad<-aligned_plot("../../tmap_alignments/iontor.tmap.map1.accuracy.txt", "Ion Torrent PGM - tmap1", 200)
ac<-aligned_plot("miseq.accuracy.txt", "MiSeq (strain 280)", 200)

doplot_Fastqc <- function(pred, obs, l, title) {
	ggplot(data=pred) +
# geom_boxplot(data = pred, aes(base, qual, group=round_any(base, 10), outlier.colour="NA")) +
		geom_boxplot(aes(x=Base, ymin=X10th.Percentile, lower=Lower.Quartile, middle=Median, upper=Upper.Quartile, ymax=X90th.Percentile), stat="identity") +
		scale_x_continuous(limits=c(0,l)) +
        	scale_y_continuous(limits=c(0,40)) +
        	geom_line(data=pred, aes(x=Base, y=Mean))
		opts(title=title)
}

doplot_qrqc <- function(fn, qual, ylim, maxlength) {
	fq<-readSeqFile(fn, quality=qual, hash=FALSE, max.length=maxlength)
	qualPlot(fq) + scale_y_continuous("Base quality (phred-scaled)", limits = c(0,40)) +
	scale_x_continuous("Base position", limits=c(0,ylim))
}

seqlenplot_qrqc <- function(fn, qual, ylim, maxlength) {
	fq<-readSeqFile(fn, quality=qual, hash=FALSE, max.length=maxlength)
	sl <- getSeqlen(fq)
	sl$density <- sl$count/sum(sl$count)
	# ggplot(sl) + geom_bar(aes(x=length, y=density), stat="identity") + scale_y_continuous("Density") + scale_x_continuous("Read length", limits=c(0,ylim))
	seqlenPlot(fq) + scale_y_continuous("Number") + scale_x_continuous("Read length", limits=c(0,ylim))
}

pdf("fig1_panelA.pdf", width=20, height=26)
#p1<-doplot(a, oa, 500, "454 GS Junior 1+2")
#p2<-doplot(b, ob, 150, "Ion Torrent PGM 1+2")
#p4<-doplot(b, od, 150, "Ion Torrent PGM 1+2")
#p3<-doplot(c, oc, 150, "MiSeq (strain 280)")
p1<-doplot_qrqc("../fastqc/454Junior-both.fastq", "sanger", 600, 2000) + opts(title="454 GS Junior 1+2")
p2<-doplot_qrqc("../fastqc/IonTorrent-both.fastq", "sanger", 200, 400) + opts(title="Ion Torrent PGM 1+2")
p3<-doplot_qrqc("../../reads/MiSeq-1_280.fastq", "illumina", 200, 400) + opts(title="MiSeq (280 strain)")

lp1<-seqlenplot_qrqc("../fastqc/454Junior-both.fastq", "sanger", 600, 2000) + opts(title="454 GS Junior 1+2")
lp2<-seqlenplot_qrqc("../fastqc/IonTorrent-both.fastq", "sanger", 200, 400) + opts(title="Ion Torrent PGM 1+2")
lp3<-seqlenplot_qrqc("../../reads/MiSeq-1_280.fastq", "illumina", 200, 400) + opts(title="MiSeq (280 strain)")
grid.arrange(p1, p2, p3, lp1, lp2, lp3, oa, ob, oc, aa, ab, ac, ncol=3, nrow=4)

#
#p2<-ggplot(b, aes(base, qual, group=round_any(base, 10))) +
#	geom_boxplot(outlier.colour="NA") +
#	scale_y_continuous(limits=c(0,40)) +
#	scale_x_continuous(limits=c(0,150)) +
#	stat_summary(fun.y=mean, geom="line", aes(group=1))
#p3<-ggplot(c, aes(base, qual, group=round_any(base, 10))) +
#	geom_boxplot(outlier.colour="NA") +
#	scale_y_continuous(limits=c(0,40)) +
#	scale_x_continuous(limits=c(0,150)) +
#	stat_summary(fun.y=mean, geom="line", aes(group=1))

#grid.arrange(p1, p2, p3, ncol=3, nrow=1)
dev.off()

pdf("fig1_panelB.pdf")

tab1<-read.table("454.accuracy.txt", header=TRUE)
tab2<-read.table("iontor.accuracy.txt", header=TRUE)
tab3<-read.table("miseq.accuracy.txt", header=TRUE)
tab3$qual<-tab3$qual - 31
tab1$"Platform" = "454 GS Junior 1+2"
tab2$"Platform" = "Ion Torrent PGM 1+2"
tab3$"Platform" = "MiSeq (280)"
a<-rbind(tab1, tab2, tab3)
b<-a
b<-cast(b, Platform + qual ~ variable, sum)
b$totals<-b$indel_mismatches + b$matches + b$subst_mismatches
b<-b[b$qual > 0,]
b<-cbind(b, "accuracy" = -10 * log10((b$indel_mismatches + b$subst_mismatches) / b$totals))

p7<-ggplot(b, aes(x = qual, y = accuracy, colour = Platform)) +
	geom_point(aes(size=totals)) + scale_area() + scale_size("Number of base calls") +
	#, breaks=c(4.5, 5, 5.5, 6, 6.5, 7.0, 7.5, 8, 8.5)) +
	geom_line() +
	geom_abline(intercept=0, slope=1) +
	#, breaks=c(4,5,6,7,8), labels=c("<1000", "<10000", "<100000", "<1000000", "<10000000")) +
	scale_x_continuous("Bin of predicted quality score (phred-scale)", limits=c(0,40)) +
	scale_y_continuous("Accuracy (%)", limits=c(0, 40), formatter='formatBack')

p7

dev.off()

