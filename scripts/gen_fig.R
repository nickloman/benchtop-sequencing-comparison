#!/usr/bin/Rscript --no-save --no-restore  
## IONTorrent.R
## Colin Hercus
## SYNOPSIS: Plot Run Length error charts
## Usage: IONTorrent.R  -b <base> -l <dL> <reportfile.pdf> < indelstats.tsv
print("Initialising")
library('getopt')
tile=100
fontsize <- 10
outfile="indelreport.pdf"
infile="indels.tsv"
fit=FALSE

opt = getopt(c(
	'help' , 'h', 0, "logical",
	'tilesize', 't', 1, "integer",
	'file', 'f', 1, "character",
	'gofit', 'g', 0, "logical",
    'report', 'r', 1, "character"))
#help was asked for.
if ( !is.null(opt$help)) {
  self = commandArgs()[1];
  cat("Usage: ",self," [options]\n",
        " [-h]                      Print thiss usage help information",
		" [-f <indelstatsfile.tsv>]	Indel hpstats file from Novoalign, default ", infile, "\n",
		" [-r <reportfile.pdf>]	    Output report file, default ", outfile, "\n",
		" [-g]                      Report regression coefficients for IonTorrent Indel model.\n" ,
		" [-t <tilesize>]           XY stats reported for tiles of this size, default", tile, "\n")
  q(status=-1);
}

if(!is.null(opt$tile)) { tile = opt$tile };
if(!is.null(opt$report)) {outfile = opt$report};
if(!is.null(opt$file)) {infile = opt$file};
fit = !is.null(opt$gofit)


maxq=50
library(gridExtra)
library(reshape)
library(ggplot2)
library(plyr)


plotter <- function(infile, platform) {

print("Read Data")
data = readLines(infile)
command = ""
print("Mark")
mark <- vector('integer', length(data))
starthp <- grep("^#HPCOUNTS", data)
if(starthp > 1) { command = sub("# (.*)","\\1",data[1]) };
startxy <- grep("^#XYCOUNTS", data)
startrl <- grep("^#READ LENGTH", data)
startid <- grep("^#INDELS", data)
mark[starthp] = 1
mark[startxy] = 1
mark[startrl] = 1
mark[startid] = 1
mark <- cumsum(mark)
mark[starthp] = 0
mark[startxy] = 0
mark[startrl] = 0
mark[startid] = 0

print("Extract")
hpcounts=read.csv(textConnection(data[mark==1]), sep='\t')
xy=read.csv(textConnection(data[mark==2]), sep='\t')
rlen=read.csv(textConnection(data[mark==3]), sep='\t')
indel=read.csv(textConnection(data[mark==4]), sep='')

indel <- indel[indel$Length < 9,]

print("RLen")
TotalReads = sum(rlen$Count)
maxrlen = 100

#max(rlen[rlen$Count > 0,]$Length)
meanrlen=weighted.mean(rlen$Length,rlen$Count)
rlen$accum = rlen$Count
rlen$cover = cumsum(rlen$Start) - cumsum(rlen$End)
for(i in 1:nrow(rlen))
	rlen[i,]$accum =sum(rlen[i:nrow(rlen),]$Count)
rlen = rlen[rlen$accum > 0,]


UsedReads=sum(xy$Reads)

#plot
print("Plots over Run Length")

hpcountsT=hpcounts[hpcounts$Q.T==1,]
hpcountsQ=hpcounts[hpcounts$Q.T==0,]
hpcountsT$Q.T=NULL
hpcountsQ$Q.T=NULL
#hpcounts$nbp = hpcounts$Length*hpcounts$Count
#hpcounts[1,]$Count = sum(hpcounts$nbp)

adj = 1.0 * TotalReads / UsedReads

indel$adjshort = 1.0 * indel$D * rlen[1,]$accum / rlen[min(c(max(rlen$Length)+1,indel$Offset+indel[2,]$Offset/2)),]$cover
indel$adjlong = 1.0 * indel$I * rlen[1,]$accum / rlen[min(c(max(rlen$Length)+1,indel$Offset+indel[2,]$Offset/2)),]$cover
t = count(indel$Offset)
nlen = nrow(rlen[rlen$accum > 0.3*rlen[1,]$accum,])/indel[2,]$Offset
indel$pshort = (-10*log10(adj*indel$adjshort/(hpcountsT[indel$Length,]$Count/nlen)))
indel$plong = (-10*log10(adj*indel$adjlong/(hpcountsQ[indel$Length,]$Count/nlen)))
indel$shortaccuracy = (adj*indel$adjshort/(hpcountsT[indel$Length,]$Count/nlen)) 
indel$longaccuracy = (adj*indel$adjlong/(hpcountsQ[indel$Length,]$Count/nlen))
indel$Lt = factor(indel$Length)
levels(indel$Lt)[length(levels(indel$Lt))]=paste(">=",levels(indel$Lt)[length(levels(indel$Lt))])
indel$Lq = indel$Lt
sig.limit=5

# maxq=max(c(max(indel[indel$D >sig.limit,]$pshort),max(indel[indel$I >sig.limit,]$plong),40))$
maxq = 50

print("Regression")

indel$P="Measured"
if(fit) {
	indel$logL = log(indel$Length)
	indel$OL = indel$Offset * indel$logL
	Olimit =max(rlen[rlen$accum > rlen[1,"accum"]*.67,"Length"]) - 25
	Olimit = max(c(Olimit, 30))

	lrlong = nls(plong ~ A + B*Offset + C*logL + D*OL, indel[indel$Offset >= 15 & indel$Offset <= Olimit & indel$I > sig.limit & indel$Length < 9 & indel$Length > 1,],
	 	start=list(A=0,B=0,C=0,D=.0))
	lrinsert = nls(plong ~ A + B*Offset, indel[indel$Offset >= 15 & indel$Offset <= Olimit & indel$I > sig.limit & indel$Length == 1,],
	 	start=list(A=0,B=.0))
	lrshort = nls(pshort ~ A + B*Offset + C*logL + D*OL, indel[indel$Offset >= 15 & indel$Offset <= Olimit & indel$D > sig.limit & indel$Length < 9,],
	 	start=list(A=0,B=0,C=0,D=0))
	indelP = indel[indel$Length < 9,]
	indelP$plong = predict(lrlong,indelP[,c("Offset","logL","OL")])
	indelP[indelP$Length == 1,]$plong = predict(lrinsert,indelP[indelP$Length == 1,c("Offset","Length")])
	indelP$pshort = predict(lrshort,indelP[,c("Offset","logL","OL")])
	indelP$P="Regression"
	combined =rbind( indel , indelP)
} else {
	combined = indel
}
combined$group = interaction(combined$P, combined$Lt)
if(fit) { 
	xl=signif(coef(lrlong),2)
	xi=signif(coef(lrinsert),2)
	xs=signif(coef(lrshort),2)
	cs=paste(xs["A"],",",xs["B"],",",xs["C"],",",xs["D"],sep="")
	cl=paste(xl["A"],",",xl["B"],",",xl["C"],",",xl["D"],",",xi["A"],",",xi["B"],sep="")
	cat("-g ",cs,",",cl,"\n",sep="") 
} else {
cs = ""
cl = ""
}

sig.limit=3

plot.short  <- ggplot(combined[combined$D > sig.limit,], aes(x=Offset, y=shortaccuracy, group=group, color=Lt)) +
	scale_y_continuous("Frequency of deletion errors", limits=c(0, 1)) + 
	scale_x_continuous("Read position", limits=c(0, maxrlen)) +
	scale_colour(name = "Homopolymer\nLength") +
	opts(title = paste(platform, "(deletions)")) +
	geom_line()
	
# plot.short <- plot.short + opts(legend.position = "bottom", legend.direction = "horizontal") 

#if (legend) {
#	plot.short <- plot.short + scale_colour(name = "Homopolymer\nLength")
#}

plot.long <- ggplot(combined[combined$I > sig.limit,], aes(x=Offset, y=longaccuracy, group=group, color=Lq)) +
	geom_line() + 
	opts(title = paste(platform, "(insertions)")) +
	scale_colour(name = "Homopolymer\nLength") +
	scale_y_continuous("Frequency of insertion error", limits=c(0, 1)) + 
	scale_x_continuous("Read position", limits=c(0, maxrlen))

# plot.long <- plot.long + opts(legend.position = "none") 

print("XY Plots")

SZ=as.integer(120/max(xy$Y/tile+2))

xy.tile = xy
xy.tile$Y = tile*floor(as.integer(xy.tile$Y)/tile)
xy.tile$X = tile*floor(as.integer(xy.tile$X)/tile)
xy.tile=ddply(xy.tile,c("X","Y"),summarise, Reads=sum(Reads), D =sum(D), D0=sum(D0), I0=sum(I0), I=sum(I))
plot.xyreads <- ggplot(xy.tile, aes(x=X, y=Y, size=Reads, color=Reads)) +
	geom_point() + 
	opts(title = "Number of Aligned Reads") +
	scale_colour_gradient(low="red", high="green",name = "Reads") +
	scale_size(to = c(1, SZ),name = "Reads")


xy.tile$short <- xy.tile$D0 + xy.tile$D

xy.tile$Ratio=-10*log10(xy.tile$short/xy.tile$Reads/meanrlen)
limits=c(50,10)
breaks = 5*(ceiling(limits[2]/5.0):floor(limits[1]/5.0))
labels=breaks
plot.xyshortvsreads<- ggplot(xy.tile[xy.tile$short > 2,], aes(x=X, y=Y, size=Ratio, color=Ratio)) +
	geom_point() + 
	opts(title = "Probability of Deletion Event, Phred Scale") +
	scale_colour_gradient(low="green", high="red",name = "-log(P)", breaks=breaks, labels =labels, limits=limits) +
	scale_size(to = c(1, SZ),name = "-log(P)", breaks=breaks, labels =labels, limits=limits)

xy.tile$Ratio=-10*log10(xy.tile$I/xy.tile$Reads/meanrlen)
plot.xylongsvsreads <- ggplot(xy.tile[xy.tile$I > 2,], aes(x=X, y=Y, size=Ratio, colour=Ratio)) +
	geom_point() + 
	opts(title = "Probability of Insertion Event, Phred Scale\n(That matches an adjacent base)") +
	scale_colour_gradient(low="green", high="red",name = "-log(P)", breaks=breaks, labels =labels, limits=limits) +
	scale_size(to = c(1, SZ),name = "-log(P)", breaks=breaks, labels =labels, limits=limits)

xy.tile$Ratio=-10*log10(xy.tile$I0/xy.tile$Reads/meanrlen)
plot.xyinsertsvsreads <- ggplot(xy.tile[xy.tile$I0 > 2,], aes(x=X, y=Y, color=Ratio, size=Ratio)) +
	geom_point() + 
	opts(title = "Probability of Insertion Event, Phred Scale\n(That does not match an adjacent base)") +
	scale_colour_gradient(low="green", high="red",name = "-log(P)", breaks=breaks, labels =labels, limits=limits) +
	scale_size(to = c(1, SZ),name = "-log(P)", breaks=breaks, labels =labels, limits=limits)

return( list(short=plot.short, long=plot.long) )

}

ion <- plotter("BHAM5_6_1.5.indel.csv", "Ion Torrent")
four <- plotter("HPA2_HPA4.indel.csv", "454 Junior")

print(paste("Writing output to",outfile))
#Open PDF output
pdf(outfile,paper="a4",pointsize=12, width=12, height=8)
leg<- ggplotGrob(ion$short + opts(keep="legend_box"))
dev.off()
pdf(outfile,paper="a4",pointsize=12, width=12, height=8)
# The theme

old_theme <- theme_update(axis.text.x=theme_text(size=fontsize),axis.text.y=theme_text(size=fontsize),
                          axis.title.x=theme_text(size=fontsize),axis.title.y=theme_text(size=fontsize,angle=90),
                          strip.text.x=theme_text(size=fontsize),strip.text.y=theme_text(size=fontsize),
                          plot.title=theme_text(size=fontsize*1.2),
                          legend.text=theme_text(size=fontsize),legend.title=theme_text(size=fontsize),
			  legend.key.size = unit(0.4, "lines"), legend.key = theme_rect(colour = "transparent")) 

legend = gTree(children=gList(leg), cl="legendGrob")
widthDetails.legendGrob<- function(x) unit (4, "cm")

#Arrange the plots skipped plot.shortvslong,
grid.arrange( 
	four$short + opts(legend.position="none"),
	four$long + opts(legend.position="none"),
	ion$short + opts(legend.position="none"),
	ion$long + opts(legend.position="none"),
	legend = legend,
	ncol=2, nrow=2)
	#, gp=gpar(font=2,fontsize=fontsize)))
dev.off()

