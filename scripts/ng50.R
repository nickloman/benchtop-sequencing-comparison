# Thanks to Lex Nederbragt for script

pdf("ng50.pdf")
al<-read.table("lengths.txt",header=TRUE,sep="\t")
 	
# ass<-unique(al$description)
ass <- unique(factor(al$description, levels=unique(as.character(al$descriptio))))


a1 <- al[al$description == ass[1],]

gxlim<-max(al$length)

ca1 = cumsum(sort(a1$length, decreasing=T))/sum(a1$length)
xlab="Minimum contig/scaffold length"
ylab="Percentage assembly covered"
plot(sort(a1$length, decreasing=T), 100*ca1, xlab=xlab, ylab=ylab, log="x", xlim = c(100, 1000000), xaxt = "n",type="l",col=1)
axis(1,at=c(1,10,100,1000,10000,100000,1000000),labels=c("1","10","100","1,000","10,000","100,000","1,000,000"))
abline(v=c(10000,100000,1000000),lty=6,col=gray(0.4))

i<-2
for (a in ass[2:length(ass)]) {
  d <- sort(al[al$description == a,]$length, decreasing=T)
  ca1 = cumsum(d) / sum(d)
  lines(d, ca1*100, type="l", col=i)
  i<- i + 1
}

legend("bottomleft", legend = ass, col=c(1:length(ass)), lty=1)
dev.off()
