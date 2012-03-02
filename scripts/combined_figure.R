library(ggplot2)
library(gridExtra)

s<-read.table("stats.txt", sep="\t", header=TRUE)
s$assembly <- factor(s$assembly, levels=unique(as.character(s$assembly)))
p1 <- ggplot(s, aes(factor(assembly))) + geom_bar(position="dodge", aes(fill=assembler)) + facet_wrap( ~ gaptype, nrow = 2, ncol = 1 ) + opts(strip.text.x = theme_text(size = 10)) + opts(axis.text.x=theme_text(angle=90, hjust=1)) + xlab('Assembly') + ylab('Count')  + opts(legend.position="left")

hp<-s[s$gaptype == "Homopolymer gap",]
filter=c("454 GS Junior (1)", "454 GS Junior (2)", "454 GS Junior (1+2)", "Ion Torrent PGM (1)", "Ion Torrent PGM (2)", "Ion Torrent PGM (1+2)")
hp <- subset(hp, hp$assembly %in% filter)
hp$assembly<-factor(hp$assembly, levels=filter)
p2 <- ggplot(data=hp, aes(x=hplen, fill=assembler)) + geom_histogram(position="dodge", binwidth=.5) + facet_wrap( ~ assembly ) + opts(strip.text.x = theme_text(size = 9)) + xlab('Homopolymer length') + ylab('Count') + opts(legend.position="none")

# + stat_bin(breaks=seq(1,10, by=1)) 

pdf("out.pdf", width=12, height=8)
grid.arrange(p1, p2, ncol=2)
grid.text(label="a", x = 0.05, y = 0.96, gp = gpar(fontface="bold"))
grid.text(label="b", x = .55, y = 0.96, gp = gpar(fontface="bold"))
dev.off()
quit()
