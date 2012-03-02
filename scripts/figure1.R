pdf("out.pdf", width=6, height=8)

plot_empirical_vs_stated_quals <- function(datafile, title) {
# code from GATK

        Qcutoff = 0
        maxQ = 45

        t=read.table(datafile, header=T)
        d.good <- t[t$nBases >= 10000 & t$Qreported >= Qcutoff,]
        d.1000 <- t[t$nBases < 1000  & t$Qreported >= Qcutoff,]
        d.10000 <- t[t$nBases < 10000 & t$nBases >= 1000  & t$Qreported >= Qcutoff,]
        f <- t[t$Qreported < Qcutoff,]
        e <- rbind(d.good, d.1000, d.10000)
        rmseGood = sqrt( sum(as.numeric((d.good$Qempirical-d.good$Qreported)^2 * d.good$nBases)) / sum(as.numeric(d.good$nBases)) ) # prevent integer overflow with as.numeric, ugh
        rmseAll = sqrt( sum(as.numeric((e$Qempirical-e$Qreported)^2 * e$nBases)) / sum(as.numeric(e$nBases)) )
        theTitle = paste("RMSE_good =", round(rmseGood,digits=3), ", RMSE_all =", round(rmseAll,digits=3))
        if( length(t$nBases) - length(f$nBases) == length(d.good$nBases) ) {
                theTitle = paste("RMSE =", round(rmseAll,digits=3));
        }

        theTitle = title
        plot(d.good$Qreported, d.good$Qempirical, type="p", col="blue", main=theTitle, xlim=c(0,maxQ), ylim=c(0,maxQ), pch=16, xlab="Reported Q score", ylab="Observed Q score", mar=c(0,2), oma=c(0,0,0,0))
        points(d.1000$Qreported, d.1000$Qempirical, type="p", col="lightblue", pch=16)
        points(d.10000$Qreported, d.10000$Qempirical, type="p", col="cornflowerblue", pch=16)
        points(f$Qreported, f$Qempirical, type="p", col="maroon1", pch=16)
        abline(0,1, lty=2)
}

read_length_plot <- function() {
        files = c("454_hpa2_counts.txt",
                "454_hpa4_counts.txt",
                "illumina_280_counts.txt",
                "illumina_all_counts.txt",
                "iontorrent_bham5_count.txt",
                "iontorrent_bham6_count.txt")
        length_labels = c("454 Junior (1)", "454 Junior (2)", "MiSeq (strain 280)", "MiSeq (all reads)", "Ion Torrent (1)", "Ion Torrent (2)")
        colour_list <- c("red", "red", "blue", "blue", "purple", "purple")
        line_styles <- c(1, 2, 1, 2, 1, 2)

        tables <- list()
        for(f in files) {
                t <- read.table(f, header=TRUE)
                brks <- cut(t$length, cutlevels, labels=la, include.lowest=TRUE)
                tables[[f]] <- tapply(t$count, brks, mean)
        }
        i <- 2
        thick <- 2
        plot(tables[[ files[1] ]],
             log="y",
             axes = FALSE,
             ylab="Number of reads",
             xlab="Read position",
             ylim=c(1, 1e09),
             type="l",
             lwd=thick,
             col = colour_list [[ 1 ]],
              lty = line_styles[i - 1],
             mar=c(0,2), oma=c(0,0,0,0))
#            main="a) Number of reads at least a certain read length")
        for ( f in files[2:length(files)] ) {
                points(tables[[ f ]], type="l", col=colour_list[[ i ]], lwd = thick, lty = line_styles[ i ])
                i<-i + 1
        }
        legend( "bottomleft", length_labels, col = colour_list, lwd = thick, lty=line_styles )
        axis(1, 1:nlevels(brks), levels(brks))
        axis(2)
	par(xpd=FALSE)
	text(38, 92562135, "a", font=2)
}

layout(matrix(c(1,1,1,2,2,2,3,4,5), 3, 3, byrow=TRUE))
par(mar=c(4.5,4,1,2))

cutlevels <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700)

la <- c( "1",     "2",     "3",     "4",     "5",     "6",
 "7",     "8",     "9-15",    "16-20",   "21-25",   "26-30",
 "31-35",   "36-40",   "41-45",   "46-50",   "51-55",   "56-60",
 "61-65",   "66-70",   "71-80",   "81-90",   "91-100",  "101-110",
"111-120", "121-130", "131-140", "141-150",
 "151-200", "201-250", "251-300", "301-350", "351-400", "401-450",
 "451-500", "501-550", "551-600" ,"601-650", "651-700" )
read_length_plot()

ilmn <- read.table("ilmn_quals.txt")
iont <- read.table("iontorrent.txt")
four <- read.table("454.txt")

brks <- cut(four$V1, cutlevels, labels=la, include.lowest=TRUE)
fourvals <- tapply(four$V2, brks, mean)

thick<-2
plot(fourvals,
     axes=FALSE,
     ylim=c(0,50),
     type="l",
     col="red",
     xlab="Read position",
     ylab="Mean Q score",
     lwd=thick,
     oma=c(0,0,0,0))

#    main="b) Empirical per-base quality scores")

ilmnbrks <- cut(ilmn$V1, cutlevels, labels=la, include.lowest=TRUE)
ilmnvals <- tapply(ilmn$V2, ilmnbrks, mean)

iontbrks <- cut(iont$V1, cutlevels, labels=la, include.lowest=TRUE)
iontvals <- tapply(iont$V2, iontbrks, mean)

points(iontvals, col="purple", type="l", lwd=thick)
points(ilmnvals, col="blue", type="l", lwd=thick)

axis(1, 1:nlevels(brks), levels(brks))

axis(2)
legend("bottomleft", c("454 Junior", "MiSeq", "Ion Torrent"), col=c("red", "blue", "purple"), lty=1, lwd=thick)
par(xpd=FALSE)
text(38, 48, "b", font=2)

plot_empirical_vs_stated_quals("HPA2_4/readgroup1.QualityScoreCovariate.dat", "c1. 454 Junior")
plot_empirical_vs_stated_quals("ILMN1_L5/readgroup1.QualityScoreCovariate.dat", "c2. MiSeq")
par(xpd=FALSE)
text(-5, 52, "c", font=2)
plot_empirical_vs_stated_quals("BHAM5_6_1.5/readgroup1.QualityScoreCovariate.dat", "c3. Ion Torrent")

dev.off()

