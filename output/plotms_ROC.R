
args <- commandArgs(trailingOnly = TRUE)

prev <- args[1]
cov <- args[2]

filename <- paste0("ROC_",prev,"_",cov,".pdf", collapse="")

if(as.numeric(prev) == 10)
  pre = "15"
if(as.numeric(prev) == 25)
  pre = "35"

#     pdf(file = ifelse(onefile, "Rplots.pdf", "Rplot%03d.pdf"),
#         width, height, onefile, family, title, fonts, version,
#         paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
#         useDingbats, useKerning, fillOddEven, compress)

pdf(file=filename, colormodel="cmyk", pointsize=4)
#jpeg(file=filename)
xax = "#False positives"
yax = "Sensitivity"
par(mfrow=c(2,2))

for(rlen in c(500,1000,2000,4000)) {
bhm1 <- read.table(paste0("0.005_0.005_",cov,"_",rlen,"_",pre,"_21.c", collapse=""), header=TRUE, sep="\t")
bhm2 <- read.table(paste0("0.01_0.01_",cov,"_",rlen,"_",pre,"_21.c", collapse=""), header=TRUE, sep="\t")
bhm3 <- read.table(paste0("0.02_0.02_",cov,"_",rlen,"_",pre,"_21.c", collapse=""), header=TRUE, sep="\t")
mt1 <- read.table(paste0("0.005_0.005_",cov,"_",rlen,"_",pre,"_21.n", collapse=""), header=TRUE, sep="\t")
mt2 <- read.table(paste0("0.01_0.01_",cov,"_",rlen,"_",pre,"_21.n", collapse=""), header=TRUE, sep="\t")
mt3 <- read.table(paste0("0.02_0.02_",cov,"_",rlen,"_",pre,"_21.n", collapse=""), header=TRUE, sep="\t")

hm1 <- bhm1[6:100,]
hm2 <- bhm2[6:100,]
hm3 <- bhm3[6:100,]

hm1[is.na(hm1$False_positives),8] <- 0
hm2[is.na(hm2$False_positives),8] <- 0
hm3[is.na(hm3$False_positives),8] <- 0
mt1[is.na(mt1$FalsePositives),3] <- 0
mt2[is.na(mt2$FalsePositives),3] <- 0
mt3[is.na(mt3$FalsePositives),3] <- 0
hm1[is.na(hm1$Sensitivity),3] <- 0
hm2[is.na(hm2$Sensitivity),3] <- 0
hm3[is.na(hm3$Sensitivity),3] <- 0
mt1[is.na(mt1$Sensitivity),2] <- 0
mt2[is.na(mt2$Sensitivity),2] <- 0
mt3[is.na(mt3$Sensitivity),2] <- 0

xrange <- range(0,hm1$False_positives,hm2$False_positives,hm3$False_positives,mt1$FalsePositives,mt2$FalsePositives,mt3$FalsePositives)
yrange <- c(0.0,1.0)
ttl = paste0("prev=",prev,";cov=",cov, ";rlen=",rlen,collapse="")

plot(xrange,yrange,type="n",xlab=xax,ylab=yax,main=ttl,cex.lab=2.0, cex.main=2.0,cex.axis=2.0)
colours <- rainbow(12)
colors <- c(colours[2],colours[3],colours[6],colours[7],colours[8],colours[9],colours[10],colours[11])
linetype <- c(3,5,6,3,5,6)
plotchar <- c(seq(21,23,1),seq(21,23,1))
lines(hm1$False_positives, hm1$Sensitivity, type="b", lwd = 1.5, lty = linetype[1], col = colors[1], pch=plotchar[1])
lines(hm2$False_positives, hm2$Sensitivity, type="b", lwd = 1.5, lty = linetype[2], col = colors[2], pch=plotchar[2])
lines(hm3$False_positives, hm3$Sensitivity, type="b", lwd = 1.5, lty = linetype[3], col = colors[3], pch=plotchar[3])
lines(mt1$FalsePositives, mt1$Sensitivity, type="b", lwd = 1.5, lty = linetype[1], col = colors[4], pch=plotchar[4])
lines(mt2$FalsePositives, mt2$Sensitivity, type="b", lwd = 1.5, lty = linetype[2], col = colors[5], pch=plotchar[5])
lines(mt3$FalsePositives, mt3$Sensitivity, type="b", lwd = 1.5, lty = linetype[3], col = colors[6], pch=plotchar[6])

if(rlen==2000) {
  legend(40,1.0,c("HapMut-1%","HapMut-2%","HapMut-4%","Mutect-1%","Mutect-2%","Mutect-4%"),cex=1.6,col=colors,pch=plotchar,lty=linetype)
}
}


