
#ttl = "False positive rate"
yax = "#False positives"
yax = ""
xax = "Coverage"

five <- read.table('_500_.txt', col.names=c("Readlen","Errrate","Coverage","Prevalence","HMTP","HMFP","MTTP","MTFP"),header=FALSE, sep="\t")
ten <- read.table('_1000_.txt', col.names = c("Readlen","Errrate","Coverage","Prevalence","HMTP","HMFP","MTTP","MTFP"),header=FALSE,sep="\t")
twenty <- read.table('_2000_.txt', col.names = c("Readlen","Errrate","Coverage","Prevalence","HMTP","HMFP","MTTP","MTFP"),header=FALSE,sep="\t")
forty <- read.table('_4000_.txt', col.names = c("Readlen","Errrate","Coverage","Prevalence","HMTP","HMFP","MTTP","MTFP"),header=FALSE,sep="\t")

xrange <- range(five[,3])
yrange <- range(0,five[,6],ten[,6],twenty[,6],forty[,6],five[,8],ten[,8],twenty[,8],forty[,8],140)
plot(xrange,yrange,type="n",xlab=xax,ylab=yax,main=ttl,cex.lab=2.0, cex.main=2.0,cex.axis=2.0)
colours <- rainbow(12)
colors <- c(colours[2],colours[3],colours[6],colours[7],colours[8],colours[9],colours[10],colours[11])
linetype <- c(3,5,6,7,3,5,6,7)
plotchar <- c(seq(21,24,1),seq(21,24,1))
lines(five$Coverage,five$HMFP, type="b", lwd = 1.5, lty = linetype[1], col = colors[1], pch=plotchar[1])
lines(ten$Coverage,ten$HMFP, type="b", lwd = 1.5, lty = linetype[2], col = colors[2], pch=plotchar[2])
lines(twenty$Coverage,twenty$HMFP, type="b", lwd = 1.5, lty = linetype[3], col = colors[3], pch=plotchar[3])
lines(forty$Coverage,forty$HMFP, type="b", lwd = 1.5, lty = linetype[4], col = colors[4], pch=plotchar[4])
lines(five$Coverage,five$MTFP, type="b", lwd = 1.5, lty = linetype[5], col = colors[5], pch=plotchar[5])
lines(ten$Coverage,ten$MTFP, type="b", lwd = 1.5, lty = linetype[6], col = colors[6], pch=plotchar[6])
lines(twenty$Coverage,twenty$MTFP, type="b", lwd = 1.5, lty = linetype[7], col = colors[7], pch=plotchar[7])
lines(forty$Coverage,forty$MTFP, type="b", lwd = 1.5, lty = linetype[8], col = colors[8], pch=plotchar[8])

#legend(xrange[1]+18.5,yrange[1]+140,c("HapMut-500bp","HapMut-1000bp","HapMut-2000bp","HapMut-4000bp","Mutect-500bp","Mutect-1000bp","Mutect-2000bp","Mutect-4000bp"),cex=1.6,col=colors,pch=plotchar,lty=linetype)

