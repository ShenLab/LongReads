
pdf(file='FP_Main.pdf', colormodel="cmyk", pointsize=4, width=9,height=3)

yax = "#False Positives"
xax = "Coverage"
par(mfrow=c(1,3))

for (prev in "10.0%") {
for (err in c("1.000%","2.00%","4.00%")) {

ttl = paste0("err=",err,", prev=",prev, collapse="")
sens_table <- read.table('table.tsv', col.names=c("Readlen","Errrate","Coverage","Prevalence","HMTP","HMFP","MTTP","MTFP"),header=FALSE, sep="\t")
five <- sens_table[which(sens_table$Prevalence==prev & sens_table$Errrate==err & sens_table$Readlen==500),]
ten <- sens_table[which(sens_table$Prevalence==prev & sens_table$Errrate==err & sens_table$Readlen==1000),]
twenty <- sens_table[which(sens_table$Prevalence==prev & sens_table$Errrate==err & sens_table$Readlen==2000),]
forty <- sens_table[which(sens_table$Prevalence==prev & sens_table$Errrate==err & sens_table$Readlen==4000),]

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

if(err=="4.00%") {
legend(30,yrange[2],c("HapMut-500bp","HapMut-1000bp","HapMut-2000bp","HapMut-4000bp","Mutect-500bp","Mutect-1000bp","Mutect-2000bp","Mutect-4000bp"),cex=1.6,col=colors,pch=plotchar,lty=linetype)
}
}
}
dev.off()



pdf(file='FP_Full.pdf', colormodel="cmyk", pointsize=4, width=6,height=9)

yax = "#False Positives"
xax = "Coverage"
par(mfrow=c(3,2))

for (err in c("1.000%","2.00%","4.00%")) {
for (prev in c("10.0%","25.0%")) {

ttl = paste0("err=",err,", prev=",prev, collapse="")
sens_table <- read.table('table.tsv', col.names=c("Readlen","Errrate","Coverage","Prevalence","HMTP","HMFP","MTTP","MTFP"),header=FALSE, sep="\t")
five <- sens_table[which(sens_table$Prevalence==prev & sens_table$Errrate==err & sens_table$Readlen==500),]
ten <- sens_table[which(sens_table$Prevalence==prev & sens_table$Errrate==err & sens_table$Readlen==1000),]
twenty <- sens_table[which(sens_table$Prevalence==prev & sens_table$Errrate==err & sens_table$Readlen==2000),]
forty <- sens_table[which(sens_table$Prevalence==prev & sens_table$Errrate==err & sens_table$Readlen==4000),]

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

if(err=="4.00%" & prev=="10.0%") {
legend(30,yrange[2],c("HapMut-500bp","HapMut-1000bp","HapMut-2000bp","HapMut-4000bp","Mutect-500bp","Mutect-1000bp","Mutect-2000bp","Mutect-4000bp"),cex=1.6,col=colors,pch=plotchar,lty=linetype)
}
}
}

