args <- commandArgs(trailingOnly = TRUE)
outputfilename <- (args[1])
h <- as.numeric(args[2])
simpletree <- (args[3])
path <- (args[4])

library("ape", lib.loc=path)

treeinfo <- read.tree(text = simpletree)
leafHeight = 40
png(file=paste(outputfilename,".png",sep =""),width=2000, height=h*leafHeight, res=300)
par(mar=c(1,1,1,1))
plot(treeinfo, cex=0.5)
garbage <- dev.off()