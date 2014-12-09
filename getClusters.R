workdir <- getwd()
myPath <- paste(workdir,"/cluster_",sep="")
setwd("tmp/")

clusters<-read.table("Clusters.txt")
clusters<-clusters[,1]
Clusters <- clusters
        namesClusters <- as.numeric(names(table(Clusters)))
        relabelledClusters <- Clusters
        for (k in namesClusters[-1]){
                index <- which(Clusters == k)
                relabelledClusters[index] <- 1:length(index)    
        }
orderedNames <- order(namesClusters) - 1
names(orderedNames) <- namesClusters


#### this is to compute thetaNew, wNew and RPKM

nReads <- read.table("nReadsPerCluster.txt")
nReads <- colSums(nReads)
K<-length(relabelledClusters)
deFinal <- numeric(K)
thetaFinal <- numeric(K)
wFinal <- numeric(K)
Dead <- which(clusters == 0)
deFinal[Dead] <- rep(0,length(Dead))
thetaFinal[Dead] <- rep(1/nReads[1],length(Dead))
wFinal[Dead] <- rep(1/nReads[2],length(Dead))
setwd("clusters/")
for(k in 1:length(list.files(pattern = "cluster_"))){
	path <- paste("cluster_",k,"/state_vector.txt",sep="");
	models <- as.numeric(read.table(path)[,1]);
	models <- models[-c(1,length(models))];
	path <- paste("cluster_",k,"/theta.txt",sep="")
	th <- as.numeric(read.table(path)[,1]);
	th <- th[-c(1,length(th))]
	path <- paste("cluster_",k,"/w.txt",sep="")
	wou <- as.numeric(read.table(path)[,1]);
	wou <- wou[-c(1,length(wou))]
	clusterName <- as.numeric(names(orderedNames[k+1]))
	transcriptIndex <- which(Clusters == clusterName)
	deFinal[transcriptIndex] <- models
	thetaFinal[transcriptIndex] <- as.numeric(th)
	wFinal[transcriptIndex] <- as.numeric(wou)
	system(paste("rm -r cluster_",k,sep=""))
	if(k%%1000 == 0){print(k)}
}
setwd("../")
system("rm -r clusters")
system("rm *.panos *.prob *.temp Clusters.txt maxK1 maxK2")
setwd("../")

thetaFinal <- thetaFinal/sum(thetaFinal)
wFinal <- wFinal/sum(wFinal)
ramones <- cbind(log(thetaFinal),log(wFinal),deFinal)

sigTranscripts1 <- sigTranscripts2 <- sigTranscripts3 <- rep("non-DE",K)

# FDR: alpha = 0.01
alpha <- 0.01
p <- ramones
perm <- order(p[,3],decreasing = TRUE)
orderedP <- p[perm,3]
K <- dim(p)[1]
myList <-  1 - orderedP[1]
k <- 1 
criterion <- myList
while (criterion < alpha){
	k <- k + 1 
	myList <- myList + 1 - orderedP[k]
	criterion <- myList/k
}
if(k > 1){
	sigTranscripts1[perm[1:(k-1)]] <- rep("DE",k-1)
}
# FDR: alpha = 0.05
alpha <- 0.05
p <- ramones
perm <- order(p[,3],decreasing = TRUE)
orderedP <- p[perm,3]
K <- dim(p)[1]
myList <-  1 - orderedP[1]
k <- 1 
criterion <- myList
while (criterion < alpha){
	k <- k + 1 
	myList <- myList + 1 - orderedP[k]
	criterion <- myList/k
}
if(k > 1){
sigTranscripts2[perm[1:(k-1)]] <- rep("DE",k-1)
}

# FDR: alpha = 0.1
alpha <- 0.1
p <- ramones
perm <- order(p[,3],decreasing = TRUE)
orderedP <- p[perm,3]
K <- dim(p)[1]
myList <-  1 - orderedP[1]
k <- 1 
criterion <- myList
while (criterion < alpha){
	k <- k + 1 
	myList <- myList + 1 - orderedP[k]
	criterion <- myList/k
}
if(k > 1){
sigTranscripts3[perm[1:(k-1)]] <- rep("DE",k-1)
}
ramones <- cbind(ramones,sigTranscripts1, sigTranscripts2, sigTranscripts3 )

write.table(ramones,col.names = c("log-Theta","log-W","Prob(D.E)","fdr_0.01","fdr_0.05","fdr_0.1"),file = "estimates.txt",quote = FALSE)


system("rm *.prob parallelGNU.bash partitionMerge2.R sparse.so")





