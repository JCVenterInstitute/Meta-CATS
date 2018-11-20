#library(gregmisc)

mc=function(mat){

nrow=nrow(mat)
ncol=ncol(mat)
result.pvalue=matrix(rep(NA,ncol*ncol),ncol=ncol)

for(ii in 1 : (ncol - 1)){
	#print(ii)
	for(jj in (ii+1):ncol){
		#print(jj)
		obs=rbind(mat[,ii],mat[,jj])
		#print(obs)
		temp1=rep(NA,nrow)
		temp2=rep(NA,nrow)

		for(ss in 1:nrow){
			temp1[ss]=sum(mat[,ii])*sum(mat[ss,ii]+mat[ss,jj])/(sum(mat[,ii])+sum(mat[,jj]))
			temp2[ss]=sum(mat[,jj])*sum(mat[ss,ii]+mat[ss,jj])/(sum(mat[,ii])+sum(mat[,jj]))
		}
		exp=rbind(temp1,temp2)
		#print(exp)

		stat=sum((obs-exp)^2/exp)
		#print(stat)

		result.pvalue[jj,ii]=pchisq(stat, df=(ncol-1)*(nrow-1), lower.tail = FALSE)

		#print(result.pvalue)

		if(result.pvalue[jj,ii] == 0) {
			result.pvalue[jj,ii] = 1
		}
		result.pvalue[ii,jj]=result.pvalue[jj,ii]
		}

}
diag(result.pvalue) = 1
return(result.pvalue)

}


residueDiversity = function(contable) {
retval <- c("")
nresidue = nrow(contable)
ngroup = ncol(contable)
for(i in 1 : ngroup) {
	retval = paste(retval, "group", i, "(", sep="")

	for(j in 1 : nresidue) {
		residue = rownames(contable)[j]
		count = contable[rownames(contable)==residue,i]
		if(count > 0) {
			retval = paste(retval, count, "_", residue, ",_", sep="")
		}
	}
	retval = paste(retval, ")", sep="")
	if(i < ngroup) {
		retval = paste(retval, "|", sep="")
	}

}
retval = gsub(",_)", ")", retval, fixed=TRUE)
return (retval)
}

resultCST<-"Chi-square Analysis Result:\n"
resultCST<-"Site Number\tChi-Square Score\tP-Value\tDegrees of Freedom\tSparse Table (i.e. <5 of any residue)\tResidue Diversity Between Groups\n"
resultMC<-"MGC Multiple Comparison Result:\n"
resultMC<-"Site Number\tP-Value for Pairwise Comparison\tGroup 1 in Comparison\tGroup 2 in Comparison\n"
sigpvals <- numeric(0)
positions <- numeric(0)
flagsparse = "N"
residueDiv = ""
pvalcutoff <- 0.05
#inFilename <- 'rMsaInput.txt'
msaTable = read.table(inFilename, sep="", as.is=TRUE, comment.char="")
temp1=c(as.matrix(msaTable))
temp2=replace(temp1, which(temp1=="TRUE"),"T")
temp2=replace(temp2, which(temp2=="FALSE"),"F")
temp3=matrix(temp2,ncol=ncol(msaTable))
msaTableNew=as.data.frame(temp3)
n=ncol(msaTableNew)-2
names(msaTableNew)=c("mgcId","group",paste("Site",1:n, sep=""))

for(pos in 1:n){
	contable = table(msaTableNew[[paste("Site",pos,sep="")]],msaTableNew[["group"]])
	if(any(rownames(contable)=='#')) {contable=contable[-(which(rownames(contable)=='#')),,drop=FALSE]}
	if(nrow(contable)==0) {next}
	if(nrow(contable)==1) {
		resultCST=paste(resultCST,pos, "\tNA\t1\tNA\tN\tNA\n")
	} 
	else{
		if(any(contable[contable[,]>0] < 5)) {flagsparse = "Y"} else {flagsparse = "N"}
		residueDiv = residueDiversity(contable)
		fit=chisq.test(contable + 0.001)
		if(fit$p.value < pvalcutoff){
   			sigpvals <- c(sigpvals, fit$p.value)
   			positions <- c(positions, pos)
		}
		resultCST=paste(resultCST,pos,"\t",fit$statistic,"\t",fit$p.value,"\t",fit$parameter,"\t",flagsparse,"\t",residueDiv, "\n")
	}

	mcMat = mc(contable + 0.001)
	contRowCount = nrow(mcMat)

	for(j in 1:(contRowCount - 1)) {
		for(k in (j + 1):contRowCount) {
     			if(j < k) {
        			resultMC=paste(resultMC,pos,"\t",mcMat[j,k],"\t",j,"\t",k,"\n")
      		}
		}
	}

}
outfilename1 <- paste(inFilename, 'rResultChisqTest.txt', sep="-")
outfilename2 <- paste(inFilename, 'rResultMGCStat.txt', sep="-")
write.table(as.data.frame(resultCST),file=outfilename1)
write.table(as.data.frame(resultMC),file=outfilename2)

rm(resultCST, resultMC, temp1, temp2, temp3)
