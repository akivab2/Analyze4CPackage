#' @export

#this is calculates the p-score from the file of number of reads for each RE site (raw data).
#this uses a window size only according to number of RE sites
#note: this is a slightly improved version and more quicker than the older pScore function
pScore_nokb <- function(filename,vp.chrom,vp.pos,erase,REs)
{
	dat<-filename
	#for now i'm not going to use the removing of the sections in cis. it isn't compatible with some other functions in the program
	#dat<-dat[dat[,1]!=vp.chrom | (dat[,1]==vp.chrom & (dat[,2]<(vp.pos-erase) | dat[,2]>(vp.pos+erase))),]
	nothing <- erase #this is just so there is a use for 'erase'. it is not really needed, and should be taken off once the line above is returned	
	sizdat <- nrow(dat)
	pVal <- rep(0,sizdat)
	chr <- unique(dat[,1])
	ind <- 1
	for (j in 1:length(chr))
	{
		tmp.data<-dat[dat[,1]==chr[j],]
		prob<-sum(tmp.data[,3]>0)/nrow(tmp.data)
		for (i in 1:nrow(tmp.data))
		{
			first<-i-REs
			last<-i+REs
			if (first<1){first<-1}
			if (last>nrow(tmp.data)){last<-nrow(tmp.data)}
						tmp<-tmp.data[first:last,]
			pVal[ind] <- binom.test(sum(tmp[,3]>0),length(tmp[,3]),prob,alternative="greater")$p.value
			ind <- ind + 1
		}
	}
	ps<- -log10(pVal)
	ps[ps == Inf] <- max(ps[ps != Inf])
	ps <- cbind(dat[,1],dat[,2],ps)
	return(ps)
}
