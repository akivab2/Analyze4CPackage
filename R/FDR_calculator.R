#' @export

FDR_calculator <- function(filename_raw,small.kb,big.kb,step.kb,small.CO,big.CO,step.CO,vp.chrom,vp.pos,erase,ratio,cyc)
{
	sz <- (floor((big.kb-small.kb)/step.kb)+1)*(round((big.CO-small.CO)/step.CO)+1); #getting the number of all the elements that will be checked
	#sz <- (floor((big.kb-small.kb)/step.kb)+1)*(floor((big.CO-small.CO)/step.CO)+1); #getting the number of all the elements that will be checked
	ind_count <- 0; #this counts the index for each column in the data frame
	window_kb <- rep(0,sz);
	window_RE <- rep(0,sz);
	CO <- rep(0,sz);
	fdr.avg.trans <- rep(0,sz);
	fdrs <- c()
	counter <- 0		
	ln <- length(unique(filename_raw[,1])); #getting the number of chromosomes
	for(i in seq(small.kb,big.kb,by=step.kb))
	{
		res_num <- ceiling(i/ratio); #getting the RE window size
		ps <- pScore(filename_raw,vp.chrom,vp.pos,erase,i,ceiling(res_num/2)); #calculating p score with current window size
		for(j in seq(small.CO,big.CO,by=step.CO))
		{
			counter <- counter + 1
			cat("\ncalculating",cyc,"cycles for round",counter,"out of",sz,"\n\n")		
			fdr.temp <- FDR_noPrint(ps,filename_raw,vp.chrom,vp.pos,erase,i,ceiling(res_num/2),cyc,j); #calculating FDR with current cutoff percent
			fdr.sum <- sum(fdr.temp)-fdr.temp[vp.chrom]; #getting the sum of fdr from trans
			fdr.avg <- fdr.sum/(ln-1); #getting the average of fdr from trans
			
			#creating the columns for the data frame
			ind_count <- ind_count + 1;
			window_kb[ind_count] <- i;
			window_RE[ind_count] <-	res_num;
			CO[ind_count] <- j;	
			fdr.avg.trans[ind_count] <- fdr.avg;
			fdrs <- rbind(fdrs,fdr.temp)
		}
	}
	#getting the chromosome names/numbers and adding them as the column names in 'fdrs'
	chrs <- unique(filename_raw[,1])
	chr_names <- rep(0,ln)
	for(t in 1:ln)
	{
		chr_names[t] <- paste("chr",chrs[t],sep="")
	}
	fdrs <- t(as.data.frame(fdrs))
	rownames(fdrs) <- chr_names
	fdrs <- t(fdrs)
	rownames(fdrs) <- seq(1:nrow(fdrs))
	
	all <- data.frame(window_kb,window_RE,CO,fdrs,fdr.avg.trans); #creating the data frame
	print(all)	
	return(all);
}
