#' @export

FDR_window_cmp_nokb <- function(filename_raw,small,big,step,vp.chrom,vp.pos,erase,cyc,CO_perc)
{
	ln <- length(unique(filename_raw[,1])); #getting the number of chromosomes
	min.FDR <- rep(2,ln+1);
	min.RE <- rep(0,ln+1);
	min.kb <- rep(0,ln+1);
	for(i in seq(small,big,by=step))
	{
		ps <- pScore_nokb(filename_raw,vp.chrom,vp.pos,erase,i/2);
		fdr.temp <- FDR_nokb_noPrint(ps,filename_raw,vp.chrom,vp.pos,erase,i/2,cyc,CO_perc);
		#getting the trans average
		fdr.sum <- sum(fdr.temp)-fdr.temp[vp.chrom];
		fdr.avg <- fdr.sum/(ln-1);
		fdr.perc <- fdr.avg/100; #turning the average into a number between 0 and 1 in order to compare to min.FDR
		
		#comparing the values of the new FDR to the minimum, if the new FDRs are lower than the minimum then they will become the new minimum
		for(k in 1:ln)
		{
			if(fdr.temp[k] < min.FDR[k])
			{
				min.FDR[k] <- fdr.temp[k];
				min.RE[k] <- i
			}
		}
		
		if(fdr.perc < min.FDR[ln+1])
		{
			min.FDR[ln+1] <- fdr.avg;
			min.RE[ln+1] <- i
		}
	}
	
	for(m in 1:ln)
	{	
		cat("\nThe smallest FDR for chromosome",m,"is: ",min.FDR[m],"\nThe window size that gives us this FDR is of size",min.RE[m],"REs\n\n");	
	}
	cat("\nThe smallest average of FDR for trans is: ",min.FDR[ln+1],"\nThe window size that gives us this average of FDRs is of size",min.RE[ln+1],"REs\n\n");
}
