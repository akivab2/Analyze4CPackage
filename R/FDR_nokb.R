#' @export

FDR_nokb <- function(filename_pscore, filename_raw, vp.chrom, vp.pos, erase, REs, cyc, perc) 
{
	#1. getting the top perc (decided by user) for each chromosome in original p score file
	chr <- unique(filename_raw[,1]);
	#cutoff for real data (for each chromosome)
	COs <- c();
	above_co <- c(); #it gets the number of RE sites above the cutoff for each chromosome (its size is the number of chromosomes)
	sum_FDRs <- rep(0,length(chr)); #starts with an array of zeros, then gets the sum of the similarities
	for (l in 1:length(chr))    
	{
		chrom_cut <- quantile(filename_pscore[filename_pscore[,1]==l,3],perc);
		above_co <- c(above_co,sum(filename_pscore[filename_pscore[,1]==l,3]>=chrom_cut));
		COs <- c(COs,chrom_cut);
	}
	COs_forprint <- t(as.data.frame(COs))
	colnames(COs_forprint) <- chr
	cat ("\nthe top ",(1-perc)*100,"% p-scores for each chromosome start at:\n",sep="")
	print(COs_forprint)
	cat("\n")
	
	#2. randomizes the raw data, calculates p scores, and gets the FDRs for each chromosome in each cycle
	#repeats by an amount chosen by user, in order to average the FDRs 
    for (k in 1:cyc)
    {
		pmixAll <- c();
		above_co_mix <- c(); #it gets the number of RE sites above the cutoff for each chromosome from the mixed data(its size is the number of chromosomes)
		#randomize (mix)
        for (r in 1:length(chr))
        {
            pmix <- sample(filename_raw[filename_raw[,1] == r,][,3]); 
            pmixAll <- c(pmixAll,pmix)
        }
        mix <- cbind(filename_raw[,1],filename_raw[,2],pmixAll);
		ps <- pScore_nokb(mix,vp.chrom,vp.pos,erase,REs);
		for(m in 1:length(chr))
		{
			above_co_mix <- c(above_co_mix,sum(ps[ps[,1]==m,3]>=COs[m]));
		}
		#fdr <- (above_co_mix/above_co)*100;
		fdr <- above_co_mix/(above_co+above_co_mix);
		cat("cycle",k,":",fdr,"\n");
		sum_FDRs <- sum_FDRs+fdr;
    }
	avg_FDRs <- sum_FDRs/cyc
	avg_FDRs <- t(as.data.frame(avg_FDRs))
	colnames(avg_FDRs) <- chr
	cat("\nthe averages of FDRs are:\n")
	print(avg_FDRs)
}
