#' @export

FDR_noPrint <- function(filename_pscore, filename_raw, vp.chrom, vp.pos, erase, wind, REs, cyc, perc) 
{
	#1. getting the top perc (decided by user) for each chromosome in original p score file
	chr <- unique(filename_raw[,1]);
	#cutoff for real data (for each chromosome)
	COs <- c();
	above_co <- c(); #it gets the number of RE sites above the cutoff for each chromosome (its size is the number of chromosomes)
	sum_FDRs <- rep(0,length(chr)); #starts with an array of zeros, then gets the sum of the FDRs
	for (l in 1:length(chr))    
	{
		chrom_cut <- quantile(filename_pscore[filename_pscore[,1]==l,3],perc);
		above_co <- c(above_co,sum(filename_pscore[filename_pscore[,1]==l,3]>=chrom_cut));
		COs <- c(COs,chrom_cut);
	}
	
	#2. randomizes the raw data, calculates p scores, and 
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
		ps <- pScore(mix,vp.chrom,vp.pos,erase,wind,REs);
		for(m in 1:length(chr))
		{
			above_co_mix <- c(above_co_mix,sum(ps[ps[,1]==m,3]>=COs[m]));
		}
		#similarity <- (above_co_mix/above_co)*100;
		FDR <- above_co_mix/(above_co+above_co_mix);
		sum_FDRs <- sum_FDRs+FDR;
    }
	return(sum_FDRs/cyc);
}
