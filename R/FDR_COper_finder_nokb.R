#' @export

#this function recieves from user a constant FDR value, a starting cutoff percentile, and p-score data
#the function calculates the FDR with the given cutoff percentile and compares to the given constant FDR
#as long as the FDR calculated are lower than the given constant FDR, the cutoff will be lowered (by a percentage given by user)
#this will continue until the given constant FDR is surpassed (by at least one chromosome)
#this is performed as many times (cycles) requested by user
#either the average or the median of the results of cutoff percentiles that give the FDRs which are closest to the given constant FDR and the FDRs that are closest themselves are printed to screen
#this function is meant for p-score files that were calculated with only RE sites (usually when a rearrangement occured)

FDR_COper_finder_nokb <- function(filename_pscore, filename_raw, vp.chrom, vp.pos, erase, REs, cyc, perc, fdr, rem=0.1) 
{
	#1. getting the top perc (decided by user) for each chromosome in original p score file
	chr <- unique(filename_raw[,1])
	ln <- length(chr)
	fdr_arr <- rep(fdr,ln)
	#cutoff for real data (for each chromosome)
	COs <- c()
	above_co <- c() #it gets the number of RE sites above the cutoff for each chromosome (its size is the number of chromosomes)
	sum_FDRs <- rep(0,ln) #starts with an array of zeros, then gets the sum of the FDRs
	sum_COs <- rep(0,ln) #starts with an array of zeros, then gets the sum of cutoff percentiles
	all_FDRs <- c() #starts empty, then gets all the final results of FDRs for each cycle
	all_COs <- 	c() #starts empty, then gets all the final results of cutoffs for each cycle
	for (l in 1:ln)    
	{
		chrom_cut <- quantile(filename_pscore[filename_pscore[,1]==l,3],perc)
		above_co <- c(above_co,sum(filename_pscore[filename_pscore[,1]==l,3]>=chrom_cut))
		COs <- c(COs,chrom_cut)
	}
	
	#2. randomizes the raw data, calculates p scores, and 
	#repeats the number of cycles give by user, in order to average the FDRs 
	fail <- 0
    for (k in 1:cyc)
    {
		pmixAll <- c()
		above_co_mix <- c() #it gets the number of RE sites above the cutoff for each chromosome from the mixed data(its size is the number of chromosomes)
		#randomize (mix)
        for (r in 1:ln)
        {
            pmix <- sample(filename_raw[filename_raw[,1] == r,][,3]) 
            pmixAll <- c(pmixAll,pmix)
        }
        mix <- cbind(filename_raw[,1],filename_raw[,2],pmixAll)
		ps <-pScore_nokb(mix,vp.chrom,vp.pos,erase,REs)
		
		for(m in 1:ln)
		{
			above_co_mix <- c(above_co_mix,sum(ps[ps[,1]==m,3]>=COs[m]))
		}
		fdrs_cur <- above_co_mix/(above_co+above_co_mix)
		
		if(all(fdrs_cur <= fdr_arr))
		{
			while(all(fdrs_cur <= fdr_arr))
			{
				COs_temp <- COs
				fdrs_temp <- fdrs_cur
				COs <- COs_temp*(1-rem)
				above_co <- c()
				above_co_mix <- c()
				for(t in 1:ln)
				{
					above_co <- c(above_co,sum(filename_pscore[filename_pscore[,1]==t,3]>=COs[t]))
					above_co_mix <- c(above_co_mix,sum(ps[ps[,1]==t,3]>=COs[t]))
				}
				fdrs_cur <- above_co_mix/(above_co+above_co_mix)
			}
			#getting the sums of FDRs and COs (meant later for the averages)
			sum_FDRs <- sum_FDRs+fdrs_temp
			sum_COs <- sum_COs+COs_temp
			#recording the final FDRs and COs (meant later for the medians)
			all_FDRs <- rbind(all_FDRs,fdrs_temp)
			all_COs <- rbind(all_COs,COs_temp)
		}
		#if this cycle produces random data p scores that straight away have at least one chromosome with an FDR score above the FDR given by user
		#this cycle will be discarded and fdrs_temp and COs_temp will get all zeros
		else
		{
			fail <- fail + 1
		}
    }
	
	#checking to see that we got at least one instance of all the chromosomes random p scores under the FDR given by user
	#if there are none (meaning that at least one chromosome is above the FDR) then we tell the user,
	#and the user should either try again (since we randomize the data so this could be random), choose a new FDR, choose a higher cutoff percentile to start with, remove that one chromosome, not use this data, or check to see that there is nothing wrong with the data
	if(cyc == fail)
	{
		cat("\nthere were no instances where the random data produced p-scores that gave an FDR below",fdr,"\n please try a new FDR or a higher p-score cutoff percentile\n\n")
	}
	else
	{
	#averages and medians:
		#averages:
		#printing the average of FDRs from all cycles, that were the closest to the users input FDR before failing and going above it
		avg_FDRs <- sum_FDRs/(cyc-fail)
		avg_FDRs <- t(as.data.frame(avg_FDRs))
		colnames(avg_FDRs) <- chr
		cat("\nthe averages of FDRs that are closest to your input FDR of ",fdr,":\n",sep="")
		print(avg_FDRs)

		#printing the lowest cutoff percentiles (the average of all cycles) that give the lowest closest values of FDR to the given FDR
		avg_COs <- sum_COs/(cyc-fail)		
		avg_COs_per <- rep(0,ln)
		for(b in 1:ln)
		{
			avg_COs_per[b] <- ecdf(filename_pscore[filename_pscore[,1]==b,3])(avg_COs[b])
		}
		avg_COs_per <- t(as.data.frame(avg_COs_per))
		colnames(avg_COs_per) <- chr
		cat("\nthe averages of cutoff percentiles that give the closest FDR to your input FDR of ",fdr,":\n",sep="")
		print(avg_COs_per)
		
		#medians:
		#printing the median of FDRs from all cycles, that were the closest to the users input FDR before failing and going above it
		med_FDRs <- apply(all_FDRs,2,median)
		med_FDRs <- t(as.data.frame(med_FDRs))
		colnames(med_FDRs) <- chr
		cat("\nthe medians of FDRs that are closest to your input FDR of ",fdr,":\n",sep="")
		print(med_FDRs)
		
		#printing the lowest cutoff percentiles (the median of all cycles) that give the lowest closest values of FDR to the given FDR
		med_COs <- apply(all_COs,2,median)
		med_COs_per <- rep(0,ln)
		for(b in 1:ln)
		{
			med_COs_per[b] <- ecdf(filename_pscore[filename_pscore[,1]==b,3])(med_COs[b])
		}
		med_COs_per <- t(as.data.frame(med_COs_per))
		colnames(med_COs_per) <- chr		
		cat("\nthe medians of cutoff percentiles that give the closest FDR to your input FDR of ",fdr,":\n",sep="")
		print(med_COs_per)			
	}
}
