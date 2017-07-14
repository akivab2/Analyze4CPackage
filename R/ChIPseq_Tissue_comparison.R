#' @export

#this function compares one ChIPseq (peaks) file with 2 different contact bands files
#it will find the difference between the comparison of each contact bands data with the ChIPseq data

#need to add an option of saving the intersection files in order to view in the browser - probably don't need this, since i could just take the contact bands file and the peaks file and apply the cutoff in the browser

ChIPseq_Tissue_comparison <- function(Experiments_4C,ChIPseqVScontacts_plots,rearranged_rawData)
{
	cat("\n**************************************************************************\n\n")
	cat("\n                        ChIPseq tissue comparison\n\n")
	cat("\n**************************************************************************\n\n")
	
	#choosing a genome size file
	ls_genomes <- system("ls ~/Analyze4C/genomes/",intern=TRUE)
	genome_file.name <- readline(prompt=cat("\nChoose the file of the organism you are dealing with:\n",ls_genomes,"",sep="\n"))
	genome_sizes <- read.table(paste("~/Analyze4C/genomes/",genome_file.name,sep=""))
	
	numOFchroms <- nrow(genome_sizes)
	chr_names <- rep(0,numOFchroms)
	for(t in 1:numOFchroms)
	{
		chr_names[t] <- paste("chr",t,sep="")
	}

	#getting or creating two contact bands files
	conts_filename <- rep(0,2)
	conts <- list(0,0)
	cis <- rep(0,2)
	sp1 <- list(0,0)
	for(p in 1:2)
	{
		#aking user if to use already existing contact band file or create a new one	
		ans1 <- as.numeric(readline(prompt=cat("\n\nfor contact bands file number",p,", would you like to:\n1) choose an already existing contact bands file\n2) create a new contact bands data file\n\n")))
		if(ans1 == 1) #taking an already existing contact bands file
		{
			ls_conts <- system("ls ~/Analyze4C/contact_bands/",intern=TRUE)
			repeat
			{
				if(length(ls_conts) != 0)
				{	
					conts_filename[p] <- readline(prompt=cat("\nchoose the file with contact bands file that you would like to use:\n",ls_conts,"",sep="\n"))
					ind_files1 <- pmatch(conts_filename[p],ls_conts,nomatch=0)
					if(ind_files1 == 0)
					{
						cat("no such file exists.\nplease try again.\n\n")
					}
					else
					{
						conts[[p]] <- read.table(paste("~/Analyze4C/contact_bands/",conts_filename[p],sep=""))
						break
					}
				}
				else
				{
					cat("\nno contact bands file exist\n")
					ans2 <- readline(prompt=cat("\nwould you like to create one?\ny/n\n\n"))
					if(ans2 == "y")
					{
						ans1 <- 2
						break
					}
				}
			}	
		}

		if(ans1 == 2) #creating new contact bands file
		{
			cat("\n**************************************************************************\n\n")
			cat("\n                        Creating contact bands\n\n")
			cat("\n**************************************************************************\n\n")		
			conts_filename[p] <- makeContactBands.mainMenu(Experiments_4C,rearranged_rawData)
			conts[[p]] <- read.table(paste("~/Analyze4C/contact_bands/",conts_filename[p],sep=""))
		}
				
		#finding the specific contact bands information in Experiments_4C
		sp1[[p]] <- unlist(strsplit(conts_filename[p],"_"))
		out <- findIn.Experiments_4C(sp1[[p]][1],sp1[[p]][2],sp1[[p]][3],Experiments_4C)

		#getting the cis chromosome number
		cis[p] <- as.numeric(out[2])
	}
	
	#get a ChIPseq peaks file:
	cat("\n\nchoosing the peaks file:\n\n")
	ls_peaks <- system("ls ~/Analyze4C/ChIPseq/peaks | grep 'tagsANDpscores.bed$'",intern=TRUE)
	
	if(length(ls_peaks) >= 1)
	{	
		repeat
		{
			peaks_filename <- readline(prompt=cat("\nchoose a peaks file that you would like to use:\n",ls_peaks,"",sep="\n"))
			ind_files2 <- pmatch(peaks_filename,ls_peaks,nomatch=0)
			if(ind_files2 == 0)
			{
				cat("no such file exists.\nplease try again.\n\n")
			}
			else
			{
				peaks <- read.table(paste("~/Analyze4C/ChIPseq/peaks/",peaks_filename,sep=""))
				peaks_tissue <- readline(prompt=cat("\nenter the tissue from where the ChIPseq data is taken:\n"))
				break
			}
		}	
	}
	else
	{
		cat("\nthere are no peaks bed files in the folder\n")
		cat("\nif you would like to execute this function, first download more ChIPseq datasets\n\n")
		return()
	}
	
	#choosing the value from which the cutoffs are created and the value which is going to be used for each peak above and equal to the cutoff
	cat("\napplying cutoffs to the peaks data:\n")
	ans9 <- as.integer(readline(prompt=cat("\nchoose an option:\n1) apply the cutoff on the p-scores, get the peak values from p-scores\n2) apply the cutoff on the p-scores, get the peak values from number of tags\n3) apply the cutoff on the number of tags, get the peak values from p-scores\n4) apply the cutoff on the number tags, get the peak values from number of tags\n\n")))
	if(ans9==1)
	{
		CO_col <- 5
		value_col <- 5
		CO_source <- "p-scores"
		value_source <- "p-scores"
	}
	else if(ans9==2)
	{
		CO_col <- 5
		value_col <- 4
		CO_source <- "p-scores"
		value_source <- "tags"			
	}
	else if(ans9==3)
	{
		CO_col <- 4
		value_col <- 5
		CO_source <- "tags"
		value_source <- "p-scores"			
	}
	else if(ans9==4)
	{
		CO_col <- 4
		value_col <- 4
		CO_source <- "tags"
		value_source <- "tags"			
	}
	
	#applying a cutoff to the peaks data
	CO_type_ans <- as.numeric(readline(prompt=cat("\nhow would you like to apply the cutoff?\n1) by percentage\n2) by a set peaks value\n\n")))
	if(CO_type_ans == 1)
	{
		peaks_COper <- as.numeric(readline(prompt=cat("\nenter a cutoff percentage (between 0 and 1) that will be applied to the peaks\nonly the peaks that are above or equal to the cutoff will be intersected to the contact bands.\nif you don't want any cutoff, enter 0:\n")))
		CO_type <- "by percentage"
		peaks_cutoff <- peaks_COper
		ans4 <- as.integer(readline(prompt=cat("\nshould the cutoff percentage be of:\n1) the whole genome\n2) all trans\n3) each chromosome separately\n\n")))
		if(ans4 == 1) #whole genome
		{
			CO <- quantile(peaks[,CO_col],peaks_COper)
			peaks_afterCO <- peaks[peaks[,CO_col]>=CO,c(1:3,value_col)]
		}
		else if(ans4 == 2) #trans
		{
			cis_ans <- as.integer(readline(prompt=cat("\nwhich cis chromosome would you like to use (enter the option number)?\n1)",conts_filename[1],"-",cis[1],"\n2)",conts_filename[2],"-",cis[2],"\n\n")))
			peaks_afterCO <- c()
			CO <- quantile(peaks[peaks[,1]!=cis[cis_ans],CO_col],peaks_COper)
			for(k in 1:numOFchroms) 
			{
				if(k!=cis[cis_ans])
				{
					peaks_temp <- peaks[peaks[,1]==k & peaks[,CO_col]>=CO,c(1:3,value_col)]
					peaks_afterCO <- rbind(peaks_afterCO,peaks_temp)
				}
				else
				{
					peaks_temp <- peaks[peaks[,1]==k,c(1:3,value_col)]
					peaks_afterCO <- rbind(peaks_afterCO,peaks_temp)
				}
			}	
		}
		else if(ans4 == 3) #specific chromosomes
		{		
			peaks_afterCO <- c()
			for(k in 1:numOFchroms) 
			{
				CO <- quantile(peaks[peaks[,1]==k,CO_col],peaks_COper)
				peaks_temp <- peaks[peaks[,1]==k & peaks[,CO_col]>=CO,c(1:3,value_col)]
				peaks_afterCO <- rbind(peaks_afterCO,peaks_temp)
			}		
		}
	}
	else
	{
		CO <- as.numeric(readline(prompt=cat("\nenter a cutoff (in p-score or number of tags) that will be applied to",peaks_filename,"\nonly the peaks that are above or equal to the cutoff will be intersected to the contact bands.\nif you don't want any cutoff, enter 0:\n")))
		CO_type <- "by value"
		peaks_cutoff <- CO
		ans4 <- as.integer(readline(prompt=cat("\nshould the cutoff be applied to:\n1) the whole genome\n2) all trans\n3) each chromosome separately\n\n")))
		if(ans4 == 1) #whole genome
		{
			peaks_afterCO <- peaks[peaks[,CO_col]>=CO,c(1:3,value_col)]
		}
		else if(ans4 == 2) #trans
		{
			cis_ans <- as.integer(readline(prompt=cat("\nwhich cis chromosome would you like to use (enter the option number)?\n1)",conts_filename[1],"-",cis[1],"\n2)",conts_filename[2],"-",cis[2],"\n\n")))
			peaks_afterCO <- c()
			for(k in 1:numOFchroms) 
			{
				if(k!=cis[cis_ans])
				{
					peaks_temp <- peaks[peaks[,1]==k & peaks[,CO_col]>=CO,c(1:3,value_col)]
					peaks_afterCO <- rbind(peaks_afterCO,peaks_temp)
				}
				else
				{
					peaks_temp <- peaks[peaks[,1]==k,c(1:3,value_col)]
					peaks_afterCO <- rbind(peaks_afterCO,peaks_temp)
				}
			}	
		}
		else if(ans4 == 3) #specific chromosomes
		{		
			peaks_afterCO <- c()
			for(k in 1:numOFchroms) 
			{
				peaks_temp <- peaks[peaks[,1]==k & peaks[,CO_col]>=CO,c(1:3,value_col)]
				peaks_afterCO <- rbind(peaks_afterCO,peaks_temp)
			}		
		}	
	}
	
	#making this into a list
	#for the rest of the function i must make sure that everytime there is 'peaks_afterCO' it should be double
	peaks_afterCO_temp <- list()
	peaks_afterCO_temp[[1]] <- peaks_afterCO
	peaks_afterCO_temp[[2]] <- peaks_afterCO
	peaks_afterCO <- peaks_afterCO_temp

	for(e in 1:2)
	{	
		#removing the translocated areas from the peaks data, the same areas that were removed from the contact bands files
		if(("rearranged" %in% sp1[[e]]) & !("removed" %in% sp1[[e]])) #if the conts filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
		{	
			#getting the line numbers that need to be removed
			for(d in 1:length(rearranged_rawData$Experiment))
			{
				sTemp <- unlist(strsplit(rearranged_rawData$Experiment[d],"_"))
				dtTemp <- unlist(strsplit(rearranged_rawData$Date_and_Time[d]," "))
				tTemp <- dtTemp[2]
				dtTemp[2] <- paste(unlist(strsplit(tTemp,":")),collapse="")
				if(all(sTemp %in% sp1[[e]]) & !is.na(rearranged_rawData$Added_Sections_Lines[d]) & all(dtTemp %in% sp1[[e]]))
				{
					lin_nums1 <- rearranged_rawData$Added_Sections_Lines[d]
					lin_nums2 <- as.integer(unlist(strsplit(lin_nums1,", ")))
					
					#getting the raw data file that the contact bands file was created from:
					if(!identical(grep("covRemoved",sp1[[e]]),integer(0)))
					{
						ls_raws <- system("ls ~/Analyze4C/rawData/coverage_removed",intern=TRUE)
					}
					else
					{
						ls_raws <- system("ls ~/Analyze4C/rawData/rearranged",intern=TRUE)
					}						
					dAndtspl <- paste(dtTemp[1],"_",dtTemp[2],sep="")
					raws_inds <- grep(dAndtspl,ls_raws)
					for(o in raws_inds)
					{
						if(identical(grep("removed",ls_raws[o]),integer(0)))
						{
							if(!identical(grep("covRemoved",sp1[[e]]),integer(0)))
							{
								#getting the raw data
								raw_dat <- read.table(paste("~/Analyze4C/rawData/coverage_removed/",ls_raws[o],sep=""))
								break
							}
							else
							{
								#getting the raw data
								raw_dat <- read.table(paste("~/Analyze4C/rawData/rearranged/",ls_raws[o],sep=""))
								break
							}
						}
					}			
					break
				}			
			}
			
			#the following actions will work under these specific conditions: 
			#1) the rearranged raw data file contains all the original chromosomes. 
			#2) the line numbers in 'rearranged_rawData' file are divided so that each pair are of one specific chromosome (meaning that 'first' and 'second' are not of two different chromosomes)

			#removing the sections from the peaks data
			all_chroms <- seq(numOFchroms)
			pointer <- 1
			peaks_fin <- c()
			for(r in seq(1,length(lin_nums2),by=2))
			{
				first <- lin_nums2[r]
				second <- lin_nums2[r+1]
				ch <- raw_dat[first,1]
				if(ch != raw_dat[second,1])
				{
					cat("\nerror:there is something wrong with the line numbers entered into the file 'rearranged_rawData'.\nthe first and last lines must be of the same chromosome\n\n")
					return()
				}
				
				#taking all the chromosomes that aren't getting the removal treatment, and adding them to the peaks data
				while(all_chroms[pointer]!=ch)
				{
					peaks_t <- peaks_afterCO[[e]][peaks_afterCO[[e]][,1]==all_chroms[pointer],]
					peaks_fin <- rbind(peaks_fin,peaks_t)
					pointer <- pointer + 1
				}
				
				peaks_t <- c() #just in case no conditions are met, this way we won't get the previous peaks_t (this shouldn't happen, it's just a precaution)
				#getting the indices to remove
				if(first==1)#if first is the 1st line in the file
				{
					if(second!=nrow(raw_dat)) #if second isn't the last of the file
					{			
						if(raw_dat[second,1]!=raw_dat[(second+1),1]) #if it covers all the chromosome (if second is the last in chromosome) (1)
						{
							peaks_t <- c()
						}
						else if(raw_dat[second,1]==raw_dat[(second+1),1]) #if second isn't the last of the chromosome (2)
						{
							#remove anything below the index at second
							peaks_t <- peaks_afterCO[[e]][(peaks_afterCO[[e]][,1]==ch)&(peaks_afterCO[[e]][,2]>=raw_dat[(second+1),2]),]
						}
					}	
					else #the line numbers cover all the data, so none of it is used (9)
					{
						peaks_t <- c()
					}				
				}
				else if(raw_dat[first,1]!=raw_dat[(first-1),1]) #if it is the first line in chromosome
				{
					if(second!=nrow(raw_dat)) #if second isn't last line
					{
						if(raw_dat[second,1]!=raw_dat[(second+1),1]) #if second is the last in chromosome (4)
						{
							peaks_t <- c()
						}
						else if(raw_dat[second,1]==raw_dat[(second+1),1]) #if second isn't last in chromosome (3)
						{
							#remove anything below the index at second
							peaks_t <- peaks_afterCO[[e]][(peaks_afterCO[[e]][,1]==ch)&(peaks_afterCO[[e]][,2]>=raw_dat[(second+1),2]),]
						}
					}
					else #if second is last line of the file (7)
					{
						peaks_t <- c()
					}	
				}
				else if(raw_dat[first,1]==raw_dat[(first-1),1])#if first is in the middle of a chromosome
				{
					if(second!=nrow(raw_dat)) #if second isn't last line
					{
						if(raw_dat[second,1]!=raw_dat[(second+1),1]) #if second is the last in chromosome (6)
						{
							#remove anything above the index at first
							peaks_t <- peaks_afterCO[[e]][(peaks_afterCO[[e]][,1]==ch)&(peaks_afterCO[[e]][,2]<=raw_dat[(first-1),2]),]
						}
						else if(raw_dat[second,1]==raw_dat[(second+1),1]) #if second isn't last in chromosome (5)
						{
							#remove anything above the index at first and below the index at second
							peaks_t <- peaks_afterCO[[e]][(peaks_afterCO[[e]][,1]==ch)&((peaks_afterCO[[e]][,2]<=raw_dat[(first-1),2])|(peaks_afterCO[[e]][,2]>=raw_dat[(second+1),2])),]
						}
					}
					else #if second is last line of the file (8)
					{
						#remove anything above the index at first
						peaks_t <- peaks_afterCO[[e]][(peaks_afterCO[[e]][,1]==ch)&(peaks_afterCO[[e]][,2]<=raw_dat[(first-1),2]),]
					}	
				}			
				peaks_fin <- rbind(peaks_fin,peaks_t)
				pointer <- pointer + 1
			}
			
			#adding the last chromosomes data (if there are any) that didn't get the section removal treatment
			while(pointer<=numOFchroms)
			{
				peaks_t <- peaks_afterCO[[e]][peaks_afterCO[[e]][,1]==all_chroms[pointer],]
				peaks_fin <- rbind(peaks_fin,peaks_t)
				pointer <- pointer + 1
			}
			
			peaks_afterCO[[e]] <- peaks_fin
			rownames(peaks_afterCO[[e]]) <- seq(length=nrow(peaks_afterCO[[e]]))
		}
	}

	write.table(peaks_afterCO[[1]],"~/Analyze4C/temp/peaks1_afterCO.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
	write.table(peaks_afterCO[[2]],"~/Analyze4C/temp/peaks2_afterCO.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
	
	#wherever there is a contact - intersect with both contact band files
	cat("\nintersecting the contact bands files and the ChIPseq data from the file that you chose\nplease wait...\n\n")
	system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[1]," -wo > ~/Analyze4C/temp/ChIPseqComparison1.bed",sep=""))
	system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[2]," -wo > ~/Analyze4C/temp/ChIPseqComparison2.bed",sep=""))
	ChIPseqComparison1 <- read.table("~/Analyze4C/temp/ChIPseqComparison1.bed")	
	ChIPseqComparison2 <- read.table("~/Analyze4C/temp/ChIPseqComparison2.bed")

	#getting the date and time in order to distinguish between file names of plots
	DandT1 <- toString(Sys.time())
	DandT2 <- gsub(" ","_",DandT1)
	DandT2 <- gsub(":","",DandT2)
	
	#compare the numbers of both - ask user if to compare all, trans, cis , or individual chromosomes
	repeat
	{
		cat("\ncomparing the intersections of the peaks:\n")
		ans3 <- as.numeric(readline(prompt=cat("\nchoose a chromosome that you would like to compare the intersections for:\n\n1) All\n2) trans\n3) specific chromosome (will include the cis chromosome)\n4) exit\n\n")))
		if(ans3 == 1) #all
		{
			#1) sum all peaks and compare
			#2) get a distribution of the peaks
			#3) compare specific contacts and the results of peaks, this is only when both peaks files intersect at the same contact band
			#4) get the average of peaks for each file and compare them

			#the mcnemar test
			#since we have two different 'peaks_afterCO' bed files (one for each removal of sections according to each contact bands file separately) we will perform the mcnemar test twice
			#this is not the most accurate way of doing it (since we the difference between the removals might be significant) but best for now
			ans.mcnemar <- readline("\nwould you like to perform a mcnemar test?\ny/n\n\n")
			if(ans.mcnemar == "y")
			{
				cat("\nmcnemar test 1:\npeaks1_afterCO - \n")
				mcnemar_tester("~/Analyze4C/temp/peaks1_afterCO.bed",paste("~/Analyze4C/contact_bands/",conts_filename[1],sep=""),paste("~/Analyze4C/contact_bands/",conts_filename[2],sep=""))
				cat("\nmcnemar test 2:\npeaks2_afterCO - \n")
				mcnemar_tester("~/Analyze4C/temp/peaks2_afterCO.bed",paste("~/Analyze4C/contact_bands/",conts_filename[1],sep=""),paste("~/Analyze4C/contact_bands/",conts_filename[2],sep=""))
			}

			ans.prop <- readline(prompt=cat("\n\nwould you like to perform a prop test on the intersected bp and the intersected peaks?\ny/n\n\n"))
			if(ans.prop == "y")
			{	
				#BPs:

				#getting the sum of all the bps of the intersections
				intersectionBPsum1 <- sum(ChIPseqComparison1[,8])
				intersectionBPsum2 <- sum(ChIPseqComparison2[,8])

				#getting the sum of all the bps of the contact bands
				all_contsBPsum1 <- sum(conts[[1]][,3]-conts[[1]][,2])
				all_contsBPsum2 <- sum(conts[[2]][,3]-conts[[2]][,2])
				
				#getting the sum of all the bps that are above the peaks cutoff
				all_peaksBPsum1 <- sum(peaks_afterCO[[1]][,3]-peaks_afterCO[[1]][,2])
				all_peaksBPsum2 <- sum(peaks_afterCO[[2]][,3]-peaks_afterCO[[2]][,2])
				
				#peaks:

				#getting the sum of all the peaks of the intersections
				ChIPseqComparison1_all <- ChIPseqComparison1[row.names(unique(ChIPseqComparison1[,c(colnames(ChIPseqComparison1)[1],colnames(ChIPseqComparison1)[2])])),]
				intersectionPeaksSum1 <- sum(ChIPseqComparison1_all[,4])
				ChIPseqComparison2_all <- ChIPseqComparison2[row.names(unique(ChIPseqComparison2[,c(colnames(ChIPseqComparison2)[1],colnames(ChIPseqComparison2)[2])])),]
				intersectionPeaksSum2 <- sum(ChIPseqComparison2_all[,4])
				
				#getting the sum of all the peaks above the peaks cutoff
				all_PeaksSum1 <- sum(peaks_afterCO[[1]][,4])
				all_PeaksSum2 <- sum(peaks_afterCO[[2]][,4])
				
				cat("\nthe peaks file chosen:\n",peaks_filename,"\n")
				cat("\nthe number of bps of peaks above the cutoff, on the peaks data after removed sections according to",conts_filename[1],"is:",all_peaksBPsum1,"\n")
				cat("\nthe number of bps of peaks above the cutoff, on the peaks data after removed sections according to",conts_filename[2],"is:",all_peaksBPsum2,"\n")
				cat("\nthe sum of peaks above the cutoff, on the peaks data after removed sections according to",conts_filename[1],"is:",all_PeaksSum1,"\n")
				cat("\nthe sum of peaks above the cutoff, on the peaks data after removed sections according to",conts_filename[2],"is:",all_PeaksSum2,"\n")
				cat("\nthe number of bps of contact bands on the whole genome of the contact bands file",conts_filename[1],"is:",all_contsBPsum1,"\n")
				cat("\nthe number of bps of contact bands on the whole genome of the contact bands file",conts_filename[2],"is:",all_contsBPsum2,"\n")
				cat("\ncontact bands file - ",conts_filename[1],":\nthe number of bps intersected are ",intersectionBPsum1," which are ",(intersectionBPsum1/all_peaksBPsum1)*100,"% of all the bps of peaks above the cutoff\n\n",sep="")
				cat("\ncontact bands file - ",conts_filename[2],":\nthe number of bps intersected are ",intersectionBPsum2," which are ",(intersectionBPsum2/all_peaksBPsum2)*100,"% of all the bps of peaks above the cutoff\n\n",sep="")
				cat("\ncontact bands file - ",conts_filename[1],":\nthe number of bps intersected are ",intersectionBPsum1," which are ",(intersectionBPsum1/all_contsBPsum1)*100,"% of all the bps of the contact bands in the whole genome\n\n",sep="")
				cat("\ncontact bands file - ",conts_filename[2],":\nthe number of bps intersected are ",intersectionBPsum2," which are ",(intersectionBPsum2/all_contsBPsum2)*100,"% of all the bps of the contact bands in the whole genome\n\n",sep="")				
				cat("\ncontact bands file - ",conts_filename[1],":\nthe sum of peaks intersected are ",intersectionPeaksSum1," which are ",(intersectionPeaksSum1/all_PeaksSum1)*100,"% of all the peaks above the cutoff\n\n",sep="")
				cat("\ncontact bands file - ",conts_filename[2],":\nthe sum of peaks intersected are ",intersectionPeaksSum2," which are ",(intersectionPeaksSum2/all_PeaksSum2)*100,"% of all the peaks above the cutoff\n\n",sep="")			
								
				#creating plot		
				contsName1 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",conts_filename[1],"on the plot:\n"))
				contsName2 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",conts_filename[2],"on the plot:\n"))
				contact_source <- c(contsName1,contsName2)
				conts_base_pairs_percentage <- c((intersectionBPsum1/all_contsBPsum1)*100,(intersectionBPsum2/all_contsBPsum2)*100)
				peaks_base_pairs_percentage <- c((intersectionBPsum1/all_peaksBPsum1)*100,(intersectionBPsum2/all_peaksBPsum2)*100)
				peaks_percentage <- c((intersectionPeaksSum1/all_PeaksSum1)*100,(intersectionPeaksSum2/all_PeaksSum2)*100)
				
				conts_bp_plotData <- data.frame(contact_source,conts_base_pairs_percentage)
				peaks_bp_plotData <- data.frame(contact_source,peaks_base_pairs_percentage)
				peaks_plotData <- data.frame(contact_source,peaks_percentage)

				#creating the percentage of bp intersected from ChIPseq peaks plot:
				plo_peaksBP <- ggplot2::ggplot(peaks_bp_plotData, ggplot2::aes(contact_source,peaks_base_pairs_percentage,fill=contact_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of bp from ChIPseq peaks intersected between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"in the whole genome")) + ggplot2::labs(x="ChIPseq source",y="intersected bp percentage ChIPseq peaks(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				if(intersectionBPsum1>0 & intersectionBPsum2>0)
				{			
					#testing the significance of difference between percentages of both data groups
					cat("\n\ntesting the significance of the difference between the percentages of bp intersected out of the bp of the peaks, for both contact bands groups...\n\n")
					#prop_peaksBP <- prop.test(peaks_bp_plotData$peaks_base_pairs_percentage,c(all_peaksBPsum1,all_peaksBPsum2),correct=FALSE)
					cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the peaks file\nthe 2nd creates a table of all the successes (bp of peaks intersections) and fails (all the bp of the peaks that didn't intersect) and uses them\n\n")					
					prop_peaksBP1 <- prop.test(c(intersectionBPsum1,intersectionBPsum2),c(all_peaksBPsum1,all_peaksBPsum2),correct=FALSE)
					print(prop_peaksBP1)
					prop_peaksBP2 <- prop.test(matrix(c(all_peaksBPsum1-intersectionBPsum1,intersectionBPsum1,all_peaksBPsum2-intersectionBPsum2,intersectionBPsum2),ncol=2),correct=FALSE)
					print(prop_peaksBP2)
					prop_ans2 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
					if(prop_ans2 == 1)
					{
						prop_peaksBP <- prop_peaksBP1
					}
					else
					{
						prop_peaksBP <- prop_peaksBP2
					}
					
					#calculating effect sizes
					cat("\neffect sizes:\n")
					compute.es::propes((intersectionBPsum1/all_peaksBPsum1),(intersectionBPsum2/all_peaksBPsum2),all_peaksBPsum1,all_peaksBPsum2)
					cat("\n")
					
					#phi
					cat("\nphi:\n")
					print(psych::phi(matrix(c(intersectionBPsum1,all_peaksBPsum1-intersectionBPsum1,intersectionBPsum2,all_peaksBPsum2-intersectionBPsum2),ncol=2,byrow=TRUE)))
					cat("\n")
					
					#creating a barplot with significance stars	
					cat("\nplotting the peaks bp percentage data\nplease wait...\n\n")					
					label.df <- data.frame(contact_source = peaks_bp_plotData[which.max(peaks_bp_plotData$peaks_base_pairs_percentage),1],peaks_base_pairs_percentage = peaks_bp_plotData[which.max(peaks_bp_plotData$peaks_base_pairs_percentage),2]+1)
					if(prop_peaksBP$p.value <= 0.0001)
					{
						print(plo_peaksBP + ggplot2::geom_text(data = label.df, label = "****",size=8))
					}
					else if(prop_peaksBP$p.value <= 0.001)
					{
						print(plo_peaksBP + ggplot2::geom_text(data = label.df, label = "***",size=8))
					}
					else if(prop_peaksBP$p.value <= 0.01)
					{
						print(plo_peaksBP + ggplot2::geom_text(data = label.df, label = "**",size=8))
					}
					else if(prop_peaksBP$p.value <= 0.05)
					{
						print(plo_peaksBP + ggplot2::geom_text(data = label.df, label = "*",size=8))
					}
					else
					{
						print(plo_peaksBP)
					}
				}
				else
				{
					if(intersectionBPsum1 == 0)
					{
						cat("\nthe",contsName1,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum2 == 0)
					{
						cat("\nthe",contsName2,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum1 == 0 & intersectionBPsum2 == 0)
					{
						cat("\nthe",contsName1,"and",contsName2,"data do not have any intersections\n\n")
					}			
					cat("\nplotting the peaks bp percentage data\nplease wait...\n\n")	
					print(plo_peaksBP)
				}

				ans5 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans5 == "y")
				{
					nm <- paste("peaksBP_peaksintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for ChIPseqVScontacts_plots.txt:
					struct <- "two contact bands files intersected with a peaks file"
					
					bpORpeaks <- "bp"
					
					percentageOf <- "ChIPseq peaks"

					if(ans4 == 1)
					{
						peaks_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						peaks_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						peaks_CO_appliedTo <- "each chromosome separately"
					}
					
					Intersection_chromosomes <- "whole genome"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to ChIPseqVScontacts_plots.txt
					ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(nm,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
					rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm ChIPseqVScontacts_plots.txt")
					write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
				}
				
				#creating the percentage of bp intersected from contact bands plot:
				plo_contsBP <- ggplot2::ggplot(conts_bp_plotData, ggplot2::aes(contact_source,conts_base_pairs_percentage,fill=contact_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of bp from contact bands intersected between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"in the whole genome")) + ggplot2::labs(x="peaks source",y="intersected bp percentage of contact bands(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				if(intersectionBPsum1>0 & intersectionBPsum2>0)
				{			
					#testing the significance of difference between percentages of both contact bands groups
					cat("\n\ntesting the significance of the difference between the percentages of bp intersected of both contact bands groups...\n")
					#prop_contsBP <- prop.test(conts_bp_plotData$conts_base_pairs_percentage,c(all_contsBPsum1,all_contsBPsum2),correct=FALSE)
					cat("\nthe 1st is the comparison of proportions of intersections out of all bps for both contact bands groups\nthe 2nd creates a table of all the successes (bp intersections) and fails (all the bp of the contact bands that didn't intersect) and uses them\n\n")
					prop_contsBP1 <- prop.test(c(intersectionBPsum1,intersectionBPsum2),c(all_contsBPsum1,all_contsBPsum2),correct=FALSE)
					print(prop_contsBP1)
					prop_contsBP2 <- prop.test(matrix(c(all_contsBPsum1-intersectionBPsum1,intersectionBPsum1,all_contsBPsum2-intersectionBPsum2,intersectionBPsum2),ncol=2),correct=FALSE)
					print(prop_contsBP2)
					prop_ans1 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
					if(prop_ans1 == 1)
					{
						prop_contsBP <- prop_contsBP1
					}
					else
					{
						prop_contsBP <- prop_contsBP2
					}
					
					#calculating effect sizes
					cat("\neffect sizes:\n")
					compute.es::propes((intersectionBPsum1/all_contsBPsum1),(intersectionBPsum2/all_contsBPsum2),all_contsBPsum1,all_contsBPsum2)
					cat("\n")
					
					#phi
					cat("\nphi:\n")
					print(psych::phi(matrix(c(intersectionBPsum1,all_contsBPsum1-intersectionBPsum1,intersectionBPsum2,all_contsBPsum2-intersectionBPsum2),ncol=2,byrow=TRUE)))
					cat("\n")
					
					#creating a barplot with significance stars
					cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")					
					label.df <- data.frame(contact_source = conts_bp_plotData[which.max(conts_bp_plotData$conts_base_pairs_percentage),1],conts_base_pairs_percentage = conts_bp_plotData[which.max(conts_bp_plotData$conts_base_pairs_percentage),2]+1)
					if(prop_contsBP$p.value <= 0.0001)
					{
						print(plo_contsBP + ggplot2::geom_text(data = label.df, label = "****",size=8))
					}
					else if(prop_contsBP$p.value <= 0.001)
					{
						print(plo_contsBP + ggplot2::geom_text(data = label.df, label = "***",size=8))
					}
					else if(prop_contsBP$p.value <= 0.01)
					{
						print(plo_contsBP + ggplot2::geom_text(data = label.df, label = "**",size=8))
					}
					else if(prop_contsBP$p.value <= 0.05)
					{
						print(plo_contsBP + ggplot2::geom_text(data = label.df, label = "*",size=8))
					}
					else
					{
						print(plo_contsBP)
					}
				}
				else
				{
					if(intersectionBPsum1 == 0)
					{
						cat("\nthe",contsName1,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum2 == 0)
					{
						cat("\nthe",contsName2,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum1 == 0 & intersectionBPsum2 == 0)
					{
						cat("\nthe",contsName1,"and",contsName2,"data do not have any intersections\n\n")
					}			
					cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")	
					print(plo_contsBP)
				}

				ans8 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans8 == "y")
				{
					nm <- paste("ContactsBP_peaksintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for ChIPseqVScontacts_plots.txt:
					struct <- "two contact bands files intersected with an peaks file"

					bpORpeaks <- "bp"
					
					percentageOf <- "Contact bands"
					
					if(ans4 == 1)
					{
						peaks_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						peaks_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						peaks_CO_appliedTo <- "each chromosome separately"
					}
					
					Intersection_chromosomes <- "whole genome"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to ChIPseqVScontacts_plots.txt
					ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(nm,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
					rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm ChIPseqVScontacts_plots.txt")
					write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
				}
								
				#creating the percentage of peaks value (p-score or number of tags) intersected of sum of the same peaks value plot:
				plopeaks <- ggplot2::ggplot(peaks_plotData, ggplot2::aes(contact_source,peaks_percentage,fill=contact_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected peaks between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"in the whole genome")) + ggplot2::labs(x="peaks source",y="intersected peaks percentage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				if(intersectionPeaksSum1>0 & intersectionPeaksSum2>0)
				{
					#testing the significance of difference between percentages of both data groups
					cat("\n\ntesting the significance of the difference between the percentages of peaks intersected of both contact bands groups...\n\n")
					#proppeaks <- prop.test(peaks_plotData$peaks_percentage,c(all_PeaksSum1,all_PeaksSum2),correct=FALSE)
					cat("\nthe 1st is the comparison of proportions of intersections out of all peaks values in the peaks\nthe 2nd creates a table of all the successes (peaks intersections) and fails (all the peaks of the peaks that didn't intersect) and uses them\n\n")					
					proppeaks1 <- prop.test(c(intersectionPeaksSum1,intersectionPeaksSum2),c(all_PeaksSum1,all_PeaksSum2),correct=FALSE)
					print(proppeaks1)
					proppeaks2 <- prop.test(matrix(c(all_PeaksSum1-intersectionPeaksSum1,intersectionPeaksSum1,all_PeaksSum2-intersectionPeaksSum2,intersectionPeaksSum2),ncol=2),correct=FALSE)
					print(proppeaks2)					
					prop_ans3 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
					if(prop_ans3 == 1)
					{
						proppeaks <- proppeaks1
					}
					else
					{
						proppeaks <- proppeaks2
					}
					
					#calculating effect sizes
					cat("\neffect sizes:\n")
					compute.es::propes((intersectionPeaksSum1/all_PeaksSum1),(intersectionPeaksSum2/all_PeaksSum2),all_PeaksSum1,all_PeaksSum2)
					cat("\n")
					
					#phi
					cat("\nphi:\n")
					print(psych::phi(matrix(c(intersectionPeaksSum1,all_PeaksSum1-intersectionPeaksSum1,intersectionPeaksSum2,all_PeaksSum2-intersectionPeaksSum2),ncol=2,byrow=TRUE)))
					cat("\n")
					
					#creating a barplot with significance stars
					cat("\nplotting the peaks percentage data\nplease wait...\n\n")					
					label.df <- data.frame(contact_source = peaks_plotData[which.max(peaks_plotData$peaks_percentage),1],peaks_percentage = peaks_plotData[which.max(peaks_plotData$peaks_percentage),2]+1)	
					if(proppeaks$p.value <= 0.0001)
					{
						print(plopeaks + ggplot2::geom_text(data = label.df, label = "****",size=8))
					}
					else if(proppeaks$p.value <= 0.001)
					{
						print(plopeaks + ggplot2::geom_text(data = label.df, label = "***",size=8))
					}
					else if(proppeaks$p.value <= 0.01)
					{
						print(plopeaks + ggplot2::geom_text(data = label.df, label = "**",size=8))
					}
					else if(proppeaks$p.value <= 0.05)
					{
						print(plopeaks + ggplot2::geom_text(data = label.df, label = "*",size=8))
					}
					else
					{
						print(plopeaks)
					}
				}	
				else
				{
					if(intersectionPeaksSum1 == 0)
					{
						cat("\nthe",contsName1,"data does not have any intersections\n\n")
					}
					else if(intersectionPeaksSum2 == 0)
					{
						cat("\nthe",contsName2,"data does not have any intersections\n\n")
					}
					else if(intersectionPeaksSum1 == 0 & intersectionPeaksSum2 == 0)
					{
						cat("\nthe",contsName1,"and",contsName2,"data do not have any intersections\n\n")
					}				
					cat("\nplotting the peaks percentage data\nplease wait...\n\n")
					print(plopeaks)
				}
				
				ans6 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans6 == "y")
				{
					nm <- paste("peaks_peaksintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for ChIPseqVScontacts_plots.txt:
					struct <- "two contact bands files intersected with an peaks file"
					
					bpORpeaks <- "peaks value"
					
					percentageOf <- "peaks value"
					
					if(ans4 == 1)
					{
						peaks_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						peaks_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						peaks_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "whole genome"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to ChIPseqVScontacts_plots.txt
					ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(nm,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
					rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm ChIPseqVScontacts_plots.txt")
					write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}
			}
			
			#venn diagrams
			if(readline(prompt=cat("\nwould you like to create venn diagrams?\ny/n\n\n")) == "y")
			{
				#creating the venn data frame for 0 percent
				venn_DF <- data.frame(rep(1,nrow(peaks_afterCO[[1]])),matrix(0,nrow(peaks_afterCO[[1]]),2))##!

				cat("\nvenn diagram parameters:\n\n")
				venn_ans1 <- readline(prompt=cat("\nwould you like to consider overlaps to be only those that intersect more than 1 bp?\ny/n\n\n"))
				if(venn_ans1 == "y")
				{
					ovlp_cov <- as.numeric(readline(prompt=cat("\nenter the minimum overlap coverage that will be considered an intersection? (between 0 and 1)\n\n")))
					
					venn_ans2 <- readline(prompt=cat("\nyou have chosen to accept an intersection only if there is at least",ovlp_cov,"overlap.\nwould you consider accepting an intersection if there is more than a certain amount of reads that intersect (e.g. if you there are more than an x amount of reads intersected, even though each one doesn't overlap more than",ovlp_cov,", altogether they are sufficient enough to be considered as an overlap)\ny/n\n\n"))
					if(venn_ans2 == "y")
					{
						ovlp_num <- as.integer(readline(prompt=cat("\nhow many reads (minimum) intersected will be considered as an overlap?\n\n")))
					}
				}
				
				#filling in the data for the venn diagrams
				if(venn_ans1 == "y")
				{
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[1]," -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[2]," -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A2.bed",sep=""))
					venn_A1 <- read.table(paste("~/Analyze4C/temp/venn_A1.bed",sep=""))
					venn_A2 <- read.table(paste("~/Analyze4C/temp/venn_A2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn_A1.bed")
					system("rm ~/Analyze4C/temp/venn_A2.bed")
					if(venn_ans2 == "y")
					{
						system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[1]," -c > ~/Analyze4C/temp/venn_B1.bed",sep=""))
						system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[2]," -c > ~/Analyze4C/temp/venn_B2.bed",sep=""))
						venn_B1 <- read.table(paste("~/Analyze4C/temp/venn_B1.bed",sep=""))
						venn_B2 <- read.table(paste("~/Analyze4C/temp/venn_B2.bed",sep=""))
						#removing files
						system("rm ~/Analyze4C/temp/venn_B1.bed")
						system("rm ~/Analyze4C/temp/venn_B2.bed")
						#if the number in the 1st is 0 but the number in the 2nd is more than ovlp_num then we will consider it to be 1
						venn_A1[venn_B1[,5]>=ovlp_num,5] <- 1
						venn_A2[venn_B2[,5]>=ovlp_num,5] <- 1
					}
					#add the correct number to each column
					venn_DF[venn_A1[,5]>=1,2] <- 1
					venn_DF[venn_A2[,5]>=1,3] <- 1
				}
				else
				{
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[1]," -c > ~/Analyze4C/temp/venn1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[2]," -c > ~/Analyze4C/temp/venn2.bed",sep=""))
					venn1 <- read.table(paste("~/Analyze4C/temp/venn1.bed",sep=""))
					venn2 <- read.table(paste("~/Analyze4C/temp/venn2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn1.bed")
					system("rm ~/Analyze4C/temp/venn2.bed")
					#add the correct number to each column
					venn_DF[venn1[,5]>=1,2] <- 1
					venn_DF[venn2[,5]>=1,3] <- 1				
				}
				
				# create venn diagram using the data frame venn_DF
				venn_flag <- 0 #if venn_flag is 1 then we save a diagram and it should be recorded in 'ChIPseqVScontacts_plots'
				vennName <- paste("venn_peaksVS2psConts_",DandT2,".jpg",sep="")
				venn_flag <- venn_creator(venn_DF,2,1,vennName)			

				#recording the data into 'ChIPseqVScontacts_plots.txt'
				if(venn_flag == 1)
				{					
					#getting parameters for ChIPseqVScontacts_plots.txt:
					struct <- "two contact bands files intersected with an peaks file - Venn Diagram"					
					bpORpeaks <- NA
					percentageOf <- NA
					if(ans4 == 1)
					{
						peaks_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						peaks_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						peaks_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "whole genome"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to ChIPseqVScontacts_plots.txt
					ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(vennName,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
					rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm ChIPseqVScontacts_plots.txt")
					write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}	
			}											
		}
		else if(ans3 == 2) #trans
		{
			#if the cis chromosomes of the contact bands chosen don't match, we can't remove both cis chromosomes from the peaks file, so trans can't be done
			if(cis[1]!=cis[2])
			{
				cat("\the cis chromosomes of both contact bands files chosen are different\na the intersections of only trans in this case cannot be done\nin order to do so, you must choose two contact bands with the same cis chromosome\n\n")
			}
			else
			{
				#the mcnemar test
				#since we have two different 'peaks_afterCO' bed files (one for each removal of sections according to each contact bands file separately) we will perform the mcnemar test twice
				#this is not the most accurate way of doing it (since we the difference between the removals might be significant) but best for now
				ans.mcnemar <- readline("\nwould you like to perform a mcnemar test?\ny/n\n\n")
				if(ans.mcnemar == "y")
				{
					peaks1_afterCO_mcnemar <- peaks_afterCO[[1]][peaks_afterCO[[1]][,1]!=cis[1],]
					peaks2_afterCO_mcnemar <- peaks_afterCO[[2]][peaks_afterCO[[2]][,1]!=cis[2],]
					
					#i need to test the conts and see how to extract the correct data since its a list
					#the reason why both cis[1] and cis[2] are written is just for comfort and aesthetic reasons
					conts1_mcnemar <- conts[[1]][conts[[1]][,1]!=cis[1],]
					conts2_mcnemar <- conts[[2]][conts[[2]][,1]!=cis[2],]
					write.table(peaks1_afterCO_mcnemar,"~/Analyze4C/temp/peaks1_afterCO_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
					write.table(peaks2_afterCO_mcnemar,"~/Analyze4C/temp/peaks2_afterCO_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
					write.table(conts1_mcnemar,"~/Analyze4C/temp/conts1_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
					write.table(conts2_mcnemar,"~/Analyze4C/temp/conts2_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

					cat("\nmcnemar test 1:\npeaks1_afterCO - \n")
					mcnemar_tester("~/Analyze4C/temp/peaks1_afterCO_mcnemar.bed","~/Analyze4C/temp/conts1_mcnemar.bed","~/Analyze4C/temp/conts2_mcnemar.bed")
					cat("\nmcnemar test 2:\npeaks2_afterCO - \n")
					mcnemar_tester("~/Analyze4C/temp/peaks2_afterCO_mcnemar.bed","~/Analyze4C/temp/conts1_mcnemar.bed","~/Analyze4C/temp/conts2_mcnemar.bed")
					
					system("rm ~/Analyze4C/temp/peaks1_afterCO_mcnemar.bed")
					system("rm ~/Analyze4C/temp/peaks2_afterCO_mcnemar.bed")
					system("rm ~/Analyze4C/temp/conts1_mcnemar.bed")
					system("rm ~/Analyze4C/temp/conts2_mcnemar.bed")													
				}
			
				ans.prop <- readline(prompt=cat("\n\nwould you like to perform a prop test on the intersected bp and the intersected peaks?\ny/n\n\n"))
				if(ans.prop == "y")
				{
					#BPs:
					
					#getting the sum of all the bps of the intersections in trans
					intersectionBPsum1 <- sum(ChIPseqComparison1[ChIPseqComparison1[,1]!=cis[1],8])
					intersectionBPsum2 <- sum(ChIPseqComparison2[ChIPseqComparison2[,1]!=cis[2],8])
					
					#getting the sum of all the bps of the contact bands in trans
					
					#here too i need to see how im getting the data from conts @
					all_contsBPsum1 <- sum(conts[[1]][conts[[1]][,1]!=cis[1],3]-conts[[1]][conts[[1]][,1]!=cis[1],2])
					all_contsBPsum2 <- sum(conts[[2]][conts[[2]][,1]!=cis[2],3]-conts[[2]][conts[[2]][,1]!=cis[2],2])
					
					#getting the sum of all the bps that are above the peaks cutoff in trans
					all_peaksBPsum1 <- sum(peaks_afterCO[[1]][peaks_afterCO[[1]][,1]!=cis[1],3]-peaks_afterCO[[1]][peaks_afterCO[[1]][,1]!=cis[1],2])
					all_peaksBPsum2 <- sum(peaks_afterCO[[2]][peaks_afterCO[[2]][,1]!=cis[2],3]-peaks_afterCO[[2]][peaks_afterCO[[2]][,1]!=cis[2],2])
					
					#peaks:

					#getting the sum of all the peaks of the intersections
					ChIPseqComparison1_trans <- ChIPseqComparison1[row.names(unique(ChIPseqComparison1[,c(colnames(ChIPseqComparison1)[1],colnames(ChIPseqComparison1)[2])])),]
					intersectionPeaksSum1 <- sum(ChIPseqComparison1_trans[ChIPseqComparison1_trans[,1]!=cis[1],4])
					ChIPseqComparison2_trans <- ChIPseqComparison2[row.names(unique(ChIPseqComparison2[,c(colnames(ChIPseqComparison2)[1],colnames(ChIPseqComparison2)[2])])),]
					intersectionPeaksSum2 <- sum(ChIPseqComparison2_trans[ChIPseqComparison2_trans[,1]!=cis[2],4])
					
					#getting the sum of all the peaks above the peaks cutoff in trans
					all_PeaksSum1 <- sum(peaks_afterCO[[1]][peaks_afterCO[[1]][,1]!=cis[1],4])
					all_PeaksSum2 <- sum(peaks_afterCO[[2]][peaks_afterCO[[2]][,1]!=cis[2],4])

					cat("\nthe peaks file chosen:\n",peaks_filename,"\n")
					cat("\nthe number of bps of peaks above the cutoff in trans, on the peaks data after removed sections according to",conts_filename[1],"is:",all_peaksBPsum1,"\n")
					cat("\nthe number of bps of peaks above the cutoff in trans, on the peaks data after removed sections according to",conts_filename[2],"is:",all_peaksBPsum2,"\n")
					cat("\nthe sum of peaks above the cutoff in trans, on the peaks data after removed sections according to",conts_filename[1],"is:",all_PeaksSum1,"\n")
					cat("\nthe sum of peaks above the cutoff in trans, on the peaks data after removed sections according to",conts_filename[2],"is:",all_PeaksSum2,"\n")
					cat("\nthe number of bps of contact bands in the trans of the contact bands file",conts_filename[1],"is:",all_contsBPsum1,"\n")
					cat("\nthe number of bps of contact bands in the trans of the contact bands file",conts_filename[2],"is:",all_contsBPsum2,"\n")
					cat("\ncontact bands file - ",conts_filename[1],":\nthe number of bps intersected in trans are ",intersectionBPsum1," which are ",(intersectionBPsum1/all_peaksBPsum1)*100,"% of all the bps of peaks above the cutoff in trans\n\n",sep="")
					cat("\ncontact bands file - ",conts_filename[2],":\nthe number of bps intersected in trans are ",intersectionBPsum2," which are ",(intersectionBPsum2/all_peaksBPsum2)*100,"% of all the bps of peaks above the cutoff in trans\n\n",sep="")
					cat("\ncontact bands file - ",conts_filename[1],":\nthe number of bps intersected in trans are ",intersectionBPsum1," which are ",(intersectionBPsum1/all_contsBPsum1)*100,"% of all the bps of the contact bands in trans\n\n",sep="")
					cat("\ncontact bands file - ",conts_filename[2],":\nthe number of bps intersected in trans are ",intersectionBPsum2," which are ",(intersectionBPsum2/all_contsBPsum2)*100,"% of all the bps of the contact bands in trans\n\n",sep="")				
					cat("\ncontact bands file - ",conts_filename[1],":\nthe sum of peaks intersected in trans are ",intersectionPeaksSum1," which are ",(intersectionPeaksSum1/all_PeaksSum1)*100,"% of all the peaks above the cutoff in trans\n\n",sep="")
					cat("\ncontact bands file - ",conts_filename[2],":\nthe sum of peaks intersected in trans are ",intersectionPeaksSum2," which are ",(intersectionPeaksSum2/all_PeaksSum2)*100,"% of all the peaks above the cutoff in trans\n\n",sep="")			
									
					#creating plot		
					contsName1 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",conts_filename[1],"on the plot:\n"))
					contsName2 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",conts_filename[2],"on the plot:\n"))
					contact_source <- c(contsName1,contsName2)
					conts_base_pairs_percentage <- c((intersectionBPsum1/all_contsBPsum1)*100,(intersectionBPsum2/all_contsBPsum2)*100)
					peaks_base_pairs_percentage <- c((intersectionBPsum1/all_peaksBPsum1)*100,(intersectionBPsum2/all_peaksBPsum2)*100)
					peaks_percentage <- c((intersectionPeaksSum1/all_PeaksSum1)*100,(intersectionPeaksSum2/all_PeaksSum2)*100)
					
					conts_bp_plotData <- data.frame(contact_source,conts_base_pairs_percentage)
					peaks_bp_plotData <- data.frame(contact_source,peaks_base_pairs_percentage)
					peaks_plotData <- data.frame(contact_source,peaks_percentage)

					#creating the percentage of bp from ChIPseq peaks plot:
					plo_peaksBP <- ggplot2::ggplot(peaks_bp_plotData, ggplot2::aes(contact_source,peaks_base_pairs_percentage,fill=contact_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of bp from ChIPseq peaks intersected between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"in trans")) + ggplot2::labs(x="peaks source",y="intersected bp percentage ChIPseq peaks(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
					if(intersectionBPsum1>0 & intersectionBPsum2>0)
					{				
						#testing the significance of difference between percentages of both data groups
						cat("\n\ntesting the significance of the difference between the percentages of bp intersected out of the bp of the peaks, for both contact bands groups...\n\n")
						#prop_peaksBP <- prop.test(peaks_bp_plotData$peaks_base_pairs_percentage,c(all_peaksBPsum1,all_peaksBPsum2),correct=FALSE)
						cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the peaks file\nthe 2nd creates a table of all the successes (bp of peaks intersections) and fails (all the bp of the peaks that didn't intersect) and uses them\n\n")					
						prop_peaksBP1 <- prop.test(c(intersectionBPsum1,intersectionBPsum2),c(all_peaksBPsum1,all_peaksBPsum2),correct=FALSE)
						print(prop_peaksBP1)
						prop_peaksBP2 <- prop.test(matrix(c(all_peaksBPsum1-intersectionBPsum1,intersectionBPsum1,all_peaksBPsum2-intersectionBPsum2,intersectionBPsum2),ncol=2),correct=FALSE)
						print(prop_peaksBP2)
						prop_ans2 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
						if(prop_ans2 == 1)
						{
							prop_peaksBP <- prop_peaksBP1
						}
						else
						{
							prop_peaksBP <- prop_peaksBP2
						}
						
						#calculating effect sizes
						cat("\neffect sizes:\n")
						compute.es::propes((intersectionBPsum1/all_peaksBPsum1),(intersectionBPsum2/all_peaksBPsum2),all_peaksBPsum1,all_peaksBPsum2)
						cat("\n")
						
						#phi
						cat("\nphi:\n")
						print(psych::phi(matrix(c(intersectionBPsum1,all_peaksBPsum1-intersectionBPsum1,intersectionBPsum2,all_peaksBPsum2-intersectionBPsum2),ncol=2,byrow=TRUE)))
						cat("\n")
						
						#creating a barplot with significance stars
						cat("\nplotting the peaks bp percentage data\nplease wait...\n\n")					
						label.df <- data.frame(contact_source = peaks_bp_plotData[which.max(peaks_bp_plotData$peaks_base_pairs_percentage),1],peaks_base_pairs_percentage = peaks_bp_plotData[which.max(peaks_bp_plotData$peaks_base_pairs_percentage),2]+1)
						if(prop_peaksBP$p.value <= 0.0001)
						{
							print(plo_peaksBP + ggplot2::geom_text(data = label.df, label = "****",size=8))
						}
						else if(prop_peaksBP$p.value <= 0.001)
						{
							print(plo_peaksBP + ggplot2::geom_text(data = label.df, label = "***",size=8))
						}
						else if(prop_peaksBP$p.value <= 0.01)
						{
							print(plo_peaksBP + ggplot2::geom_text(data = label.df, label = "**",size=8))
						}
						else if(prop_peaksBP$p.value <= 0.05)
						{
							print(plo_peaksBP + ggplot2::geom_text(data = label.df, label = "*",size=8))
						}
						else
						{
							print(plo_peaksBP)
						}
					}
					else
					{
						if(intersectionBPsum1 == 0)
						{
							cat("\nthe",contsName1,"data does not have any intersections\n\n")
						}
						else if(intersectionBPsum2 == 0)
						{
							cat("\nthe",contsName2,"data does not have any intersections\n\n")
						}
						else if(intersectionBPsum1 == 0 & intersectionBPsum2 == 0)
						{
							cat("\nthe",contsName1,"and",contsName2,"data do not have any intersections\n\n")
						}			
						cat("\nplotting the peaks bp percentage data\nplease wait...\n\n")
						print(plo_peaksBP)
					}

					ans5 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
					if(ans5 == "y")
					{
						nm <- paste("peaksBP_peaksintersectedContacts_",DandT2,".png",sep="")
						wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
						ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
						ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
						
						#getting parameters for ChIPseqVScontacts_plots.txt:
						struct <- "two contact bands files intersected with an peaks file"
						
						bpORpeaks <- "bp"
						
						percentageOf <- "ChIPseq peaks"

						if(ans4 == 1)
						{
							peaks_CO_appliedTo <- "whole genome"
						}
						else if(ans4 == 2)
						{
							peaks_CO_appliedTo <- "trans"
						}
						else if(ans4 == 3)
						{
							peaks_CO_appliedTo <- "each chromosome separately"
						}
						
						Intersection_chromosomes <- "trans"
						
						ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
						if(ans7 == "y")
						{
							notes <- readline(prompt=cat("\nenter the notes:\n\n"))
						}
						else
						{
							notes <- NA
						}
						
						#saving the parameters and details to ChIPseqVScontacts_plots.txt
						ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(nm,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
						#sorting the list of experiments by bait alphabetically (and sorting the row indices)
						ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
						rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
						#adding the new data to the file (by erasing the old one and creating a new one)
						system("rm ChIPseqVScontacts_plots.txt")
						write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
					}
					
					#creating the percentage of bp intersected from contact bands plot:
					plo_contsBP <- ggplot2::ggplot(conts_bp_plotData, ggplot2::aes(contact_source,conts_base_pairs_percentage,fill=contact_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of bp from contact bands intersected between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"in trans")) + ggplot2::labs(x="peaks source",y="intersected bp percentage of contact bands(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
					if(intersectionBPsum1>0 & intersectionBPsum2>0)
					{		
						#testing the significance of difference between percentages of both data groups
						cat("\n\ntesting the significance of the difference between the percentages of bp intersected of both contact bands groups...\n")
						#prop_contsBP <- prop.test(conts_bp_plotData$conts_base_pairs_percentage,c(all_contsBPsum1,all_contsBPsum2),correct=FALSE)
						cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the contact bands\nthe 2nd creates a table of all the successes (bp intersections) and fails (all the bp of the contact bands that didn't intersect) and uses them\n\n")
						prop_contsBP1 <- prop.test(c(intersectionBPsum1,intersectionBPsum2),c(all_contsBPsum1,all_contsBPsum2),correct=FALSE)
						print(prop_contsBP1)
						prop_contsBP2 <- prop.test(matrix(c(all_contsBPsum1-intersectionBPsum1,intersectionBPsum1,all_contsBPsum2-intersectionBPsum2,intersectionBPsum2),ncol=2),correct=FALSE)
						print(prop_contsBP2)
						prop_ans <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
						if(prop_ans == 1)
						{
							prop_contsBP <- prop_contsBP1
						}
						else
						{
							prop_contsBP <- prop_contsBP2
						}
						
						#calculating effect sizes
						cat("\neffect sizes:\n")
						compute.es::propes((intersectionBPsum1/all_contsBPsum1),(intersectionBPsum2/all_contsBPsum2),all_contsBPsum1,all_contsBPsum2)
						cat("\n")
						
						#phi
						cat("\nphi:\n")
						print(psych::phi(matrix(c(intersectionBPsum1,all_contsBPsum1-intersectionBPsum1,intersectionBPsum2,all_contsBPsum2-intersectionBPsum2),ncol=2,byrow=TRUE)))
						cat("\n")
						
						#creating a barplot with significance stars
						cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")					
						label.df <- data.frame(contact_source = conts_bp_plotData[which.max(conts_bp_plotData$conts_base_pairs_percentage),1],conts_base_pairs_percentage = conts_bp_plotData[which.max(conts_bp_plotData$conts_base_pairs_percentage),2]+1)					
						if(prop_contsBP$p.value <= 0.0001)
						{
							print(plo_contsBP + ggplot2::geom_text(data = label.df, label = "****",size=8))
						}
						else if(prop_contsBP$p.value <= 0.001)
						{
							print(plo_contsBP + ggplot2::geom_text(data = label.df, label = "***",size=8))
						}
						else if(prop_contsBP$p.value <= 0.01)
						{
							print(plo_contsBP + ggplot2::geom_text(data = label.df, label = "**",size=8))
						}
						else if(prop_contsBP$p.value <= 0.05)
						{
							print(plo_contsBP + ggplot2::geom_text(data = label.df, label = "*",size=8))
						}
						else
						{
							print(plo_contsBP)
						}
					}
					else
					{
						if(intersectionBPsum1 == 0)
						{
							cat("\nthe",contsName1,"data does not have any intersections\n\n")
						}
						else if(intersectionBPsum2 == 0)
						{
							cat("\nthe",contsName2,"data does not have any intersections\n\n")
						}
						else if(intersectionBPsum1 == 0 & intersectionBPsum2 == 0)
						{
							cat("\nthe",contsName1,"and",contsName2,"data do not have any intersections\n\n")
						}			
						cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")	
						print(plo_contsBP)
					}

					ans8 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
					if(ans8 == "y")
					{
						nm <- paste("ContactsBP_peaksintersectedContacts_",DandT2,".png",sep="")
						wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
						ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
						ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
						
						#getting parameters for ChIPseqVScontacts_plots.txt:
						struct <- "two contact bands files intersected with an peaks file"
						
						bpORpeaks <- "bp"
						
						percentageOf <- "Contact bands"
						
						if(ans4 == 1)
						{
							peaks_CO_appliedTo <- "whole genome"
						}
						else if(ans4 == 2)
						{
							peaks_CO_appliedTo <- "trans"
						}
						else if(ans4 == 3)
						{
							peaks_CO_appliedTo <- "each chromosome separately"
						}
						
						Intersection_chromosomes <- "trans"
						
						ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
						if(ans7 == "y")
						{
							notes <- readline(prompt=cat("\nenter the notes:\n\n"))
						}
						else
						{
							notes <- NA
						}
						
						#saving the parameters and details to ChIPseqVScontacts_plots.txt
						ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(nm,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
						#sorting the list of experiments by bait alphabetically (and sorting the row indices)
						ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
						rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
						#adding the new data to the file (by erasing the old one and creating a new one)
						system("rm ChIPseqVScontacts_plots.txt")
						write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
					}
										
					#creating the percentage of peaks value (p-score or number of tags) intersected of sum of the same peaks value plot:
					plopeaks <- ggplot2::ggplot(peaks_plotData, ggplot2::aes(contact_source,peaks_percentage,fill=contact_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected peaks between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"in trans")) + ggplot2::labs(x="peaks source",y="intersected peaks percentage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
					if(intersectionPeaksSum1>0 & intersectionPeaksSum2>0)
					{
						#testing the significance of difference between percentages of both data groups
						cat("\n\ntesting the significance of the difference between the percentages of peaks intersected of both contat bands groups...\n\n")
						#proppeaks <- prop.test(peaks_plotData$peaks_percentage,c(all_PeaksSum1,all_PeaksSum2),correct=FALSE)
						cat("\nthe 1st is the comparison of proportions of intersections out of all peaks values in the peaks\nthe 2nd creates a table of all the successes (peaks intersections) and fails (all the peaks of the peaks that didn't intersect) and uses them\n\n")					
						proppeaks1 <- prop.test(c(intersectionPeaksSum1,intersectionPeaksSum2),c(all_PeaksSum1,all_PeaksSum2),correct=FALSE)
						print(proppeaks1)
						proppeaks2 <- prop.test(matrix(c(all_PeaksSum1-intersectionPeaksSum1,intersectionPeaksSum1,all_PeaksSum2-intersectionPeaksSum2,intersectionPeaksSum2),ncol=2),correct=FALSE)
						print(proppeaks2)					
						prop_ans3 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
						if(prop_ans3 == 1)
						{
							proppeaks <- proppeaks1
						}
						else
						{
							proppeaks <- proppeaks2
						}
						
						#calculating effect sizes
						cat("\neffect sizes:\n")
						compute.es::propes((intersectionPeaksSum1/all_PeaksSum1),(intersectionPeaksSum2/all_PeaksSum2),all_PeaksSum1,all_PeaksSum2)
						cat("\n")
						
						#phi
						cat("\nphi:\n")
						print(psych::phi(matrix(c(intersectionPeaksSum1,all_PeaksSum1-intersectionPeaksSum1,intersectionPeaksSum2,all_PeaksSum2-intersectionPeaksSum2),ncol=2,byrow=TRUE)))
						cat("\n")
						
						#creating a barplot with significance stars					
						cat("\nplotting the peaks percentage data\nplease wait...\n\n")
						label.df <- data.frame(contact_source = peaks_plotData[which.max(peaks_plotData$peaks_percentage),1],peaks_percentage = peaks_plotData[which.max(peaks_plotData$peaks_percentage),2]+1)
						if(proppeaks$p.value <= 0.0001)
						{
							print(plopeaks + ggplot2::geom_text(data = label.df, label = "****",size=8))
						}
						else if(proppeaks$p.value <= 0.001)
						{
							print(plopeaks + ggplot2::geom_text(data = label.df, label = "***",size=8))
						}
						else if(proppeaks$p.value <= 0.01)
						{
							print(plopeaks + ggplot2::geom_text(data = label.df, label = "**",size=8))
						}
						else if(proppeaks$p.value <= 0.05)
						{
							print(plopeaks + ggplot2::geom_text(data = label.df, label = "*",size=8))
						}
						else
						{
							print(plopeaks)
						}
					}	
					else
					{
						if(intersectionPeaksSum1 == 0)
						{
							cat("\nthe",contsName1,"data does not have any intersections\n\n")
						}
						else if(intersectionPeaksSum2 == 0)
						{
							cat("\nthe",contsName2,"data does not have any intersections\n\n")
						}
						else if(intersectionPeaksSum1 == 0 & intersectionPeaksSum2 == 0)
						{
							cat("\nthe",contsName1,"and",contsName2,"data do not have any intersections\n\n")
						}				
						cat("\nplotting the peaks percentage data\nplease wait...\n\n")
						print(plopeaks)
					}
					
					ans6 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
					if(ans6 == "y")
					{
						nm <- paste("peaks_peaksintersectedContacts_",DandT2,".png",sep="")
						wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
						ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
						ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
						
						#getting parameters for ChIPseqVScontacts_plots.txt:
						struct <- "two contact bands files intersected with an peaks file"
						
						bpORpeaks <- "peaks value"
						
						percentageOf <- "peaks value"
						
						if(ans4 == 1)
						{
							peaks_CO_appliedTo <- "whole genome"
						}
						else if(ans4 == 2)
						{
							peaks_CO_appliedTo <- "trans"
						}
						else if(ans4 == 3)
						{
							peaks_CO_appliedTo <- "each chromosome separately"
						}

						Intersection_chromosomes <- "trans"
						
						ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
						if(ans7 == "y")
						{
							notes <- readline(prompt=cat("\nenter the notes:\n\n"))
						}
						else
						{
							notes <- NA
						}
						
						#saving the parameters and details to ChIPseqVScontacts_plots.txt
						ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(nm,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
						#sorting the list of experiments by bait alphabetically (and sorting the row indices)
						ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
						rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
						#adding the new data to the file (by erasing the old one and creating a new one)
						system("rm ChIPseqVScontacts_plots.txt")
						write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
					}
				}
			}

			#venn diagrams
			if(readline(prompt=cat("\nwould you like to create venn diagrams?\ny/n\n\n")) == "y")
			{
				#creating the venn data frame for 0 percent
				venn_DF <- data.frame(rep(1,nrow(peaks_afterCO[[1]])),matrix(0,nrow(peaks_afterCO[[1]]),2))##!

				cat("\nvenn diagram parameters:\n\n")
				venn_ans1 <- readline(prompt=cat("\nwould you like to consider overlaps to be only those that intersect more than 1 bp?\ny/n\n\n"))
				if(venn_ans1 == "y")
				{
					ovlp_cov <- as.numeric(readline(prompt=cat("\nenter the minimum overlap coverage that will be considered an intersection? (between 0 and 1)\n\n")))
					
					venn_ans2 <- readline(prompt=cat("\nyou have chosen to accept an intersection only if there is at least",ovlp_cov,"overlap.\nwould you consider accepting an intersection if there is more than a certain amount of reads that intersect (e.g. if you there are more than an x amount of reads intersected, even though each one doesn't overlap more than",ovlp_cov,", altogether they are sufficient enough to be considered as an overlap)\ny/n\n\n"))
					if(venn_ans2 == "y")
					{
						ovlp_num <- as.integer(readline(prompt=cat("\nhow many reads (minimum) intersected will be considered as an overlap?\n\n")))
					}
				}
				
				#filling in the data for the venn diagrams
				if(venn_ans1 == "y")
				{
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[1]," -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[2]," -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A2.bed",sep=""))
					venn_A1 <- read.table(paste("~/Analyze4C/temp/venn_A1.bed",sep=""))
					venn_A2 <- read.table(paste("~/Analyze4C/temp/venn_A2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn_A1.bed")
					system("rm ~/Analyze4C/temp/venn_A2.bed")
					if(venn_ans2 == "y")
					{
						system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[1]," -c > ~/Analyze4C/temp/venn_B1.bed",sep=""))
						system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[2]," -c > ~/Analyze4C/temp/venn_B2.bed",sep=""))
						venn_B1 <- read.table(paste("~/Analyze4C/temp/venn_B1.bed",sep=""))
						venn_B2 <- read.table(paste("~/Analyze4C/temp/venn_B2.bed",sep=""))
						#removing files
						system("rm ~/Analyze4C/temp/venn_B1.bed")
						system("rm ~/Analyze4C/temp/venn_B2.bed")
						#if the number in the 1st is 0 but the number in the 2nd is more than ovlp_num then we will consider it to be 1
						venn_A1[venn_B1[,5]>=ovlp_num,5] <- 1
						venn_A2[venn_B2[,5]>=ovlp_num,5] <- 1
					}
					#add the correct number to each column
					venn_DF[venn_A1[,5]>=1 & venn_A1[,1]!=cis[1],2] <- 1
					venn_DF[venn_A2[,5]>=1 & venn_A2[,1]!=cis[2],3] <- 1					
				}
				else
				{
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[1]," -c > ~/Analyze4C/temp/venn1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[2]," -c > ~/Analyze4C/temp/venn2.bed",sep=""))
					venn1 <- read.table(paste("~/Analyze4C/temp/venn1.bed",sep=""))
					venn2 <- read.table(paste("~/Analyze4C/temp/venn2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn1.bed")
					system("rm ~/Analyze4C/temp/venn2.bed")
					#add the correct number to each column
					venn_DF[venn1[,5]>=1 & venn1[,1]!=cis[1],2] <- 1
					venn_DF[venn2[,5]>=1 & venn2[,1]!=cis[2],3] <- 1					
				}
				
				# create venn diagram using the data frame venn_DF
				venn_flag <- 0 #if venn_flag is 1 then we save a diagram and it should be recorded in 'ChIPseqVScontacts_plots'
				vennName <- paste("venn_peaksVS2psConts_",DandT2,".jpg",sep="")
				venn_flag <- venn_creator(venn_DF,2,1,vennName)			

				#recording the data into 'ChIPseqVScontacts_plots.txt'
				if(venn_flag == 1)
				{					
					#getting parameters for ChIPseqVScontacts_plots.txt:
					struct <- "two contact bands files intersected with an peaks file - Venn Diagram"					
					bpORpeaks <- NA
					percentageOf <- NA
					if(ans4 == 1)
					{
						peaks_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						peaks_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						peaks_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "trans"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to ChIPseqVScontacts_plots.txt
					ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(vennName,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
					rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm ChIPseqVScontacts_plots.txt")
					write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}	
			}			
		}	
		else if(ans3 == 3) #chromosomes separately
		{
			#the mcnemar test
			#since we have two different 'peaks_afterCO' bed files (one for each removal of sections according to each contact bands file separately) we will perform the mcnemar test twice
			#this is not the most accurate way of doing it (since we the difference between the removals might be significant) but best for now			
			ans.mcnemar <- readline("\nwould you like to perform a mcnemar test?\ny/n\n\n")
			if(ans.mcnemar == "y")
			{
				for(u in 1:numOFchroms)
				{	
					cat("\nmcnemar test - chromosome ",u,":\n\n",sep="")
					peaks1_afterCO_mcnemar <- peaks_afterCO[[1]][peaks_afterCO[[1]][,1]==u,]
					peaks2_afterCO_mcnemar <- peaks_afterCO[[2]][peaks_afterCO[[2]][,1]==u,]
					
					#i need to test the conts and see how to extract the correct data since its a list
					#the reason why both cis[1] and cis[2] are written is just for comfort and aesthetic reasons
					conts1_mcnemar <- conts[[1]][conts[[1]][,1]==u,]
					conts2_mcnemar <- conts[[2]][conts[[2]][,1]==u,]
					write.table(peaks1_afterCO_mcnemar,"~/Analyze4C/temp/peaks1_afterCO_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
					write.table(peaks2_afterCO_mcnemar,"~/Analyze4C/temp/peaks2_afterCO_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
					write.table(conts1_mcnemar,"~/Analyze4C/temp/conts1_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
					write.table(conts2_mcnemar,"~/Analyze4C/temp/conts2_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
		
					cat("\nmcnemar test 1:\npeaks1_afterCO - \n")
					mcnemar_tester("~/Analyze4C/temp/peaks1_afterCO_mcnemar.bed","~/Analyze4C/temp/conts1_mcnemar.bed","~/Analyze4C/temp/conts2_mcnemar.bed")
					cat("\nmcnemar test 2:\npeaks2_afterCO - \n")
					mcnemar_tester("~/Analyze4C/temp/peaks2_afterCO_mcnemar.bed","~/Analyze4C/temp/conts1_mcnemar.bed","~/Analyze4C/temp/conts2_mcnemar.bed")
					
					system("rm ~/Analyze4C/temp/peaks1_afterCO_mcnemar.bed")
					system("rm ~/Analyze4C/temp/peaks2_afterCO_mcnemar.bed")
					system("rm ~/Analyze4C/temp/conts1_mcnemar.bed")
					system("rm ~/Analyze4C/temp/conts2_mcnemar.bed")
				}	
			}
		
			ans.prop <- readline(prompt=cat("\n\nwould you like to perform a prop test on the intersected bp and the intersected peaks?\ny/n\n\n"))
			if(ans.prop == "y")
			{
				#intitiating these vectors with the size of the number of chromosomes, in order to put the sums per chromosome in them
				intersectionBPsum1 <- rep(0,numOFchroms)
				intersectionBPsum2 <- rep(0,numOFchroms)
				all_contsBPsum1 <- rep(0,numOFchroms)
				all_contsBPsum2 <- rep(0,numOFchroms)
				all_peaksBPsum1 <- rep(0,numOFchroms)
				all_peaksBPsum2 <- rep(0,numOFchroms)
				intersectionPeaksSum1 <- rep(0,numOFchroms)
				intersectionPeaksSum2 <- rep(0,numOFchroms)
				all_PeaksSum1 <- rep(0,numOFchroms)
				all_PeaksSum2 <- rep(0,numOFchroms)
				contsBP_intersectedPercent1 <- rep(0,numOFchroms)
				contsBP_intersectedPercent2 <- rep(0,numOFchroms)
				peaksBP_intersectedPercent1 <- rep(0,numOFchroms)
				peaksBP_intersectedPercent2 <- rep(0,numOFchroms)			
				peaks_intersectedPercent1 <- rep(0,numOFchroms)
				peaks_intersectedPercent2 <- rep(0,numOFchroms)
				
				for(w in 1:numOFchroms)
				{
					#BPs:
					
					#getting the sum of all the bps of the intersections per chromosome
					intersectionBPsum1[w] <- sum(ChIPseqComparison1[ChIPseqComparison1[,1]==w,8])
					intersectionBPsum2[w] <- sum(ChIPseqComparison2[ChIPseqComparison2[,1]==w,8])
					
					#getting the sum of all the bps of the contact bands per chromosome
					
					#here too i need to see how im getting the data from conts @
					#maybe ill just have conts be conts1 and conts2 instead of a list
					all_contsBPsum1[w] <- sum(conts[[1]][conts[[1]][,1]==w,3]-conts[[1]][conts[[1]][,1]==w,2])
					all_contsBPsum2[w] <- sum(conts[[2]][conts[[2]][,1]==w,3]-conts[[2]][conts[[2]][,1]==w,2])
					
					#getting the sum of all the bps that are above the peaks cutoff per chromosome
					all_peaksBPsum1[w] <- sum(peaks_afterCO[[1]][peaks_afterCO[[1]][,1]==w,3]-peaks_afterCO[[1]][peaks_afterCO[[1]][,1]==w,2])
					all_peaksBPsum2[w] <- sum(peaks_afterCO[[2]][peaks_afterCO[[2]][,1]==w,3]-peaks_afterCO[[2]][peaks_afterCO[[2]][,1]==w,2])
					
					#peaks:

					#getting the sum of all the peaks of the intersections per chromosome
					ChIPseqComparison1_chr <- ChIPseqComparison1[row.names(unique(ChIPseqComparison1[,c(colnames(ChIPseqComparison1)[1],colnames(ChIPseqComparison1)[2])])),]
					intersectionPeaksSum1[w] <- sum(ChIPseqComparison1_chr[ChIPseqComparison1_chr[,1]==w,4])
					ChIPseqComparison2_chr <- ChIPseqComparison2[row.names(unique(ChIPseqComparison2[,c(colnames(ChIPseqComparison2)[1],colnames(ChIPseqComparison2)[2])])),]
					intersectionPeaksSum2[w] <- sum(ChIPseqComparison2_chr[ChIPseqComparison2_chr[,1]==w,4])
					
					#getting the sum of all the peaks above the peaks cutoff per chromosome
					all_PeaksSum1[w] <- sum(peaks_afterCO[[1]][peaks_afterCO[[1]][,1]==w,4])
					all_PeaksSum2[w] <- sum(peaks_afterCO[[2]][peaks_afterCO[[2]][,1]==w,4])

					contsBP_intersectedPercent1[w] <- (intersectionBPsum1[w]/all_contsBPsum1[w])*100
					contsBP_intersectedPercent2[w] <- (intersectionBPsum2[w]/all_contsBPsum2[w])*100				
					peaksBP_intersectedPercent1[w] <- (intersectionBPsum1[w]/all_peaksBPsum1[w])*100
					peaksBP_intersectedPercent2[w] <- (intersectionBPsum2[w]/all_peaksBPsum2[w])*100
					peaks_intersectedPercent1[w] <- (intersectionPeaksSum1[w]/all_PeaksSum1[w])*100
					peaks_intersectedPercent2[w] <- (intersectionPeaksSum2[w]/all_PeaksSum2[w])*100
				}
				
				#printing the results:
				
				cat("\nthe peaks file chosen:\n",peaks_filename,"\n")
				for(s in 1:numOFchroms)
				{
					cat("\nthe number of bps of peaks above the cutoff in chromosome",s,"on the peaks data after removed sections according to",conts_filename[1],"is:",all_peaksBPsum1[s],"\n")
					cat("\nthe number of bps of peaks above the cutoff in chromosome",s,"on the peaks data after removed sections according to",conts_filename[2],"is:",all_peaksBPsum2[s],"\n")
					cat("\nthe sum of peaks above the cutoff in chromosome",s,"on the peaks data after removed sections according to",conts_filename[1],"is:",all_PeaksSum1[s],"\n")
					cat("\nthe sum of peaks above the cutoff in chromosome",s,"on the peaks data after removed sections according to",conts_filename[2],"is:",all_PeaksSum2[s],"\n")	
					cat("\nthe number of bps of contact bands in chromosome",s,"of the contact bands file",conts_filename[1],"is:",all_contsBPsum1[s],"\n")
					cat("\nthe number of bps of contact bands in chromosome",s,"of the contact bands file",conts_filename[2],"is:",all_contsBPsum2[s],"\n")					
				}

				#perecentage of all bps of ChIPseq peaks
				cat("\ncontact bands file - ",conts_filename[1],":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis[1])
					{
						if(cis[1]==cis[2])
						{
							cat("chr",q," (cis of both contact bands files): the number of bps intersected are ",intersectionBPsum1[q]," which are ",peaksBP_intersectedPercent1[q],"% of all the bps of peaks above the cutoff in the chromosome\n",sep="")													
						}
						else
						{
							cat("chr",q," (cis of ",conts_filename[1],"): the number of bps intersected are ",intersectionBPsum1[q]," which are ",peaksBP_intersectedPercent1[q],"% of all the bps of peaks above the cutoff in the chromosome\n",sep="")													
						}
					}
					else if(q == cis[2])
					{
						cat("chr",q," (cis of ",conts_filename[2],"): the number of bps intersected are ",intersectionBPsum1[q]," which are ",peaksBP_intersectedPercent1[q],"% of all the bps of peaks above the cutoff in the chromosome\n",sep="")																		
					}
					else
					{
						cat("chr",q,": the number of bps intersected are ",intersectionBPsum1[q]," which are ",peaksBP_intersectedPercent1[q],"% of all the bps of peaks above the cutoff in the chromosome\n",sep="")				
					}
				}
				cat("\n")
				
				cat("\ncontact bands file - ",conts_filename[2],":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis[1])
					{
						if(cis[1]==cis[2])
						{
							cat("chr",q," (cis of both contact bands files): the number of bps intersected are ",intersectionBPsum2[q]," which are ",peaksBP_intersectedPercent2[q],"% of all the bps of peaks above the cutoff in the chromosome\n",sep="")													
						}
						else
						{
							cat("chr",q," (cis of ",conts_filename[1],"): the number of bps intersected are ",intersectionBPsum2[q]," which are ",peaksBP_intersectedPercent2[q],"% of all the bps of peaks above the cutoff in the chromosome\n",sep="")													
						}	
					}
					else if(q == cis[2])
					{
						cat("chr",q," (cis of ",conts_filename[2],"): the number of bps intersected are ",intersectionBPsum2[q]," which are ",peaksBP_intersectedPercent2[q],"% of all the bps of peaks above the cutoff in the chromosome\n",sep="")													
					}
					else
					{
						cat("chr",q,": the number of bps intersected are ",intersectionBPsum2[q]," which are ",peaksBP_intersectedPercent2[q],"% of all the bps of peaks above the cutoff in the chromosome\n",sep="")				
					}
				}
				cat("\n")
				
				#percentage of all bps in contact bands
				cat("\ncontact bands file - ",conts_filename[1],":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis[1])
					{
						if(cis[1]==cis[2])
						{
							cat("chr",q," (cis of both contact bands files): the number of bps intersected are ",intersectionBPsum1[q]," which are ",contsBP_intersectedPercent1[q],"% of all the bps of contact bands in the chromosome\n",sep="")													
						}
						else
						{
							cat("chr",q," (cis of ",conts_filename[1],"): the number of bps intersected are ",intersectionBPsum1[q]," which are ",contsBP_intersectedPercent1[q],"% of all the bps of contact bands in the chromosome\n",sep="")													
						}
					}
					else if(q == cis[2])
					{
						cat("chr",q," (cis of ",conts_filename[2],"): the number of bps intersected are ",intersectionBPsum1[q]," which are ",contsBP_intersectedPercent1[q],"% of all the bps of contact bands in the chromosome\n",sep="")							
					}
					else
					{
						cat("chr",q,": the number of bps intersected are ",intersectionBPsum1[q]," which are ",contsBP_intersectedPercent1[q],"% of all the bps contact bands in the chromosome\n",sep="")				
					}
				}
				cat("\n")
				
				cat("\ncontact bands file - ",conts_filename[2],":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis[1])
					{
						if(cis[1]==cis[2])
						{
							cat("chr",q," (cis of both contact bands files): the number of bps intersected are ",intersectionBPsum2[q]," which are ",contsBP_intersectedPercent2[q],"% of all the bps of contact bands in the chromosome\n",sep="")																			
						}
						else
						{
							cat("chr",q," (cis of ",conts_filename[1],"): the number of bps intersected are ",intersectionBPsum2[q]," which are ",contsBP_intersectedPercent2[q],"% of all the bps of contact bands in the chromosome\n",sep="")													
						}
					}
					else if(q == cis[2])
					{
						cat("chr",q," (cis of ",conts_filename[2],"): the number of bps intersected are ",intersectionBPsum2[q]," which are ",contsBP_intersectedPercent2[q],"% of all the bps of contact bands in the chromosome\n",sep="")																		
					}
					else
					{
						cat("chr",q,": the number of bps intersected are ",intersectionBPsum2[q]," which are ",contsBP_intersectedPercent2[q],"% of all the bps of contact bands in the chromosome\n",sep="")				
					}
				}
				cat("\n")
				
				#percentage of sum of peaks
				cat("\ncontact bands file - ",conts_filename[1],":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis[1])
					{
						if(cis[1]==cis[2])
						{
							cat("chr",q," (cis of both contact bands files): the number of peaks intersected are ",intersectionPeaksSum1[q]," which are ",peaks_intersectedPercent1[q],"% of all the peaks above the cutoff in the chromosome\n",sep="")														
						}
						else
						{
							cat("chr",q," (cis of ",conts_filename[1],"): the number of peaks intersected are ",intersectionPeaksSum1[q]," which are ",peaks_intersectedPercent1[q],"% of all the peaks above the cutoff in the chromosome\n",sep="")																				
						}
					}
					else if(q == cis[2])
					{
						cat("chr",q," (cis of ",conts_filename[2],"): the number of peaks intersected are ",intersectionPeaksSum1[q]," which are ",peaks_intersectedPercent1[q],"% of all the peaks above the cutoff in the chromosome\n",sep="")																									
					}
					else
					{
						cat("chr",q,": the number of peaks intersected are ",intersectionPeaksSum1[q]," which are ",peaks_intersectedPercent1[q],"% of all the peaks above the cutoff in the chromosome\n",sep="")				
					}	
				}
				cat("\n")	

				cat("\ncontact bands file - ",conts_filename[2],":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis[1])
					{
						if(cis[1]==cis[2])
						{
							cat("chr",q," (cis of both contact bands files): the number of peaks intersected are ",intersectionPeaksSum2[q]," which are ",peaks_intersectedPercent2[q],"% of all the peaks above the cutoff in the chromosome\n",sep="")														
						}
						else
						{
							cat("chr",q," (cis of ",conts_filename[1],"): the number of peaks intersected are ",intersectionPeaksSum2[q]," which are ",peaks_intersectedPercent2[q],"% of all the peaks above the cutoff in the chromosome\n",sep="")																				
						}
					}
					else if(q == cis[2])
					{
						cat("chr",q," (cis of ",conts_filename[2],"): the number of peaks intersected are ",intersectionPeaksSum2[q]," which are ",peaks_intersectedPercent2[q],"% of all the peaks above the cutoff in the chromosome\n",sep="")																									
					}
					else
					{
						cat("chr",q,": the number of peaks intersected are ",intersectionPeaksSum2[q]," which are ",peaks_intersectedPercent2[q],"% of all the peaks above the cutoff in the chromosome\n",sep="")				
					}	
				}
				cat("\n")

				#creating plot		
				contsName1 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",conts_filename[1],"on the plot:\n"))
				contsName2 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",conts_filename[2],"on the plot:\n"))
				contact_source <- c(replicate(numOFchroms,contsName1),replicate(numOFchroms,contsName2))
				chromosomes <- c(chr_names,chr_names)
				conts_base_pairs_percentage <- c(contsBP_intersectedPercent1,contsBP_intersectedPercent2)
				peaks_base_pairs_percentage <- c(peaksBP_intersectedPercent1,peaksBP_intersectedPercent2)
				peaks_percentage <- c(peaks_intersectedPercent1,peaks_intersectedPercent2)
				
				conts_bp_plotData <- data.frame(contact_source,chromosomes,conts_base_pairs_percentage)
				peaks_bp_plotData <- data.frame(contact_source,chromosomes,peaks_base_pairs_percentage)
				peaks_plotData <- data.frame(contact_source,chromosomes,peaks_percentage)

				#creating the percentage of bp from ChIPseq peaks plot:
				
				#cat("\nplotting the peaks bp percentage data\nplease wait...\n\n")
				#print(ggplot2::ggplot(peaks_bp_plotData, ggplot2::aes(contact_source,peaks_base_pairs_percentage,fill=chromosomes)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected bp from ChIPseq peaks between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"for each chromosome separately")) + ggplot2::labs(x="peaks source",y="intersected bp percentage(%)")) + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				
				plo_peaksBP <- ggplot2::ggplot(peaks_bp_plotData, ggplot2::aes(chromosomes,peaks_base_pairs_percentage,fill=contact_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected bp from ChIPseq peaks between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"for each chromosome separately")) + ggplot2::labs(x="peaks source",y="intersected bp percentage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the peaks\nthe 2nd creates a table of all the successes (bp of peaks intersections) and fails (all the bp of the peaks peaks that didn't intersect) and uses them\n\n")
				for(r in 1:length(unique(chromosomes)))
				{
					#prop_peaksBP <- prop.test(peaks_bp_plotData[peaks_bp_plotData$chromosomes==unique(chromosomes)[r],]$peaks_base_pairs_percentage,c(all_peaksBPsum1[r],all_peaksBPsum2[r]),correct=FALSE)
					if(intersectionBPsum1[r]>0 & intersectionBPsum2[r]>0)
					{
						cat("\nchromosome ",r,":\n",sep="")
						prop_peaksBP1 <- prop.test(c(intersectionBPsum1[r],intersectionBPsum2[r]),c(all_peaksBPsum1[r],all_peaksBPsum2[r]),correct=FALSE)
						print(prop_peaksBP1)
						prop_peaksBP2 <- prop.test(matrix(c(all_peaksBPsum1[r]-intersectionBPsum1[r],intersectionBPsum1[r],all_peaksBPsum2[r]-intersectionBPsum2[r],intersectionBPsum2[r]),ncol=2),correct=FALSE)
						print(prop_peaksBP2)
						prop_ans2 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
						if(prop_ans2 == 1)
						{
							prop_peaksBP <- prop_peaksBP1
						}
						else
						{
							prop_peaksBP <- prop_peaksBP2
						}
						
						#calculating effect sizes
						cat("\neffect sizes:\n")
						compute.es::propes((intersectionBPsum1[r]/all_peaksBPsum1[r]),(intersectionBPsum2[r]/all_peaksBPsum2[r]),all_peaksBPsum1[r],all_peaksBPsum2[r])
						cat("\n")
						
						#phi
						cat("\nphi:\n")
						print(psych::phi(matrix(c(intersectionBPsum1[r],all_peaksBPsum1[r]-intersectionBPsum1[r],intersectionBPsum2[r],all_peaksBPsum2[r]-intersectionBPsum2[r]),ncol=2,byrow=TRUE)))
						cat("\n")
						
						#creating a barplot with significance stars
						if(prop_peaksBP$p.value <= 0.05)
						{
							exps <- peaks_bp_plotData[peaks_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(peaks_bp_plotData[peaks_bp_plotData$chromosomes==unique(chromosomes)[r],]$peaks_base_pairs_percentage),]$contact_source		
							chroms <- peaks_bp_plotData[peaks_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(peaks_bp_plotData[peaks_bp_plotData$chromosomes==unique(chromosomes)[r],]$peaks_base_pairs_percentage),]$chromosomes
							perces <- peaks_bp_plotData[peaks_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(peaks_bp_plotData[peaks_bp_plotData$chromosomes==unique(chromosomes)[r],]$peaks_base_pairs_percentage),]$peaks_base_pairs_percentage+1
							label.df <- data.frame(contact_source=exps,chromosomes=chroms,peaks_base_pairs_percentage=perces)

							if(prop_peaksBP$p.value <= 0.0001)
							{
								plo_peaksBP <- plo_peaksBP + ggplot2::geom_text(data = label.df, label = "****",size=8)
							}
							else if(prop_peaksBP$p.value <= 0.001)
							{
								plo_peaksBP <- plo_peaksBP + ggplot2::geom_text(data = label.df, label = "***",size=8)
							}
							else if(prop_peaksBP$p.value <= 0.01)
							{
								plo_peaksBP <- plo_peaksBP + ggplot2::geom_text(data = label.df, label = "**",size=8)
							}
							else if(prop_peaksBP$p.value <= 0.05)
							{
								plo_peaksBP <- plo_peaksBP + ggplot2::geom_text(data = label.df, label = "*",size=8)
							}
						}
					}	
				}
				cat("\nplotting the peaks bp percentage data\nplease wait...\n\n")
				print(plo_peaksBP)
							
				ans5 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans5 == "y")
				{
					nm <- paste("peaksBP_peaksintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for ChIPseqVScontacts_plots.txt:
					struct <- "two contact bands files intersected with an peaks file"
					
					bpORpeaks <- "bp"
					
					percentageOf <- "ChIPseq peaks"

					if(ans4 == 1)
					{
						peaks_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						peaks_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						peaks_CO_appliedTo <- "each chromosome separately"
					}
					
					Intersection_chromosomes <- "each chromosome separately"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to ChIPseqVScontacts_plots.txt
					ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(nm,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
					rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm ChIPseqVScontacts_plots.txt")
					write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
				}
				
				#creating the percentage of bp intersected from contact bands plot:
				
				#cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")
				#print(ggplot2::ggplot(conts_bp_plotData, ggplot2::aes(contact_source,conts_base_pairs_percentage,fill=chromosomes)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected bp from contact bands between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"for each chromosome separately")) + ggplot2::labs(x="peaks source",y="intersected bp percentage(%)")) + ggplot2::theme(plot.title = ggplot2::element_text(size=13))

				plo_contsBP <- ggplot2::ggplot(conts_bp_plotData, ggplot2::aes(chromosomes,conts_base_pairs_percentage,fill=contact_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected bp from contact bands between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"for each chromosome separately")) + ggplot2::labs(x="peaks source",y="intersected bp percentage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				cat("\nperforming prop tests\nthe 1st is the comparison of proportions of intersections out of all bps in the contact bands\nthe 2nd creates a table of all the successes (bp intersections) and fails (all the bp of the contact bands that didn't intersect) and uses them\n\n")
				for(r in 1:length(unique(chromosomes)))
				{
					#prop_contsBP <- prop.test(conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],]$conts_base_pairs_percentage,c(all_contsBPsum1[r],all_contsBPsum2[r]),correct=FALSE)
					if(intersectionBPsum1[r]>0 & intersectionBPsum2[r]>0)
					{
						cat("\nchromosome ",r,":\n",sep="")
						prop_contsBP1 <- prop.test(c(intersectionBPsum1[r],intersectionBPsum2[r]),c(all_contsBPsum1[r],all_contsBPsum2[r]),correct=FALSE)
						print(prop_contsBP1)
						prop_contsBP2 <- prop.test(matrix(c(all_contsBPsum1[r]-intersectionBPsum1[r],intersectionBPsum1[r],all_contsBPsum2[r]-intersectionBPsum2[r],intersectionBPsum2[r]),ncol=2),correct=FALSE)
						print(prop_contsBP2)
						prop_ans <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
						if(prop_ans == 1)
						{
							prop_contsBP <- prop_contsBP1
						}
						else
						{
							prop_contsBP <- prop_contsBP2
						}
						
						#calculating effect sizes
						cat("\neffect sizes:\n")
						compute.es::propes((intersectionBPsum1[r]/all_contsBPsum1[r]),(intersectionBPsum2[r]/all_contsBPsum2[r]),all_contsBPsum1[r],all_contsBPsum2[r])
						cat("\n")
						
						#phi
						cat("\nphi:\n")
						print(psych::phi(matrix(c(intersectionBPsum1[r],all_contsBPsum1[r]-intersectionBPsum1[r],intersectionBPsum2[r],all_contsBPsum2[r]-intersectionBPsum2[r]),ncol=2,byrow=TRUE)))
						cat("\n")
						
						#creating a barplot with significance stars
						if(prop_contsBP$p.value <= 0.05)
						{
							exps <- conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],]$conts_base_pairs_percentage),]$contact_source		
							chroms <- conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],]$conts_base_pairs_percentage),]$chromosomes
							perces <- conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],]$conts_base_pairs_percentage),]$conts_base_pairs_percentage+1
							label.df <- data.frame(contact_source=exps,chromosomes=chroms,conts_base_pairs_percentage=perces)

							if(prop_contsBP$p.value <= 0.0001)
							{
								plo_contsBP <- plo_contsBP + ggplot2::geom_text(data = label.df, label = "****",size=8)
							}
							else if(prop_contsBP$p.value <= 0.001)
							{
								plo_contsBP <- plo_contsBP + ggplot2::geom_text(data = label.df, label = "***",size=8)
							}
							else if(prop_contsBP$p.value <= 0.01)
							{
								plo_contsBP <- plo_contsBP + ggplot2::geom_text(data = label.df, label = "**",size=8)
							}
							else if(prop_contsBP$p.value <= 0.05)
							{
								plo_contsBP <- plo_contsBP + ggplot2::geom_text(data = label.df, label = "*",size=8)
							}
						}
					}	
				}
				cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")
				print(plo_contsBP)
							
				ans8 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans8 == "y")
				{
					nm <- paste("contactsBP_peaksintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for ChIPseqVScontacts_plots.txt:
					struct <- "two contact bands files intersected with an peaks file"
					
					bpORpeaks <- "bp"
					
					percentageOf <- "Contact bands"

					if(ans4 == 1)
					{
						peaks_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						peaks_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						peaks_CO_appliedTo <- "each chromosome separately"
					}
					
					Intersection_chromosomes <- "each chromosome separately"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to ChIPseqVScontacts_plots.txt
					ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(nm,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
					rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm ChIPseqVScontacts_plots.txt")
					write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
				}
				
				#creating the percentage of peaks value (p-score or number of tags) intersected of sum of the same peaks value plot:
				
				#cat("\nplotting the peaks percentage data\nplease wait...\n\n")
				#print(ggplot2::ggplot(peaks_plotData, ggplot2::aes(contact_source,peaks_percentage,fill=chromosomes)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected peaks between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"for each chromosome separately")) + ggplot2::labs(x="peaks source",y="intersected peaks percentage(%)")) + ggplot2::theme(plot.title = ggplot2::element_text(size=13))

				plopeaks <- ggplot2::ggplot(peaks_plotData, ggplot2::aes(chromosomes,peaks_percentage,fill=contact_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected peaks between peaks of",peaks_tissue,"and contacts of",sp1[[1]][1],sp1[[1]][2],sp1[[1]][3],"and",sp1[[2]][1],sp1[[2]][2],sp1[[2]][3],"for each chromosome separately")) + ggplot2::labs(x="peaks source",y="intersected peaks percentage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				cat("\nthe 1st is the comparison of proportions of intersections out of all peaks values in the peaks\nthe 2nd creates a table of all the successes (peaks intersections) and fails (all the peaks of the peaks that didn't intersect) and uses them\n\n")					
				for(r in 1:length(unique(chromosomes)))
				{
					if(intersectionPeaksSum1[r]>0 & intersectionPeaksSum2[r]>0)
					{
						cat("\nchromosome ",r,":\n",sep="")
						#proppeaks <- prop.test(peaks_plotData[peaks_plotData$chromosomes==unique(chromosomes)[r],]$peaks_percentage,c(all_PeaksSum1[r],all_PeaksSum2[r]),correct=FALSE)
						proppeaks1 <- prop.test(c(intersectionPeaksSum1[r],intersectionPeaksSum2[r]),c(all_PeaksSum1[r],all_PeaksSum2[r]),correct=FALSE)
						print(proppeaks1)
						proppeaks2 <- prop.test(matrix(c(all_PeaksSum1[r]-intersectionPeaksSum1[r],intersectionPeaksSum1[r],all_PeaksSum2[r]-intersectionPeaksSum2[r],intersectionPeaksSum2[r]),ncol=2),correct=FALSE)
						print(proppeaks2)					
						prop_ans3 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
						if(prop_ans3 == 1)
						{
							proppeaks <- proppeaks1
						}
						else
						{
							proppeaks <- proppeaks2
						}
						
						#calculating effect sizes
						cat("\neffect sizes:\n")
						compute.es::propes((intersectionPeaksSum1[r]/all_PeaksSum1[r]),(intersectionPeaksSum2[r]/all_PeaksSum2[r]),all_PeaksSum1[r],all_PeaksSum2[r])
						cat("\n")
						
						#phi
						cat("\nphi:\n")
						print(psych::phi(matrix(c(intersectionPeaksSum1[r],all_PeaksSum1[r]-intersectionPeaksSum1[r],intersectionPeaksSum2[r],all_PeaksSum2[r]-intersectionPeaksSum2[r]),ncol=2,byrow=TRUE)))
						cat("\n")
						
						#creating a barplot with significance stars
						if(proppeaks$p.value <= 0.05)
						{
							exps <- peaks_plotData[peaks_plotData$chromosomes==unique(chromosomes)[r],][which.max(peaks_plotData[peaks_plotData$chromosomes==unique(chromosomes)[r],]$peaks_percentage),]$contact_source		
							chroms <- peaks_plotData[peaks_plotData$chromosomes==unique(chromosomes)[r],][which.max(peaks_plotData[peaks_plotData$chromosomes==unique(chromosomes)[r],]$peaks_percentage),]$chromosomes
							perces <- peaks_plotData[peaks_plotData$chromosomes==unique(chromosomes)[r],][which.max(peaks_plotData[peaks_plotData$chromosomes==unique(chromosomes)[r],]$peaks_percentage),]$peaks_percentage+1
							label.df <- data.frame(contact_source=exps,chromosomes=chroms,peaks_percentage=perces)

							if(proppeaks$p.value <= 0.0001)
							{
								plopeaks <- plopeaks + ggplot2::geom_text(data = label.df, label = "****",size=8)
							}
							else if(proppeaks$p.value <= 0.001)
							{
								plopeaks <- plopeaks + ggplot2::geom_text(data = label.df, label = "***",size=8)
							}
							else if(proppeaks$p.value <= 0.01)
							{
								plopeaks <- plopeaks + ggplot2::geom_text(data = label.df, label = "**",size=8)
							}
							else if(proppeaks$p.value <= 0.05)
							{
								plopeaks <- plopeaks + ggplot2::geom_text(data = label.df, label = "*",size=8)
							}
						}
					}	
				}
				cat("\nplotting the peaks percentage data\nplease wait...\n\n")
				print(plopeaks)
				
				ans6 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans6 == "y")
				{
					nm <- paste("peaks_peaksintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for ChIPseqVScontacts_plots.txt:
					struct <- "two contact bands files intersected with an peaks file"
					
					bpORpeaks <- "peaks value"
					
					percentageOf <- "peaks value"
					
					if(ans4 == 1)
					{
						peaks_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						peaks_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						peaks_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "each chromosome separately"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to ChIPseqVScontacts_plots.txt
					ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(nm,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
					rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm ChIPseqVScontacts_plots.txt")
					write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}
			}

			#venn diagrams
			if(readline(prompt=cat("\nwould you like to create venn diagrams?\ny/n\n\n")) == "y")
			{
				#creating the venn data frame for 0 percent
				venn_DF <- data.frame(rep(1,nrow(peaks_afterCO[[1]])),matrix(0,nrow(peaks_afterCO[[1]]),2))##!

				cat("\nvenn diagram parameters:\n\n")
				venn_ans1 <- readline(prompt=cat("\nwould you like to consider overlaps to be only those that intersect more than 1 bp?\ny/n\n\n"))
				if(venn_ans1 == "y")
				{
					ovlp_cov <- as.numeric(readline(prompt=cat("\nenter the minimum overlap coverage that will be considered an intersection? (between 0 and 1)\n\n")))
					
					venn_ans2 <- readline(prompt=cat("\nyou have chosen to accept an intersection only if there is at least",ovlp_cov,"overlap.\nwould you consider accepting an intersection if there is more than a certain amount of reads that intersect (e.g. if you there are more than an x amount of reads intersected, even though each one doesn't overlap more than",ovlp_cov,", altogether they are sufficient enough to be considered as an overlap)\ny/n\n\n"))
					if(venn_ans2 == "y")
					{
						ovlp_num <- as.integer(readline(prompt=cat("\nhow many reads (minimum) intersected will be considered as an overlap?\n\n")))
					}
				}
				
				venn_flag <- 0 #if venn_flag is 1 then we save a diagram and it should be recorded in 'ChIPseqVScontacts_plots'				
				#filling in the data for the venn diagrams
				if(venn_ans1 == "y")
				{
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[1]," -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[2]," -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A2.bed",sep=""))
					venn_A1 <- read.table(paste("~/Analyze4C/temp/venn_A1.bed",sep=""))
					venn_A2 <- read.table(paste("~/Analyze4C/temp/venn_A2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn_A1.bed")
					system("rm ~/Analyze4C/temp/venn_A2.bed")
					if(venn_ans2 == "y")
					{
						system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[1]," -c > ~/Analyze4C/temp/venn_B1.bed",sep=""))
						system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[2]," -c > ~/Analyze4C/temp/venn_B2.bed",sep=""))
						venn_B1 <- read.table(paste("~/Analyze4C/temp/venn_B1.bed",sep=""))
						venn_B2 <- read.table(paste("~/Analyze4C/temp/venn_B2.bed",sep=""))
						#removing files
						system("rm ~/Analyze4C/temp/venn_B1.bed")
						system("rm ~/Analyze4C/temp/venn_B2.bed")
						#if the number in the 1st is 0 but the number in the 2nd is more than ovlp_num then we will consider it to be 1
						venn_A1[venn_B1[,5]>=ovlp_num,5] <- 1
						venn_A2[venn_B2[,5]>=ovlp_num,5] <- 1
					}
					#add the correct number to each column and create venn diagram using the data frame venn_DF
					for(w in 1:numOFchroms)
					{
						cat("\nchromosome ",w,":\n\n",sep="")
						venn_DF[venn_A1[,5]>=1 & venn_A1[,1]==w,2] <- 1
						venn_DF[venn_A2[,5]>=1 & venn_A2[,1]==w,3] <- 1
						vennName <- paste("venn_peaksVS2psConts_chr",w,"_",DandT2,".jpg",sep="")
						venn_flag <- venn_creator(venn_DF,2,1,vennName)						
					}					
				}
				else
				{
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[1]," -c > ~/Analyze4C/temp/venn1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename[2]," -c > ~/Analyze4C/temp/venn2.bed",sep=""))
					venn1 <- read.table(paste("~/Analyze4C/temp/venn1.bed",sep=""))
					venn2 <- read.table(paste("~/Analyze4C/temp/venn2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn1.bed")
					system("rm ~/Analyze4C/temp/venn2.bed")
					#add the correct number to each column and create venn diagram using the data frame venn_DF
					for(w in 1:numOFchroms)
					{
						cat("\nchromosome ",w,":\n\n",sep="")
						venn_DF[venn1[,5]>=1 & venn1[,1]==w,2] <- 1
						venn_DF[venn2[,5]>=1 & venn2[,1]==w,3] <- 1
						vennName <- paste("venn_peaksVS2psConts_chr",w,"_",DandT2,".jpg",sep="")
						venn_flag <- venn_creator(venn_DF,2,1,vennName)	
					}					
				}
				
				#recording the data into 'ChIPseqVScontacts_plots.txt'
				if(venn_flag == 1)
				{					
					#getting parameters for ChIPseqVScontacts_plots.txt:
					struct <- "two contact bands files intersected with an peaks file - Venn Diagram"					
					bpORpeaks <- NA
					percentageOf <- NA
					if(ans4 == 1)
					{
						peaks_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						peaks_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						peaks_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "each chromosome separately"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the ChIPseqVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to ChIPseqVScontacts_plots.txt
					ChIPseqVScontacts_plots[nrow(ChIPseqVScontacts_plots)+1,] <- c(vennName,struct,peaks_filename,conts_filename[1],conts_filename[2],bpORpeaks,percentageOf,DandT1,CO_col,CO_type,peaks_cutoff,peaks_CO_appliedTo,value_col,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					ChIPseqVScontacts_plots <- ChIPseqVScontacts_plots[order(ChIPseqVScontacts_plots$Plotfile_name),]
					rownames(ChIPseqVScontacts_plots) <- seq(length=nrow(ChIPseqVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm ChIPseqVScontacts_plots.txt")
					write.table(ChIPseqVScontacts_plots,"ChIPseqVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}	
			}						
		}		
		else if(ans3 == 4) #exit
		{
			#removing the files in temp
			system("rm ~/Analyze4C/temp/peaks1_afterCO.bed")
			system("rm ~/Analyze4C/temp/peaks2_afterCO.bed")
			system("rm ~/Analyze4C/temp/ChIPseqComparison1.bed")
			system("rm ~/Analyze4C/temp/ChIPseqComparison2.bed")
			#system("rm ~/Analyze4C/temp/*")
			break
		}				
	}
}
