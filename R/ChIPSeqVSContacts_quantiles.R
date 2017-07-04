#' @export

#anywhere we ask 'if(z == 1)' its in order to calculate the new p-score data but avoid multiple calculations of peaks

#there are a couple of ways to compare the data:
#1) divide the data into quantiles and compare (intesect) the different quantiles p-scores (turned into contacts) with the same quantiles of peaks - 
# this method will show if the highest expression of a certain tissue correlates with the most significant contact of a tissue
#2) do the same as above just without quantiles, just count to see where the contacts and peaks intersect
#3) create windows, sum the values of contacts and peaks and correlate between them
#4) get the reads under the contacts and compare the number of reads with the peaks under the contacts
#5) look at the contact areas and at the peaks under them, compare between the root and leaf to see if where there is a root contact the expression is higher for root and vice versa

#i need to test this way more. it seems that the results for multiple files are the same. this is wrong probably

#if the user wants to test the whole genome without dividing it into quantiles, they just need to enter the quantile of 0%

#this function also returns a venn diagram of the contact bands that intersect with the peaks data

ChIPSeqVSContacts_quantiles <- function(Experiments_4C,ChIPseqVScontacts_sumOFintersections_plots,rearranged_rawData)
{
	#starting by clearing out all files from 'temp'
	if(system("ls ~/Analyze4C/temp/") != 0)
	{
		system("rm ~/Analyze4C/temp/*")
	}
	
	#if the user chooses the quantile method, then we need to create contact bands according to the p-score quantiles
	#otherwise just take a contact band file

	#i should probably have an option to do all the chromosomes separately, since the p scores are calculated separately the quantiles are decieving
	#i could get quantiles for each chromosome and then add the contacts together into one file or data frame
	cat("\n**************************************************************************\n\n")
	cat("\n                    ChIPseq Vs. Contacts - by quantile\n\n")
	cat("\n**************************************************************************\n\n")
	cat("\nyou can choose to test multiple experiments.\nmeaning the same ChIPseq data (tags or p-scores) will be compared to the number of experiments you choose.\n")
	numOFexperiments <- as.integer(readline(prompt=cat("\nHow many experiments would you like to use?\n\n")))
	
	#choosing a genome size file
	ls_genomes <- system("ls ~/Analyze4C/genomes/",intern=TRUE)
	genome_file.name <- readline(prompt=cat("\nchoose the file with genome sizes that you would like to use:\n",ls_genomes,"",sep="\n"))
	genome_sizes <- read.table(paste("~/Analyze4C/genomes/",genome_file.name,sep=""))
	
	numOFchroms <- nrow(genome_sizes)
	
	choice2 <- as.integer(readline(prompt=cat("\nwhat chromosomes should the intersection be applied to? (enter the option number)\n1) all chromosomes\n2) trans\n3) cis\n4) specific chromosome\n\n")))
	summer_all <- list() #this will contain all the entries (from each experiment data) of summing the values
	summer_peaks_self_all <- list()
	summer_ps_self_all <- list()
	all_file_names <- rep(0,numOFexperiments) #this will contain the names of the chosen p-score files
	translocated_flag <- -1 #this has 2 purposes: to let the function know if the current p-score file chosen is a translocation switched file (==1), and if there has previously (in the current run) been a usage of a p-score file that wasn't a translocation switched
	first <- -1 #as long as first equals -1 we know there was no translocation removed
	all_rem <- list() #this will contain all the details of the sections removed
	rem_counter <- 1 #this will counted how many sections were removed
	venn_col_num <- 1 #this is the index of the column for each file,column 1 is occupied by the ChIPseq peaks
	venn_DF <- list(0) #creating the empty list with the data frames for the venn diagrams

	for(z in 1:numOFexperiments)
	{
		cat("\n******************************************************************************************\n")
		cat("\nexperiment number ",z,":\n",sep="")

		#initializing the venn index to 1
		venn_ind <- 1
		venn_col_num <- venn_col_num + 1

		#choosing p-score and peaks files, and getting the specific data from those files
		{
			#choose a p-score file:
			#getting the name of the p score file the user wants to use
			repeat
			{
				choice3 <- as.integer(readline(prompt=cat("\nwhat folder of p-score files would you like to choose from? (enter the number)\n1) no bp windows\n2) with bp windows\n\n")))
				if(choice3 == 1) #no bp windows
				{
					ps_folder <- "no_bp"
				}
				else if(choice3 == 2) #with bp windows
				{
					ps_folder <- "with_bp"
				}	

				ps_files <- system(paste("ls ~/Analyze4C/pScores/",ps_folder,sep=""),intern=TRUE)
				ps.file <- readline(prompt=cat("These are the p-score files available\nwhich would you like to use?\n",ps_files,"",sep="\n"))
				ind_files2 <- pmatch(ps.file,ps_files,nomatch=0)
				if(ind_files2 == 0)
				{
					cat("no such file exists.\nplease try again.\n\n")
				}
				else
				{
					#here we check to see if the p-score file chosen is of a translocation switched type, and if so that they aren't removed
					if(identical(grep("translocatedSwitched",ps.file),integer(0)) == "FALSE")
					{
						if(translocated_flag == 0)
						{
							cat("\nerror: you cannot choose this file since you have already chosen a file that was not with a translocation switched\n")
							cat("this will 'mess up' the calculation when you will need to remove the translocation sections\n")
							ans3 <- as.integer(readline(prompt=cat("would you like to:\n1) choose a different file\n2) start the analysis from the beginning - this will delete the calculation of the files already chosen. if you want the same files we advise you choose the translocation switched files first, then the others\n\n")))
							if(ans3 == 2)
							{
								summer_all <- list() #this will contain all the entries (from each experiment data) of summing the values
								summer_peaks_self_all <- list()
								summer_ps_self_all <- list()
								all_file_names <- rep(0,numOFexperiments) #this will contain the names of the chosen p-score files
								translocated_flag <- -1 #this has 2 purposes: to let the function know if the current p-score file chosen is a translocation switched file (==1), and if there has previously (in the current run) been a usage of a p-score file that wasn't a translocation switched
								first <- -1 #as long as first equals -1 we know there was no translocation removed
								all_rem <- list() #this will contain all the details of the sections removed
								rem_counter <- 1 #this will counted how many sections were removed
								z <- 1 #initializing z again, since we are starting the loop all over again
								cat("\nexperiment number ",z,":\n\n",sep="")
								
								#i might need to add here an option to remove some files -  im not sure about this, it could be that i dont need it
							}
						}
						else if(identical(grep("removed",ps.file),integer(0)) == "FALSE")
						{
							cat("\nerror: the p-score chosen had the translocation areas switched and removed\nyou must choose the p-score file which still contain the translocated areas\n\n")								
						}
						else
						{
							translocated_flag <- 1
							all_file_names[z] <- ps.file
							break
						}
					}
					else
					{
						translocated_flag <- 0
						all_file_names[z] <- ps.file
						break
					}
				}
			}
			
			if(z==1)
			{
				exps <- ps.file
			}
			else
			{
				exps <- paste(exps,ps.file)
			}
			
			#importing the p-score data from the file 	
			ps <- read.table(paste("~/Analyze4C/pScores/",ps_folder,"/",ps.file,sep=""))

			#finding the specific raw data files information in Experiments_4C
			sp3 <- unlist(strsplit(ps.file,"[.]"))
			sp4 <- strsplit(sp3[[1]][1],"_")
			
			#getting the size of the window per RE sites
			sp5 <- sp4[[1]][as.integer(grep("RE$",sp4[[1]]))]
			RE_wind <- strsplit(sp5[[1]][1],"RE")
			RE_wind <- as.numeric(RE_wind)
			
			if(choice3 == 2)
			{
				#getting the size of the window per bp sites
				sp6 <- sp4[[1]][as.integer(grep("bp$",sp4[[1]]))]
				bp_wind <- strsplit(sp6[[1]][1],"bp")
				bp_wind <- as.numeric(bp_wind)
			}

			#finding the specific raw data files information in Experiments_4C			
			out <- findIn.Experiments_4C(sp4[[1]][1],sp4[[1]][2],sp4[[1]][3],Experiments_4C)
			
			#getting the cis chromosome number
			cis <- as.numeric(out[2])
			
			if(z == 1)
			{
				#getting the peaks file	
				#ls_files_peaks <- system("ls ~/Analyze4C/ChIPseq | grep '.bed$'",intern=TRUE)
				ls_files_peaks <- system("ls ~/Analyze4C/ChIPseq/peaks | grep 'tagsANDpscores.bed$'",intern=TRUE)
				repeat
				{	
					file.name_peaks <- readline(prompt=cat("\nThese are the ChIPseq data (tags or p-scores) files available\nwhich would you like to use?\n",ls_files_peaks,"",sep="\n"))
					ind_files_peaks <- pmatch(file.name_peaks,ls_files_peaks,nomatch=0)
					if(ind_files_peaks == 0)
					{
						cat("no such file exists.\nplease try again.\n\n")
					}
					else
					{
						#importing the data from the file
						#peaks <- read.table(paste("~/Analyze4C/ChIPseq/",file.name_peaks,sep=""))
						peaks <- read.table(paste("~/Analyze4C/ChIPseq/peaks/",file.name_peaks,sep=""))
						break
					}

					#here i might need to edit and sort the data in peaks
					#peaks <- peaks[peaks[,1]!="Pt" & peaks[,1]!="Mt",]
					#peaks <- peaks[order(peaks[,1],peaks[,2]),]
					#names(peaks) <- c("chr","first","last","peaks")
				}
			}

			#getting the p-score and peaks data by chromosome according to the request of the user, and also the chromosomes chosen
			if(choice2 == 1) #all chromosomes
			{
				ps_dat <- ps
				peaks_dat <- peaks
				chroms <- c(1:numOFchroms)
				int_chroms <- "all" #this will be entered later in 'ChIPseqVScontacts_sumOFintersections_plots.txt'
			}
			else if(choice2 == 2) #trans
			{
				ps_dat <- ps[ps[,1] != cis,]
				peaks_dat <- peaks[peaks[,1] != cis,]
				chroms <- Filter(function(x) x!=cis,1:numOFchroms) #this here excludes cis chromosome number in the sequence, another way would be - (1:numOFchroms)[-cis]
				int_chroms <- "trans" #this will be entered later in 'ChIPseqVScontacts_sumOFintersections_plots.txt'
			}
			else if(choice2 == 3) #cis
			{
				ps_dat <- ps[ps[,1] == cis,]
				peaks_dat <- peaks[peaks[,1] == cis,]			
				chroms <- cis
				int_chroms <- "cis" #this will be entered later in 'ChIPseqVScontacts_sumOFintersections_plots.txt'
			}
			else if(choice2 == 4) #by chromosome
			{
				chr_chosen <- as.integer(readline(prompt=cat("\nchoose a chromosome you would like to use:\n\n")))
				ps_dat <- ps[ps[,1] == chr_chosen,]
				peaks_dat <- peaks[peaks[,1] == chr_chosen,]			
				chroms <- chr_chosen
				int_chroms <- paste("chr",chr_chosen,sep="") #this will be entered later in 'ChIPseqVScontacts_sumOFintersections_plots.txt'
			}
			write.table(ps_dat,"~/Analyze4C/temp/ps_dat_temp.sgr",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t") #creating a temp file of the p-score data chosen
		}	
		
		#choosing the value from which the qunatiles are created and the value which is going to be used for each peak in the quantile
		cat("\ndividing the ChIPseq data into quantiles:\n")
		ans4 <- as.integer(readline(prompt=cat("choose an option:\n1) get quantiles from p-scores, peak values from p-scores\n2) get quantiles from p-scores, peak values from number of tags\n3) get quantiles from number of tags, peak values from p-scores\n4) get quantiles from number of tags, peak values from number of tags\n\n")))
		if(ans4==1)
		{
			quant_col <- 5
			value_col <- 5
			quant_source <- "p-scores"
			value_source <- "p-scores"
		}
		else if(ans4==2)
		{
			quant_col <- 5
			value_col <- 4
			quant_source <- "p-scores"
			value_source <- "tags"			
		}
		else if(ans4==3)
		{
			quant_col <- 4
			value_col <- 5
			quant_source <- "tags"
			value_source <- "p-scores"			
		}
		else if(ans4==4)
		{
			quant_col <- 4
			value_col <- 4
			quant_source <- "tags"
			value_source <- "tags"			
		}
		
		if(z == 1)
		{
			#getting the quantile percentages
			{
				quant_flag <- 0
				quant_counter <- 0
				quant_percentages <- c() #this will contain all the quantiles
				cat("\nchoosing the quantiles:\n\n")
				while(quant_flag == 0)
				{
					quant_counter <- quant_counter + 1
					quant_percentages_temp <- -1
					while(quant_percentages_temp < 0 | quant_percentages_temp > 1)
					{
						quant_percentages_temp <- as.numeric(readline(prompt=cat("\nplease enter quantile number ",quant_counter,": (must be between 0 and 1)\n\n",sep="")))
					}	
					quant_percentages <- c(quant_percentages,quant_percentages_temp)
					cat("\nthe current quantiles are: ",quant_percentages,"\n")
					ans1 <- readline(prompt=cat("would you like to add another quantile?\ny/n\n\n"))
					if(ans1 == "n")
					{
						quant_flag <- 1
					}
				}
				quant_percentages <- sort(quant_percentages) #this will make sure the quantiles are in order, this makes things easier
				numOfquants <- length(quant_percentages) #contains the number of quantiles
			}
		}
		
		#calculating the quantiles for the p-score and peaks data
		{
			if(z == 1)
			{
				#getting parameters for the contact band making
				RE_gap <- as.numeric(readline(prompt=cat("\nplease choose up to how many negative RE sites in between positive RE sites it is not considered a break in the contact band:\n\n")))
				bp_gap <- as.numeric(readline(prompt=cat("\nplease choose the max number of bp that are allowed to be in between RE sites and not be considered a gap:\n\n")))
			}	
			
			quants_ps <- c() #this will contain the quantiles of the p score data
			maxes_ps <- c()
			mins_ps <- c()
			
			if(z == 1)
			{
				quants_peaks <- c() #this will contain the quantiles of the peaks data
				maxes_peaks <- c()
				mins_peaks <- c()
			}
			
			#getting the quantiles from the p-score and peaks data, by chromosome
			for(k in 1:length(chroms))
			{
				#p-score data
				quants_ps_bychrom <- quantile(ps_dat[ps_dat[,1]==chroms[k],3],quant_percentages)
				quants_ps <- rbind(quants_ps,quants_ps_bychrom)
				#getting the minimum and maximum values of p-score for each chromosome
				max_ps_temp <- max(ps_dat[ps_dat[,1]==chroms[k],3])
				maxes_ps <- c(maxes_ps,max_ps_temp)
				min_ps_temp <- min(ps_dat[ps_dat[,1]==chroms[k],3])
				mins_ps <- c(mins_ps,min_ps_temp)
				
				if(z == 1)
				{
					#peaks data
					quants_peaks_bychrom <- quantile(peaks_dat[peaks_dat[,1]==chroms[k],quant_col],quant_percentages)
					quants_peaks <- rbind(quants_peaks,quants_peaks_bychrom)
					#getting the minimum and maximum values of peaks for each chromosome
					max_peaks_temp <- max(peaks_dat[peaks_dat[,1]==chroms[k],quant_col])
					maxes_peaks <- c(maxes_peaks,max_peaks_temp)
					min_peaks_temp <- min(peaks_dat[peaks_dat[,1]==chroms[k],quant_col])
					mins_peaks <- c(mins_peaks,min_peaks_temp)
				}	
			}
			row.names(quants_ps) <- chroms #changing the row names so that they fit each chromosome
			if(z == 1)
			{
				row.names(quants_peaks) <- chroms #changing the row names so that they fit each chromosome
			}	
		}

		#creating the contact bands for p-score data
		{
			#removing the translocated areas if a file of this type was chosen:
			if(translocated_flag == 1)
			{
				cat("\nremoving translocated areas:\nthe purpose of this is in order to be able and intersect the contact bands with peaks bands\n")
				cat("the translocated areas indices mess up the order and will cause an error\n")
				cat("make sure that no translocated area is left otherwise it will crash the program\n")
				cat("it is suggested that you view the p-score file - 'ps_dat_temp.sgr', in order to see what sections to remove\n\n")
				
				len2 <- nrow(ps_dat)
				ps_removed <- cbind(ps_dat,rep(1,len2))
				rem_temp <- 1 #indicates if the user wants to remove another section, when equals 0 they do
				while(rem_temp == 1)
				{
					ind_equal <- 0 #tells us if the indices of the first and second are on the same chromosome or not
					while(ind_equal == 0)
					{
						cat("\nenter the indices of the section that you would like be removed (make sure the first and second are on the same chromosome)\n\n")
						first <- as.integer(readline(prompt="\nenter the first index of the section:\n"))
						second <- as.integer(readline(prompt="\nenter the second index of the section:\n"))
						if(ps_removed[first,1] == ps_removed[second,1])
						{
							ind_equal <- 1
						}
						else
						{
							cat("\nerror: the first and second indices should be on the same chromosome\nplease enter them again\n\n")
						}
					}	
					all_rem[[rem_counter]] <- rbind(ps_removed[first,],ps_removed[second,]) #entering the removed section indices into all_rem
					ps_removed[first:second,ncol(ps_removed)] <- 0 #the chosen removed section will get a 0 and will later be excluded
					ans2 <- readline(prompt="\nwould you like to remove another section?\ny/n\n")
					if(ans2 != "y")
					{
						rem_temp <- 0
						rem_counter <- rem_counter + 1
					}
				}
				ps_dat <- ps_removed[ps_removed[,ncol(ps_removed)]==1,1:3]
				rownames(ps_dat) <- seq(length=nrow(ps_dat))		
			}
			else if(translocated_flag == 0 & first != -1) #if there was a translocation removal already on a different file, we must remove the same sections here
			{
				for(b in 1:length(all_rem))
				{
					len2 <- nrow(ps_dat)
					ps_removed <- cbind(ps_dat,rep(1,len2))
					ps_removed[ps_removed[,1]==all_rem[[b]][1,1] & (ps_removed[,2]>=all_rem[[b]][1,2] | ps_removed[,1]<=all_rem[[b]][2,2]),ncol(ps_removed)] <- 0 #this section will get a 0 and will later be excluded							
				}
				ps_dat <- ps_removed[ps_removed[,ncol(ps_removed)]==1,1:3]
				rownames(ps_dat) <- seq(length=nrow(ps_dat))
			}
			
			#getting the p scores for each quantile and creating bed files of this data:
			
			#getting the quantile p-scores and peaks of the range of 0% till the first quantile percent the user entered. This will not be performed if the first quantile is 0% anyhow. and then creating contact bands of this range.
			if(quant_percentages[1] != 0)
			{
				#getting the p-score and peaks data of the first range (for all chromosomes)
				ps_range <- c() #this will contain the p-score range data
				peaks_range <- c() #this will contain the peaks range data	
				
				for(h in 1:length(chroms))
				{
					#p-score data
					ps_chrom <- ps_dat[ps_dat[,1]==chroms[h],] #getting the p score of a specific chromosome
					ps_temp <- ps_chrom[ps_chrom[,3]>=mins_ps[h] & ps_chrom[,3]<quants_ps[h,1],] #here we get the p-scores in a range of a quant
					ps_range <- rbind(ps_range,ps_temp)

					#peaks data
					peaks_chrom <- peaks_dat[peaks_dat[,1]==chroms[h],] #getting the peaks of a specific chromosome
					peaks_temp <- peaks_chrom[peaks_chrom[,quant_col]>=mins_peaks[h] & peaks_chrom[,quant_col]<quants_peaks[h,1],c(1:3,value_col)] #here we get the peaks in a range of a quant
					peaks_range <- rbind(peaks_range,peaks_temp)
				}

				#creating contact bands for the p-score data of the first range
				cat("\ncalculating the contact bands for the quantile range of 0% and ",quant_percentages[1]*100,"%\nplease wait...\n\n",sep="")	
				ps_conts <- contactBands_byChrom(ps_range,RE_gap,bp_gap,0,genome_sizes,rearranged_rawData)
				ps_conts <- ps_conts[[1]]
				
				#creating the bed files for the p-score and peaks data in the first range
				write.table(ps_conts,paste("~/Analyze4C/temp/ps_conts_0per.bed",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
				write.table(peaks_range,paste("~/Analyze4C/temp/peaks_0per.bed",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")					

				#creating the venn data frame for 0 percent
				if(z == 1)
				{
					#venn_DF[[venn_ind]] <- list(0)
					#venn_DF[[venn_ind]][[1]] <- as.character(paste(peaks_range[,1],peaks_range[,2],peaks_range[,3]))
					venn_DF[[venn_ind]] <- data.frame(rep(1,nrow(peaks_range)),matrix(0,nrow(peaks_range),numOFexperiments))
					venn_ind <- venn_ind + 1
				}
			}
			
			#getting the quantile p-scores of the ranges of quantile percent the user entered (except for the last). and then creating contact bands of this range. this will only be performed if there are more than two quantile percentages entered.
			if(numOfquants > 1)
			{
				for(j in 1:(numOfquants-1))
				{
					#getting the p-score and peaks data of the first range (for all chromosomes)
					ps_range <- c() #this will contain the p-score range data
					peaks_range <- c() #this will contain the peaks range data
					
					for(h in 1:length(chroms))
					{
						#p-score data			
						ps_chrom <- ps_dat[ps_dat[,1]==chroms[h],] #getting the p score of a specific chromosome
						ps_temp <- ps_chrom[ps_chrom[,3]>=quants_ps[h,j] & ps_chrom[,3]<quants_ps[h,j+1],] #here we get the p-scores in a range of a quant
						ps_range <- rbind(ps_range,ps_temp)
						
						#peaks data
						peaks_chrom <- peaks_dat[peaks_dat[,1]==chroms[h],] #getting the peaks of a specific chromosome
						peaks_temp <- peaks_chrom[peaks_chrom[,quant_col]>=quants_peaks[h,j] & peaks_chrom[,quant_col]<quants_peaks[h,j+1],c(1:3,value_col)] #here we get the peaks in a range of a quant
						peaks_range <- rbind(peaks_range,peaks_temp)
					}
					
					#creating contact bands for the p-score data of the range			
					quant <- quant_percentages[j]*100 #getting the specific quantile in percentage (no fraction)
					cat("\ncalculating the contact bands for the quantile range of ",quant,"% and ",quant_percentages[j+1]*100,"%\nplease wait...\n\n",sep="")
					#create contact bands (the 0 is because we are already taking the data in the range of the quantiles, so there is no need for the cutoff)			
					ps_conts <- contactBands_byChrom(ps_range,RE_gap,bp_gap,0,genome_sizes,rearranged_rawData)
					ps_conts <- ps_conts[[1]]
					
					#creating the bed files for the p-score and peaks data in the range			
					write.table(ps_conts,paste("~/Analyze4C/temp/ps_conts_",quant,"per.bed",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
					write.table(peaks_range,paste("~/Analyze4C/temp/peaks_",quant,"per.bed",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
					
					#creating the venn data frame for all but 0 percent
					if(z == 1)
					{
						#venn_DF[[venn_ind]] <- list(0)
						#venn_DF[[venn_ind]][[1]] <- as.character(paste(peaks_range[,1],peaks_range[,2],peaks_range[,3]))
						venn_DF[[venn_ind]] <- data.frame(rep(1,nrow(peaks_range)),matrix(0,nrow(peaks_range),numOFexperiments))
						venn_ind <- venn_ind + 1
					}
				}
			}
			
			#getting the quantile p-scores and peaks of the last quantile entered till the max p-scores. and then creating contact bands of this range.
			{
				ps_range <- c() #this will contain the p-score range data
				peaks_range <- c() #this will contain the peaks range data
				
				for(h in 1:length(chroms))
				{
					#p-score data		
					ps_chrom <- ps_dat[ps_dat[,1]==chroms[h],] #getting the p score of a specific chromosome
					ps_temp <- ps_chrom[ps_chrom[,3]>=quants_ps[h,numOfquants] & ps_chrom[,3]<=maxes_ps[h],] #here we get the p-scores in a range of a quant
					ps_range <- rbind(ps_range,ps_temp)
					
					#peaks data
					peaks_chrom <- peaks_dat[peaks_dat[,1]==chroms[h],] #getting the peaks of a specific chromosome
					peaks_temp <- peaks_chrom[peaks_chrom[,quant_col]>=quants_peaks[h,numOfquants] & peaks_chrom[,quant_col]<=maxes_peaks[h],c(1:3,value_col)] #here we get the peaks in a range of a quant
					peaks_range <- rbind(peaks_range,peaks_temp)
				}

				#creating contact bands for the p-score data of the range					
				quant <- quant_percentages[numOfquants]*100 #getting the specific quantile in percentage (no fraction)
				cat("\ncalculating the contact bands for the quantile range of ",quant,"% and 100%\nplease wait...\n\n",sep="")
				#create contact bands (the 0 is because we are already taking the data in the range of the quantiles, so there is no need for the cutoff)
				ps_conts <- contactBands_byChrom(ps_range,RE_gap,bp_gap,0,genome_sizes,rearranged_rawData)
				ps_conts <- ps_conts[[1]]
				
				#creating the bed files for the p-score and peaks data in the last range					
				write.table(ps_conts,paste("~/Analyze4C/temp/ps_conts_",quant,"per.bed",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
				write.table(peaks_range,paste("~/Analyze4C/temp/peaks_",quant,"per.bed",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
				
				#creating the venn data frame for the last quantile percent
				if(z == 1)
				{
					#venn_DF[[venn_ind]] <- list(0)
					#venn_DF[[venn_ind]][[1]] <- as.character(paste(peaks_range[,1],peaks_range[,2],peaks_range[,3]))
					venn_DF[[venn_ind]] <- data.frame(rep(1,nrow(peaks_range)),matrix(0,nrow(peaks_range),numOFexperiments))					
				}				
			}		
		}
		
		#intersecting and comparing the data from the contacts and peaks ranges of the same quantile
		{					
			#intersecting the rest of the quantile ranges
			intersections <- list() #this will contain all the intersections. each member of list is a different quantile
			peaks_self <- list() #this will contain all the self of peaks. each member of list is a different quantile
			ps_self <- list() #this will contain all the self of the p-scores. each member of list is a different quantile
			intersections_names <- c() #this will contain the numbers of the quantiles which will eventually be the names of the members of the list

			#starting with the range of 0% till the first quantile
			system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks_0per.bed -b ~/Analyze4C/temp/ps_conts_0per.bed -wo > ~/Analyze4C/temp/peaks_VS_contacts_0per.bed",sep="")) #maybe i should put the intersection files in a designated folder					
			
			#getting the sum of self
			system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks_0per.bed -b ~/Analyze4C/temp/peaks_0per.bed -wo > ~/Analyze4C/temp/peaks_self_0per.bed",sep=""))	
			system(paste("bedtools intersect -a ~/Analyze4C/temp/ps_conts_0per.bed -b ~/Analyze4C/temp/ps_conts_0per.bed -wo > ~/Analyze4C/temp/contacts_self_0per.bed",sep=""))					
			
			#importing the data from the intersection files			
			if(file.info("~/Analyze4C/temp/peaks_VS_contacts_0per.bed")$size != 0)
			{
				intersections[[1]] <- read.table("~/Analyze4C/temp/peaks_VS_contacts_0per.bed")
			}
			else
			{
				intersections[[1]] <- data.frame(0,0,0,0,0,0,0,0)
			}
			
			if(file.info("~/Analyze4C/temp/peaks_self_0per.bed")$size != 0)
			{
				peaks_self[[1]] <- read.table("~/Analyze4C/temp/peaks_self_0per.bed")
			}
			else
			{
				peaks_self[[1]] <- data.frame(0,0,0,0,0,0,0,0)
			}
			
			if(file.info("~/Analyze4C/temp/contacts_self_0per.bed")$size != 0)
			{
				ps_self[[1]] <- read.table("~/Analyze4C/temp/contacts_self_0per.bed")
			}
			else
			{
				ps_self[[1]] <- data.frame(0,0,0,0,0,0,0,0)
			}
			
			quant_name <- "0%"
			intersections_names <- c(intersections_names,quant_name)

			if(z == 1)
			{
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
			}
			
			for(n in 1:numOfquants)
			{
				quant <- quant_percentages[n]*100 #getting the specific quantile in percentage (no fraction)
				system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks_",quant,"per.bed -b ~/Analyze4C/temp/ps_conts_",quant,"per.bed -wo > ~/Analyze4C/temp/peaks_VS_contacts_",quant,"per.bed",sep="")) #maybe i should put the intersection files in a designated folder

				#getting the sum of self
				system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks_",quant,"per.bed -b ~/Analyze4C/temp/peaks_",quant,"per.bed -wo > ~/Analyze4C/temp/peaks_self_",quant,"per.bed",sep=""))	
				system(paste("bedtools intersect -a ~/Analyze4C/temp/ps_conts_",quant,"per.bed -b ~/Analyze4C/temp/ps_conts_",quant,"per.bed -wo > ~/Analyze4C/temp/contacts_self_",quant,"per.bed",sep=""))					
				
				#importing the data from the intersection files			
				if(file.info(paste("~/Analyze4C/temp/peaks_VS_contacts_",quant,"per.bed",sep=""))$size != 0)
				{
					intersections[[n+1]] <- read.table(paste("~/Analyze4C/temp/peaks_VS_contacts_",quant,"per.bed",sep=""))
				}
				else
				{
					intersections[[n+1]] <- data.frame(0,0,0,0,0,0,0,0)
				}
				
				if(file.info(paste("~/Analyze4C/temp/peaks_self_",quant,"per.bed",sep=""))$size != 0)
				{
					peaks_self[[n+1]] <- read.table(paste("~/Analyze4C/temp/peaks_self_",quant,"per.bed",sep=""))
				}
				else
				{
					peaks_self[[n+1]] <- data.frame(0,0,0,0,0,0,0,0)
				}
				
				if(file.info(paste("~/Analyze4C/temp/contacts_self_",quant,"per.bed",sep=""))$size != 0)
				{
					ps_self[[n+1]] <- read.table(paste("~/Analyze4C/temp/contacts_self_",quant,"per.bed",sep=""))
				}
				else
				{
					ps_self[[n+1]] <- data.frame(0,0,0,0,0,0,0,0)
				}
				
				quant_name <- paste(quant,"%",sep="")
				intersections_names <- c(intersections_names,quant_name)
			}
			names(intersections) <- intersections_names #giving each member of list the corresponding quantile
			names(peaks_self) <- intersections_names #giving each member of list the corresponding quantile
			names(ps_self) <- intersections_names #giving each member of list the corresponding quantile
			
			#ask if to sum the values per whole genome, per trans, or per chromosome
			if(z == 1)
			{
				cat("\nwhat actions would you like to do with the intersections data?\n")
				choice4 <- readline(prompt=cat("\n1) sum the data\n2) other actions\n\n")) #i need to add here other actions i could do
			}

			if(choice4 == 1) #summing the data
			{
				summer <- matrix(NA,nrow=1,ncol=(numOfquants+1))
				summer_peaks_self <- matrix(NA,nrow=1,ncol=(numOfquants+1))
				summer_ps_self <- matrix(NA,nrow=1,ncol=(numOfquants+1))
				
				colnames(summer) <- intersections_names
				colnames(summer_peaks_self) <- intersections_names
				colnames(summer_ps_self) <- intersections_names
				
				if(z == 1)
				{
					choice5 <- as.integer(readline(prompt=cat("\nwhat chromosomes should be summed? (enter the option number)\n1) whole genome\n2) trans\n3) cis\n4) each chromosome separately\n\n")))
				}
				
				if(choice5 == 1) #whole genome
				{
					for(g in 1:(numOfquants+1))
					{
						summer[g] <- sum(intersections[[g]][8])
						summer_peaks_self[g] <- sum(peaks_self[[g]][9])
						summer_ps_self[g] <- sum(ps_self[[g]][7])
					}
					summed_chr <- "all" #this will be entered later in 'ChIPseqVScontacts_sumOFintersections_plots.txt'
				}
				else if(choice5 == 2) #trans
				{
					for(g in 1:(numOfquants+1))
					{
						summer[g] <- sum(intersections[[g]][intersections[[g]][,1]!=cis,8])
						summer_peaks_self[g] <- sum(peaks_self[[g]][peaks_self[[g]][,1]!=cis,9])
						summer_ps_self[g] <- sum(ps_self[[g]][ps_self[[g]][,1]!=cis,7])
					}
					summed_chr <- "trans" #this will be entered later in 'ChIPseqVScontacts_sumOFintersections_plots.txt'							
				}
				else if(choice5 == 3) #cis
				{
					for(g in 1:(numOfquants+1))
					{
						summer[g] <- sum(intersections[[g]][intersections[[g]][,1]==cis,8])
						summer_peaks_self[g] <- sum(peaks_self[[g]][peaks_self[[g]][,1]==cis,9])
						summer_ps_self[g] <- sum(ps_self[[g]][ps_self[[g]][,1]==cis,7])
					}
					summed_chr <- "cis" #this will be entered later in 'ChIPseqVScontacts_sumOFintersections_plots.txt'		
				}
				else if(choice5 == 4) #each chromosome separately
				{
					summer <- matrix(NA,nrow=numOFchroms,ncol=(numOfquants+1))
					summer_peaks_self <- matrix(NA,nrow=numOFchroms,ncol=(numOfquants+1))
					summer_ps_self <- matrix(NA,nrow=numOFchroms,ncol=(numOfquants+1))

					for(d in 1:numOFchroms)
					{
						for(g in 1:(numOfquants+1))
						{
							summer[d,g] <- sum(intersections[[g]][intersections[[g]][,1]==d,8])
							summer_peaks_self[d,g] <- sum(peaks_self[[g]][peaks_self[[g]][,1]==d,9])
							summer_ps_self[d,g] <- sum(ps_self[[g]][ps_self[[g]][,1]==d,7])
						}
					}
					summer <- data.frame(summer)
					colnames(summer) <- intersections_names
					rownames(summer) <- c(1:numOFchroms)						
					print(summer)
					
					summer_peaks_self <- data.frame(summer_peaks_self)
					colnames(summer_peaks_self) <- intersections_names
					rownames(summer_peaks_self) <- c(1:numOFchroms)						
					print(summer_peaks_self)
					
					summer_ps_self <- data.frame(summer_ps_self)
					colnames(summer_ps_self) <- intersections_names
					rownames(summer_ps_self) <- c(1:numOFchroms)						
					print(summer_ps_self)
					
					summed_chr <- "each separate" #this will be entered later in 'ChIPseqVScontacts_sumOFintersections_plots.txt'
				}
				
				summer_all[[z]] <- summer
				summer_peaks_self_all[[z]] <- summer_peaks_self
				summer_ps_self_all[[z]] <- summer_ps_self
			}
			else if(choice4 == 2) #i need to fill this with other actions
			{
			
			}	
		}

	#filling in the data for the venn diagrams	
	{
		if(quant_percentages[1] != 0)
		{
			venn_quants <- c(0,quant_percentages)
		}
		else
		{
			venn_quants <- quant_percentages
		}

		ind_count <- 0
		for(m in (venn_quants)*100)
		{	
			ind_count <- ind_count + 1
			if(venn_ans1 == "y")
			{
				system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks_",m,"per.bed -b ~/Analyze4C/temp/ps_conts_",m,"per.bed -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_",m,"per_A.bed",sep=""))
				venn_A <- read.table(paste("~/Analyze4C/temp/venn_",m,"per_A.bed",sep=""))
				if(venn_ans2 == "y")
				{
					system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks_",m,"per.bed -b ~/Analyze4C/temp/ps_conts_",m,"per.bed -c > ~/Analyze4C/temp/venn_",m,"per_B.bed",sep=""))
					venn_B <- read.table(paste("~/Analyze4C/temp/venn_",m,"per_B.bed",sep=""))
					#if the number in the 1st is 0 but the number in the 2nd is more than ovlp_num then we will consider it to be 1
					venn_A[venn_B[,5]>=ovlp_num,5] <- 1
				}
				#add the correct number to venn_col_num rows
#i need to collect the data differently into venn_DF, collect the actual indices that have 1, and turn them into a one string, and then list them together (instead of columns in a data frame)
# and then use calculate.overlap
				#venn_DF[[ind_count]][[venn_col_num]] <- as.character(paste(venn_A[venn_A[,5] == 1,1:3][,1],venn_A[venn_A[,5] == 1,1:3][,2],venn_A[venn_A[,5] == 1,1:3][,3]))			
				venn_DF[[ind_count]][venn_A[,5]>=1,venn_col_num] <- 1
			}
			else
			{
				system(paste("bedtools intersect -a ~/Analyze4C/temp/peaks_",m,"per.bed -b ~/Analyze4C/temp/ps_conts_",m,"per.bed -c > ~/Analyze4C/temp/venn_",m,"per.bed",sep=""))
				venn <- read.table(paste("~/Analyze4C/temp/venn_",m,"per.bed",sep=""))
				#venn_DF[[ind_count]][[venn_col_num]] <- as.character(paste(venn[venn[,5] == 1,1:3][,1],venn[venn[,5] == 1,1:3][,2],venn[venn[,5] == 1,1:3][,3]))			
				#add the correct number to venn_col_num rows
				venn_DF[[ind_count]][venn[,5]>=1,venn_col_num] <- 1	
			}			
		}
	}
		
		#removing the files from 'temp' folder
		{
			#removing first the range that begins with 0%
			#system(paste("rm ~/Analyze4C/temp/ps_conts_0per.bed",sep=""))
			#system(paste("rm ~/Analyze4C/temp/peaks_0per.bed",sep=""))
			#system(paste("rm ~/Analyze4C/temp/peaks_VS_contacts_0per.bed",sep=""))
			
			#removing the rest of the range files
			#for(m in 1:numOfquants)
			#{
			#	quant <- quant_percentages[m]*100 #getting the specific quantile in percentage (no fraction)
			#	system(paste("rm ~/Analyze4C/temp/ps_conts_",quant,"per.bed",sep=""))
			#	system(paste("rm ~/Analyze4C/temp/peaks_",quant,"per.bed",sep=""))
			#	#system(paste("rm ~/Analyze4C/temp/peaks_VS_contacts_",quant,"per.bed",sep=""))
			#}
			
			#removing the intersection files
			#system(paste("rm ~/Analyze4C/temp/peaks_VS_contacts_0per.bed",sep=""))
			#for(n in 1:numOfquants)
			#{
			#	quant <- quant_percentages[n]*100 #getting the specific quantile in percentage (no fraction)
			#	system(paste("rm ~/Analyze4C/temp/peaks_VS_contacts_",quant,"per.bed",sep=""))
			#}
			
			#system("rm ~/Analyze4C/temp/ps_dat_temp.sgr")
			system("rm ~/Analyze4C/temp/*")
		}	
		#i might want to take the last column of the intersection files (column 7, these are the number of bp that intersect), and sum them.
		#then compare this number and the number we get with other experiments
	}		

	names(summer_all) <- all_file_names #adding the names of the files to the list, each 'summer' will get the name of the corresponding file
	names(summer_peaks_self_all) <- all_file_names #adding the names of the files to the list, each 'summer_peaks_self_all' will get the name of the corresponding file
	names(summer_ps_self_all) <- all_file_names #adding the names of the files to the list, each 'summer_ps_self_all' will get the name of the corresponding file

	#create plots:		
	peaks_dataframe <- c()
	ps_dataframe <- c()
	exp_names <- c() #contains the names given by user to each example
	
	for(h in 1:z)
	{
		exp_name <- readline(prompt=cat("\nenter the name you would like to use in the graph for experiment",h,":\n"))
		exp_names <- c(exp_names,exp_name)
		peaks_per <- summer_all[[h]]/summer_peaks_self_all[[h]]
		ps_per <- summer_all[[h]]/summer_ps_self_all[[h]]
		
		percentages_peaks <- c()
		for(a in 1:(length(colnames(peaks_per))-1))
		{
			percentages_peaks_temp <- paste(colnames(peaks_per)[a],"-",colnames(peaks_per)[a+1],sep="")
			percentages_peaks <- c(percentages_peaks,percentages_peaks_temp)
		}
		percentages_peaks_temp <- paste(colnames(peaks_per)[length(colnames(peaks_per))],"-100%",sep="")
		percentages_peaks <- c(percentages_peaks,percentages_peaks_temp)
		
		peaks_dataframe_temp <- data.frame(percentages_peaks,as.numeric(t(peaks_per)),rep(exp_name,nrow(peaks_per)))
		names(peaks_dataframe_temp) <- c("quantile_range","per","experiment")
		peaks_dataframe <- rbind(peaks_dataframe,peaks_dataframe_temp)
		
		percentages_ps <- c()
		for(a in 1:(length(colnames(ps_per))-1))
		{
			percentages_ps_temp <- paste(colnames(ps_per)[a],"-",colnames(ps_per)[a+1],sep="")
			percentages_ps <- c(percentages_ps,percentages_ps_temp)
		}
		percentages_ps_temp <- paste(colnames(ps_per)[length(colnames(ps_per))],"-100%",sep="")
		percentages_ps <- c(percentages_ps,percentages_ps_temp)
		
		ps_dataframe_temp <- data.frame(percentages_ps,as.numeric(t(ps_per)),rep(exp_name,nrow(ps_per)))
		names(ps_dataframe_temp) <- c("quantile_range","per","experiment")
		ps_dataframe <- rbind(ps_dataframe,ps_dataframe_temp)
	}

	#getting the date and time
	DandT1 <- toString(Sys.time())
	DandT2 <- gsub(" ","_",DandT1)
	DandT2 <- gsub(":","",DandT2)

	#creating a plot
	peaks_plot <- ggplot2::ggplot(peaks_dataframe, ggplot2::aes(x=quantile_range,y=per,fill=experiment)) + ggplot2::geom_bar(position = "dodge",stat = "identity") + ggplot2::ggtitle("ratio of contact intersections with peaks and all peaks") + ggplot2::labs(x="quantile (%)",y="intersections(bp)/peaks(bp)") + ggplot2::theme(plot.title = ggplot2::element_text(size=15))	
	if(numOFexperiments == 2) #this comparison will only work if there are 2 experiments, if there is a different number then we can plot the data but not show significance
	{	
		ORs1 <- c() #combines all odd ratio of this calculation
		for(e in 1:(numOfquants+1))
		{
			cat("\n****************************************************\n\n")
			cat("\nthe range of",percentages_peaks[e],":")
			if(summer_all[[1]][e]>0 & summer_all[[2]][e]>0)
			{			
				#testing the significance of difference between percentages of both data groups
				cat("\n\ntesting the significance of the difference between the percentages of bp intersected out of the bp of the peaks bands, for both contact bands groups...\n\n")
				cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the peaks file\nthe 2nd creates a table of all the successes (bp of peaks intersections) and fails (all the bp of the peaks bands that didn't intersect) and uses them\n\n")					
				prop_peaks1 <- prop.test(c(summer_all[[1]][e],summer_all[[2]][e]),c(summer_peaks_self_all[[1]][e],summer_peaks_self_all[[2]][e]),correct=FALSE)
				print(prop_peaks1)
				prop_peaks2 <- prop.test(matrix(c(summer_peaks_self_all[[1]][e]-summer_all[[1]][e],summer_all[[1]][e],summer_peaks_self_all[[2]][e]-summer_all[[2]][e],summer_all[[2]][e]),ncol=2),correct=FALSE)
				print(prop_peaks2)
				prop_ans1 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
				if(prop_ans1 == 1)
				{
					prop_peaks <- prop_peaks1
				}
				else
				{
					prop_peaks <- prop_peaks2
				}
				
				#creating a barplot with significance stars	
				label.df <- data.frame(quantile_range = percentages_peaks[e],experiment = peaks_dataframe[peaks_dataframe$quantile_range == percentages_peaks[e],]$experiment[which.max(peaks_dataframe[peaks_dataframe$quantile_range == percentages_peaks[e],]$per)],per = peaks_dataframe[peaks_dataframe$quantile_range == percentages_peaks[e],]$per[which.max(peaks_dataframe[peaks_dataframe$quantile_range == percentages_peaks[e],]$per)]+(0.1*peaks_dataframe[peaks_dataframe$quantile_range == percentages_peaks[e],]$per[which.max(peaks_dataframe[peaks_dataframe$quantile_range == percentages_peaks[e],]$per)]))
				if(prop_peaks$p.value <= 0.0001)
				{
					peaks_plot <- peaks_plot + ggplot2::geom_text(data = label.df, label = "****",size=8)
				}
				else if(prop_peaks$p.value <= 0.001)
				{
					peaks_plot <- peaks_plot + ggplot2::geom_text(data = label.df, label = "***",size=8)
				}
				else if(prop_peaks$p.value <= 0.01)
				{
					peaks_plot <- peaks_plot + ggplot2::geom_text(data = label.df, label = "**",size=8)
				}
				else if(prop_peaks$p.value <= 0.05)
				{
					peaks_plot <- peaks_plot + ggplot2::geom_text(data = label.df, label = "*",size=8)
				}
			}
			else
			{
				if(summer_all[[1]][e] == 0)
				{
					cat("\nthe",exp_names[1],"data does not have any intersections in the range of",percentages_peaks[e],"\n\n")
				}
				else if(summer_all[[2]][e] == 0)
				{
					cat("\nthe",exp_names[2],"data does not have any intersections in the range of",percentages_peaks[e],"\n\n")
				}
				else if(summer_all[[1]][e] == 0 & summer_all[[2]][e] == 0)
				{
					cat("\nthe",exp_names[1],"and",exp_names[2],"data do not have any intersections in the range of",percentages_peaks[e],"\n\n")
				}			
			}
			
			#calculating effect sizes
			cat("\neffect sizes:\n")
			OR1 <- compute.es::propes(peaks_dataframe[peaks_dataframe$experiment==exp_names[1],2][e],peaks_dataframe[peaks_dataframe$experiment==exp_names[2],2][e],summer_peaks_self_all[[1]][1],summer_peaks_self_all[[2]][1])
			cat("\n")
			ORs1 <- c(ORs1,OR1$OR)
			
			#phi
			cat("\nphi:\n")
			print(psych::phi(matrix(c(peaks_dataframe[peaks_dataframe$experiment==exp_names[1],2][e],summer_peaks_self_all[[1]][1]-peaks_dataframe[peaks_dataframe$experiment==exp_names[1],2][e],peaks_dataframe[peaks_dataframe$experiment==exp_names[2],2][e],summer_peaks_self_all[[2]][1]-peaks_dataframe[peaks_dataframe$experiment==exp_names[2],2][e]),ncol=2,byrow=TRUE)))
			cat("\n")
		}
		names(ORs1) <- percentages_ps
		cat("All odds ratios:\n")
		print(ORs1)
		cat("\n\n")		
	}
	
	cat("\nplotting the peaks bands bp percentage data\nplease wait...\n\n")	
	print(peaks_plot)
	choice6 <- readline("\nwould you like to save these plots?\ny/n\n\n")
	if(choice6 == "y")
	{	
		wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
		ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))	
		peaks_plot_name <- paste("~/Analyze4C/plots/expressionVScontacts_sumOFintersections_peaks_plot_",DandT2,".png",sep="")
		ggplot2::ggsave(peaks_plot_name,width=wd,height=ht)
	}
	
	#creating a plot for the p-score data
	ps_plot <- 	ggplot2::ggplot(ps_dataframe, ggplot2::aes(x=quantile_range,y=per,fill=experiment)) + ggplot2::geom_bar(position = "dodge",stat = "identity") + ggplot2::ggtitle("ratio of contact intersections with peaks and all contacts") + ggplot2::labs(x="quantile (%)",y="intersections(bp)/contacts(bp)") + ggplot2::theme(plot.title = ggplot2::element_text(size=15))	
	if(numOFexperiments == 2) #this comparison will only work if there are 2 experiments, if there is a different number then we can plot the data but not show significance
	{	
		ORs2 <- c() #combines all odd ratio of this calculation
		for(e in 1:(numOfquants+1))
		{
			cat("\n****************************************************\n\n")
			cat("\nthe range of",percentages_ps[e],":")
			if(summer_all[[1]][e]>0 & summer_all[[2]][e]>0)
			{			
				#testing the significance of difference between percentages of both data groups
				cat("\n\ntesting the significance of the difference between the percentages of bp intersected out of the bp of the contact bands, for both contact bands groups...\n\n")
				cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the contact bands file\nthe 2nd creates a table of all the successes (bp of contacts intersections) and fails (all the bp of the contact bands that didn't intersect) and uses them\n\n")					
				prop_ps1 <- prop.test(c(summer_all[[1]][e],summer_all[[2]][e]),c(summer_ps_self_all[[1]][e],summer_ps_self_all[[2]][e]),correct=FALSE)
				print(prop_ps1)
				prop_ps2 <- prop.test(matrix(c(summer_ps_self_all[[1]][e]-summer_all[[1]][e],summer_all[[1]][e],summer_ps_self_all[[2]][e]-summer_all[[2]][e],summer_all[[2]][e]),ncol=2),correct=FALSE)
				print(prop_ps2)
				prop_ans1 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
				if(prop_ans1 == 1)
				{
					prop_ps <- prop_ps1
				}
				else
				{
					prop_ps <- prop_ps2
				}
				
				#creating a barplot with significance stars	
				label.df <- data.frame(quantile_range = percentages_ps[e],experiment = ps_dataframe[ps_dataframe$quantile_range == percentages_ps[e],]$experiment[which.max(ps_dataframe[ps_dataframe$quantile_range == percentages_ps[e],]$per)],per = ps_dataframe[ps_dataframe$quantile_range == percentages_ps[e],]$per[which.max(ps_dataframe[ps_dataframe$quantile_range == percentages_ps[e],]$per)]+(0.1*ps_dataframe[ps_dataframe$quantile_range == percentages_ps[e],]$per[which.max(ps_dataframe[ps_dataframe$quantile_range == percentages_ps[e],]$per)]))
				if(prop_ps$p.value <= 0.0001)
				{
					ps_plot <- ps_plot + ggplot2::geom_text(data = label.df, label = "****",size=8)
				}
				else if(prop_ps$p.value <= 0.001)
				{
					ps_plot <- ps_plot + ggplot2::geom_text(data = label.df, label = "***",size=8)
				}
				else if(prop_ps$p.value <= 0.01)
				{
					ps_plot <- ps_plot + ggplot2::geom_text(data = label.df, label = "**",size=8)
				}
				else if(prop_ps$p.value <= 0.05)
				{
					ps_plot <- ps_plot + ggplot2::geom_text(data = label.df, label = "*",size=8)
				}
			}
			else
			{
				if(summer_all[[1]][e] == 0)
				{
					cat("\nthe",exp_names[1],"data does not have any intersections in the range of",percentages_ps[e],"\n\n")
				}
				else if(summer_all[[2]][e] == 0)
				{
					cat("\nthe",exp_names[2],"data does not have any intersections in the range of",percentages_ps[e],"\n\n")
				}
				else if(summer_all[[1]][e] == 0 & summer_all[[2]][e] == 0)
				{
					cat("\nthe",exp_names[1],"and",exp_names[2],"data do not have any intersections in the range of",percentages_ps[e],"\n\n")
				}			
			}

			#calculating effect sizes
			cat("\neffect sizes:\n")
			OR2 <- compute.es::propes(ps_dataframe[ps_dataframe$experiment==exp_names[1],2][e],ps_dataframe[ps_dataframe$experiment==exp_names[2],2][e],summer_ps_self_all[[1]][1],summer_ps_self_all[[2]][1])
			cat("\n")
			ORs2 <- c(ORs2,OR2$OR)

			#phi
			cat("\nphi:\n")
			print(psych::phi(matrix(c(ps_dataframe[ps_dataframe$experiment==exp_names[1],2][e],summer_ps_self_all[[1]][1]-ps_dataframe[ps_dataframe$experiment==exp_names[1],2][e],ps_dataframe[ps_dataframe$experiment==exp_names[2],2][e],summer_ps_self_all[[2]][1]-ps_dataframe[ps_dataframe$experiment==exp_names[2],2][e]),ncol=2,byrow=TRUE)))
			cat("\n")			
		}
		names(ORs2) <- percentages_ps	
		cat("All odds ratios:\n")
		print(ORs2)
		cat("\n\n")
	}
	
	cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")
	print(ps_plot)
	choice7 <- readline("\nwould you like to save these plots?\ny/n\n\n")
	if(choice7 == "y")
	{
		wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
		ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))		
		ps_plot_name <- paste("~/Analyze4C/plots/expressionVScontacts_sumOFintersections_pScore_plot_",DandT2,".png",sep="")
		ggplot2::ggsave(ps_plot_name,width=wd,height=ht)		
	}

	# create venn diagram using the data frame venn_DF[[ind_count]]	
	venn_flag <- 0 #if venn_flag is 1 then we save a diagram and it should be recorded in 'ChIPseqVScontacts_sumOFintersections_plots'
	for(w in 1:length(venn_quants))
	{
		#depending on how many files numOFexperiments is, the venn diagram function will be different
		cat("\nvenn diagram for quantile",as.numeric(venn_quants[w]),":\n\n")
		vennName <- paste("venn_quantiles_",as.numeric(venn_quants[w])*100,"per_",DandT2,".jpg",sep="")
		venn_flag <- venn_creator(venn_DF[[w]],numOFexperiments,1,vennName)			
	}
	
	#recording the data into 'ChIPseqVScontacts_sumOFintersections_plots.txt'
	if(choice6 == "y" | choice7 == "y" | venn_flag == 1)
	{
		#saving the parameters and details to ChIPseqVScontacts_sumOFintersections_plots.txt
		ChIPseqVScontacts_sumOFintersections_plots[nrow(ChIPseqVScontacts_sumOFintersections_plots)+1,] <- c(DandT1,file.name_peaks,numOFexperiments,exps,int_chroms,summed_chr,RE_gap,bp_gap,quant_source,value_source)
		#sorting the list of experiments by bait alphabetically (and sorting the row indices)
		ChIPseqVScontacts_sumOFintersections_plots <- ChIPseqVScontacts_sumOFintersections_plots[order(ChIPseqVScontacts_sumOFintersections_plots$Date_and_Time),]
		rownames(ChIPseqVScontacts_sumOFintersections_plots) <- seq(length=nrow(ChIPseqVScontacts_sumOFintersections_plots))
		#adding the new data to the file (by erasing the old one and creating a new one)
		system("rm ChIPseqVScontacts_sumOFintersections_plots.txt")
		write.table(ChIPseqVScontacts_sumOFintersections_plots,"ChIPseqVScontacts_sumOFintersections_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
	}
	choice6 <- "n"
	choice7 <- "n"
	venn_flag <- 0	
	#intersect the peaks of a quartile and the contact bands of the same quartile
	#the methods of intersection could be: 1) sum of all intersected (whole genome, cis, trans, or by chromosome). 2) sum all bps of intersections (whole genome, cis, trans, or by chromosome). 3) create windows and do 1 or 2, or count the reads of the RE sites where the contacts are and the peaks.
}
