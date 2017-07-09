#' @export

#comparing FPKM data with the number p-scores of intersected contacts

rnaSeqCorContacts <- function(Experiments_4C,rearranged_rawData,expressionVScontacts_correlation_plots)
{
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

	translocated_flag <- 0
	files_flag <- 0 #this marks if the user does or doesn't want to choose more files to use (0 choose more, 1 don't choose more)
	z <- 1 #counts how many files are chosen
	secs_rmv <- "no" #this says if there were any sections removed
	CO_ps_type <- c()
	CO_ps <- c()
	CO_ps_chroms <- c()
	CO_FPKM_type <- c()
	CO_FPKM <- c()
	CO_FPKM_chroms <- c()
	wind_flag <- 0 #this marks if the user wants to change the window sizes
	by_winds <- "no"
	winds <- "NA"
	psORreads <- "pScores" #are p-score or reads correlated with the FPKM data
	save_flag <- 0 #marks if there was any plot saved, if yes then the details are entered into 'expressionVScontacts_correlation_plots'
		
	#the user decides if to correlate individual RE sites or by windows
	choice4 <- as.integer(readline(prompt=cat("\nchoose a method of comparison:\n\n1) intersect FPKM bands with contact bands\n2) create windows, add the values of FPKM and contacts, and correlate them\n\n")))
	if(choice4 == 2)
	{
		by_winds <- "yes"
		winds <- c()
	}
	
	#asking if to compare multiple FPKMs against one p-score or one FPKM file against multiple p-score files
	choice1 <- as.integer(readline(prompt=cat("\nchoose how you would like to compare the files:\n1) compare one FPKM file against multiple p-score files\n2) compare one p-score file against multiple FPKM files\n\n")))
	if(choice1 == 1)
	{
	#choosing one FPKM file
	#choosing multiple p-scores
	
		#getting the FPKM file	
		ls_files_FPKM <- system("ls ~/Analyze4C/RNAseq/FPKM | grep '.bed$'",intern=TRUE)
		repeat
		{	
			file.name_FPKM <- readline(prompt=cat("\nThese are the expression data (FPKM) files available\nwhich would you like to use?\n",ls_files_FPKM,"",sep="\n"))
			ind_files_FPKM <- pmatch(file.name_FPKM,ls_files_FPKM,nomatch=0)
			if(ind_files_FPKM == 0)
			{
				cat("no such file exists.\nplease try again.\n\n")
			}
			else
			{
				#asking if to not use all the FPKMs that are 0, this is important because when getting the quantiles it considers the 0 examples as well
				#and since FPKMs of 0 could be considered as if there is no expression at all in those areas, then we might be able to remove them
				#the user decided
				remZeroFPKM <- readline(prompt=cat("\nshould FPKM's of 0 be removed?\ny/n\n\n"))

				#importing the data from the file
				FPKM <- read.table(paste("~/Analyze4C/RNAseq/FPKM/",file.name_FPKM,sep=""))
				if(remZeroFPKM == "y") #removing all FPKMs that are 0 (since they are basically rna-seq reads that don't exist)
				{
					FPKM <- FPKM[FPKM[,4]>0,]
				}				
				break
			}

			#here i might need to edit and sort the data in FPKM
			#FPKM <- FPKM[FPKM[,1]!="Pt" & FPKM[,1]!="Mt",]
			#FPKM <- FPKM[order(FPKM[,1],FPKM[,2]),]
			#names(FPKM) <- c("chr","first","last","FPKM")
		}
		FPKM_Experiments <- file.name_FPKM
		num_of_FPKM_files <- 1
		
		ps_datName_all <- c() #this will contain all the names the user gives the p-score data for the plot
	}
	else
	{

	#choosing one p-score file
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
				#here we check to see if the p-score file chosen is of a translocation switched type
				if(identical(grep("translocatedSwitched",ps.file),integer(0)) == "FALSE") #the p-score file is of translocation switched and needs to be removed
				{
					#the translocation area needs to be removed both from the p-score and from the FPKM file
					translocated_flag <- 1
				}
				break
			}
		}
		ps_Experiments <- ps.file
		num_of_pScore_files <- 1
		
		#importing the p-score data from the file 	
		ps <- read.table(paste("~/Analyze4C/pScores/",ps_folder,"/",ps.file,sep=""))

		#finding the specific p-score files information in Experiments_4C
		sp1 <- unlist(strsplit(ps.file,"[.]"))
		sp2 <- unlist(strsplit(sp1[1],"_"))
		
		#getting the size of the window per RE sites
		sp3 <- sp2[as.integer(grep("RE$",sp2))]
		RE_wind <- as.numeric(strsplit(sp3,"RE"))
		
		if(choice3 == 2)
		{
			#getting the size of the window per bp sites
			sp4 <- sp2[as.integer(grep("bp$",sp2))]
			bp_wind <- as.numeric(strsplit(sp4,"bp"))
		}

		#finding the specific p-scores files information in Experiments_4C			
		out <- findIn.Experiments_4C(sp2[1],sp2[2],sp2[3],Experiments_4C)
		
		#getting the cis chromosome number
		cis <- as.numeric(out[2])

		FPKM_datName_all <- c() #this will contain all the names the user gives the FPKM data for the plot	
	}
	
	intersected_pairs_all <- c() #this will contain all the intersections data from all the files
	intersected_reads_pairs_all <- c() #this will contain all the intersections data from all the files, for the reads

	while(files_flag == 0)
	{
		if(choice1 == 1)
		{
			cat("\np-score file number ",z,":\n",sep="")
			#choosing a p-score file:
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
					#here we check to see if the p-score file chosen is of a translocation switched type
					if(identical(grep("translocatedSwitched",ps.file),integer(0)) == "FALSE") #the p-score file is of translocation switched and needs to be removed
					{
						#the translocation area needs to be removed both from the p-score and from the FPKM file
						translocated_flag <- 1
					}
					break
				}
			}	
				
			#gets the names of the experiments to later be entered into the 'expressionVScontacts_correlation_plots' file
			if(z==1)
			{
				ps_Experiments <- ps.file
			}
			else
			{
				ps_Experiments <- paste(ps_Experiments," ; ",ps.file,sep="")
			}

			#importing the p-score data from the file 	
			ps <- read.table(paste("~/Analyze4C/pScores/",ps_folder,"/",ps.file,sep=""))

			#finding the specific p-scores files information in Experiments_4C
			sp1 <- unlist(strsplit(ps.file,"[.]"))
			sp2 <- unlist(strsplit(sp1[1],"_"))
			
			#getting the size of the window per RE sites
			sp3 <- sp2[as.integer(grep("RE$",sp2))]
			RE_wind <- as.numeric(strsplit(sp3,"RE"))
			
			if(choice3 == 2)
			{
				#getting the size of the window per bp sites
				sp4 <- sp2[as.integer(grep("bp$",sp2))]
				bp_wind <- as.numeric(strsplit(sp4,"bp"))
			}

			#finding the specific p-scores files information in Experiments_4C			
			out <- findIn.Experiments_4C(sp2[1],sp2[2],sp2[3],Experiments_4C)
			
			#getting the cis chromosome number
			cis <- as.numeric(out[2])
		}
		else
		{
			#choosing an FPKM file
			cat("\nFPKM file number ",z,":\n",sep="")
			ls_files_FPKM <- system("ls ~/Analyze4C/RNAseq/FPKM | grep '.bed$'",intern=TRUE)
			repeat
			{	
				file.name_FPKM <- readline(prompt=cat("\nThese are the expression data (FPKM) files available\nwhich would you like to use?\n",ls_files_FPKM,"",sep="\n"))
				ind_files_FPKM <- pmatch(file.name_FPKM,ls_files_FPKM,nomatch=0)
				if(ind_files_FPKM == 0)
				{
					cat("no such file exists.\nplease try again.\n\n")
				}
				else
				{
					#asking if to not use all the FPKMs that are 0, this is important because when getting the quantiles it considers the 0 examples as well
					#and since FPKMs of 0 could be considered as if there is no expression at all in those areas, then we might be able to remove them
					#the user decided
					remZeroFPKM <- readline(prompt=cat("\nshould FPKM's of 0 be removed?\ny/n\n\n"))
				
					#importing the data from the file
					FPKM <- read.table(paste("~/Analyze4C/RNAseq/FPKM/",file.name_FPKM,sep=""))
					if(remZeroFPKM == "y") #removing all FPKMs that are 0 (since they are basically rna-seq reads that don't exist)
					{
						FPKM <- FPKM[FPKM[,4]>0,]
					}					
					break
				}

				#here i might need to edit and sort the data in FPKM
				#FPKM <- FPKM[FPKM[,1]!="Pt" & FPKM[,1]!="Mt",]
				#FPKM <- FPKM[order(FPKM[,1],FPKM[,2]),]
				#names(FPKM) <- c("chr","first","last","FPKM")
			}	
			
			#gets the names of the experiments to later be entered into the 'expressionVScontacts_correlation_plots' file
			if(z==1)
			{
				FPKM_Experiments <- file.name_FPKM
			}
			else
			{
				FPKM_Experiments <- paste(FPKM_Experiments," ; ",file.name_FPKM,sep="")
			}		
		}

		if(choice4 == 2 & wind_flag == 0)
		{
			genome_sizes_temp <- as.data.frame(cbind(gsub("^chr","",genome_sizes[,1]),genome_sizes[,2])) #removing 'chr' from the chromosome numbers, if there are any with 'chr'
			write.table(genome_sizes_temp,"~/Analyze4C/temp/genome_sizes_temp.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t") #creating a temp file of the genome sizes			
			winSize <- as.integer(readline(prompt=cat("\nplease enter a window size (bp) that the whole genome will divided into:\n\n")))
			system(paste("bedtools makewindows -g ~/Analyze4C/temp/genome_sizes_temp.txt -w",winSize,"> ~/Analyze4C/temp/windows.bed"))
		}
		
		#applying a cutoff to the p-score data
		ans8 <- as.integer(readline(prompt=cat("\nenter the number of what type of cutoff method you would like to do for the p-score:\n\n1) all chromosomes with the same cutoff (in p-score)\n2) all chromosomes with the same cutoff percentage\n\n")))
		if(ans8 == 1)
		{
			CO1 <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff. Anything above the cutoff will be considered a positive RE site:\n\n")))
			CO_ps_type <- paste(CO_ps_type,";","fixed")
			CO_ps <- paste(CO_ps,";",CO1)
			
			ans4 <- as.integer(readline(prompt=cat("\nshould the cutoff be applied to:\n1) the whole genome\n2) all trans\n3) cis\n\n")))
			if(ans4 == 1) #whole genome
			{
				ps_afterCO <- ps[ps[,3]>=CO1,]
				CO_ps_chroms <- paste(CO_ps_chroms,";","all")
			}
			else if(ans4 == 2) #trans
			{
				ps_afterCO <- c()
				for(k in 1:numOFchroms) 
				{
					if(k!=cis)
					{
						ps_temp <- ps[ps[,1]==k & ps[,3]>=CO1,]
						ps_afterCO <- rbind(ps_afterCO,ps_temp)
					}
					else
					{
						ps_temp <- ps[ps[,1]==k,]
						ps_afterCO <- rbind(ps_afterCO,ps_temp)
					}
				}
				CO_ps_chroms <- paste(CO_ps_chroms,";","trans")
			}
			else if(ans4 == 3) #each chromosome separately
			{		
				ps_afterCO <- c()
				for(k in 1:numOFchroms) 
				{
					if(k==cis)
					{
						ps_temp <- ps[ps[,1]==k & ps[,3]>=CO1,]
						ps_afterCO <- rbind(ps_afterCO,ps_temp)
					}
					else
					{
						ps_temp <- ps[ps[,1]==k,]
						ps_afterCO <- rbind(ps_afterCO,ps_temp)
					}
				}
				CO_ps_chroms <- paste(CO_ps_chroms,";","cis")
			}				
		}
		else
		{
			ps_COper <- as.numeric(readline(prompt=cat("\nenter a cutoff percentage (between 0 and 1) that will be applied to the p-score\nonly the p-scores that are above or equal to the cutoff will be intersected with the FPKM data.\nif you don't want any cutoff, enter 0:\n")))			
			CO_ps_type <- paste(CO_ps_type,";","top percent per chromosome")
			CO_ps <- paste(CO_ps,";",ps_COper)

			ans4 <- as.integer(readline(prompt=cat("\nshould the cutoff percentage be of:\n1) the whole genome\n2) all trans\n3) each chromosome separately\n\n")))
			if(ans4 == 1) #whole genome
			{
				CO1 <- quantile(ps[,3],ps_COper)
				ps_afterCO <- ps[ps[,3]>=CO1,]
				CO_ps_chroms <- paste(CO_ps_chroms,";","all")
			}
			else if(ans4 == 2) #trans
			{
				ps_afterCO <- c()
				CO1 <- quantile(ps[ps[,1]!=cis,3],ps_COper)
				for(k in 1:numOFchroms) 
				{
					if(k!=cis)
					{
						ps_temp <- ps[ps[,1]==k & ps[,3]>=CO1,]
						ps_afterCO <- rbind(ps_afterCO,ps_temp)
					}
					else
					{
						ps_temp <- ps[ps[,1]==k,]
						ps_afterCO <- rbind(ps_afterCO,ps_temp)
					}
				}
				CO_ps_chroms <- paste(CO_ps_chroms,";","trans")
			}
			else if(ans4 == 3) #each chromosome separately
			{		
				ps_afterCO <- c()
				for(k in 1:numOFchroms) 
				{
					CO1 <- quantile(ps[ps[,1]==k,3],ps_COper)
					ps_temp <- ps[ps[,1]==k & ps[,3]>=CO1,]
					ps_afterCO <- rbind(ps_afterCO,ps_temp)
				}
				CO_ps_chroms <- paste(CO_ps_chroms,";","each chromosome separately")	
			}				
		}
					
		#adjusting the row numbers
		rownames(ps_afterCO) <- seq(length=nrow(ps_afterCO))
		
		#applying a cutoff to the FPKM data
		ans9 <- as.integer(readline(prompt=cat("\nenter the number of what type of cutoff method you would like to do for the FPKM:\n\n1) all chromosomes with the same cutoff (in FPKM)\n2) all chromosomes with the same cutoff percentage\n\n")))
		if(ans9 == 1)
		{
			CO2 <- as.numeric(readline(prompt=cat("\nplease enter an FPKM cutoff. Anything above the cutoff will be considered a positive RE site:\n\n")))
			CO_FPKM_type <- paste(CO_FPKM_type,";","fixed")
			CO_FPKM <- paste(CO_FPKM,";",CO2)
			
			ans5 <- as.integer(readline(prompt=cat("\nshould the cutoff be applied to:\n1) the whole genome\n2) all trans\n3) cis\n\n")))
			if(ans5 == 1) #whole genome
			{
				FPKM_afterCO <- FPKM[FPKM[,4]>=CO2,]
				CO_FPKM_chroms <- paste(CO_FPKM_chroms,";","all")
			}
			else if(ans5 == 2) #trans
			{
				FPKM_afterCO <- c()
				for(k in 1:numOFchroms) 
				{
					if(k!=cis)
					{
						FPKM_temp <- FPKM[FPKM[,1]==k & FPKM[,4]>=CO2,]
						FPKM_afterCO <- rbind(FPKM_afterCO,FPKM_temp)
					}
					else
					{
						FPKM_temp <- FPKM[FPKM[,1]==k,]
						FPKM_afterCO <- rbind(FPKM_afterCO,FPKM_temp)
					}
				}
				CO_FPKM_chroms <- paste(CO_FPKM_chroms,";","trans")
			}
			else if(ans5 == 3) #cis
			{		
				FPKM_afterCO <- c()
				for(k in 1:numOFchroms) 
				{
					if(k==cis)
					{
						FPKM_temp <- FPKM[FPKM[,1]==k & FPKM[,4]>=CO2,]
						FPKM_afterCO <- rbind(FPKM_afterCO,FPKM_temp)
					}
					else
					{
						FPKM_temp <- FPKM[FPKM[,1]==k,]
						FPKM_afterCO <- rbind(FPKM_afterCO,FPKM_temp)
					}
				}
				CO_FPKM_chroms <- paste(CO_FPKM_chroms,";","cis")
			}
		}
		else
		{	
			FPKM_COper <- as.numeric(readline(prompt=cat("\nenter a cutoff percentage (between 0 and 1) that will be applied to the FPKM\nonly the FPKMs that are above or equal to the cutoff will be intersected to the p-scores.\nif you don't want any cutoff, enter 0:\n")))
			CO_FPKM_type <- paste(CO_FPKM_type,";","top percent per chromosome")
			CO_FPKM <- paste(CO_FPKM,";",FPKM_COper)

			ans5 <- as.integer(readline(prompt=cat("\nshould the cutoff percentage be of:\n1) the whole genome\n2) all trans\n3) each chromosome separately\n\n")))
			if(ans5 == 1) #whole genome
			{
				CO2 <- quantile(FPKM[,4],FPKM_COper)
				FPKM_afterCO <- FPKM[FPKM[,4]>=CO2,]
				CO_FPKM_chroms <- paste(CO_FPKM_chroms,";","all")
			}
			else if(ans5 == 2) #trans
			{
				FPKM_afterCO <- c()
				CO2 <- quantile(FPKM[FPKM[,1]!=cis,4],FPKM_COper)
				for(k in 1:numOFchroms) 
				{
					if(k!=cis)
					{
						FPKM_temp <- FPKM[FPKM[,1]==k & FPKM[,4]>=CO2,]
						FPKM_afterCO <- rbind(FPKM_afterCO,FPKM_temp)
					}
					else
					{
						FPKM_temp <- FPKM[FPKM[,1]==k,]
						FPKM_afterCO <- rbind(FPKM_afterCO,FPKM_temp)
					}
				}
				CO_FPKM_chroms <- paste(CO_FPKM_chroms,";","trans")					
			}
			else if(ans5 == 3) #each chromosome separately
			{		
				FPKM_afterCO <- c()
				for(k in 1:numOFchroms) 
				{
					CO2 <- quantile(FPKM[FPKM[,1]==k,4],FPKM_COper)
					FPKM_temp <- FPKM[FPKM[,1]==k & FPKM[,4]>=CO2,]
					FPKM_afterCO <- rbind(FPKM_afterCO,FPKM_temp)
				}
				CO_FPKM_chroms <- paste(CO_FPKM_chroms,";","each chromosome separately")					
			}
		}
		
		#adjusting the row numbers			
		rownames(FPKM_afterCO) <- seq(length=nrow(FPKM_afterCO))
		
		#removing the translocated areas from 'ps_afterCO' and 'FPKM_afterCO' if they were switched
		if(translocated_flag == 1)
		{
			cat("\nremoving translocated areas:\nthe purpose of this is in order to be able and intersect the p-scores with FPKM bands\n")
			cat("the translocated areas indices mess up the order and will cause an error\n")
					
			#getting the line numbers that need to be removed
			for(d in 1:length(rearranged_rawData$Experiment))
			{
				sTemp <- unlist(strsplit(rearranged_rawData$Experiment[d],"_"))
				dtTemp <- unlist(strsplit(rearranged_rawData$Date_and_Time[d]," "))
				tTemp <- dtTemp[2]
				dtTemp[2] <- paste(unlist(strsplit(tTemp,":")),collapse="")
				if(all(sTemp %in% sp2) & !is.na(rearranged_rawData$Added_Sections_Lines[d]) & all(dtTemp %in% sp2))
				{
					lin_nums1 <- rearranged_rawData$Added_Sections_Lines[d]
					lin_nums2 <- as.integer(unlist(strsplit(lin_nums1,", ")))
					
					#getting the raw data file that the p-scores file was created from:
					if(!identical(grep("covRemoved",sp1),integer(0)))
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
							if(!identical(grep("covRemoved",sp1),integer(0)))
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
					secs_rmv <- "yes"
					break
				}			
			}

			#the following actions will work under these specific conditions: 
			#1) the rearranged raw data file contains all the original chromosomes. 
			#2) the line numbers in 'rearranged_rawData' file are divided so that each pair are of one specific chromosome (meaning that 'first' and 'second' are not of two different chromosomes)

			#removing the sections from the FPKM and ps data
			all_chroms <- seq(numOFchroms)
			pointer <- 1
			FPKM_fin <- c()
			ps_fin <- c()
			rawData_removed <- c() 
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
				
				#taking all the chromosomes that aren't getting the removal treatment, and adding them to the FPKM and ps data
				while(all_chroms[pointer]!=ch)
				{
					ps_t <- ps_afterCO[ps_afterCO[,1]==all_chroms[pointer],]
					ps_fin <- rbind(ps_fin,ps_t)
					FPKM_t <- FPKM_afterCO[FPKM_afterCO[,1]==all_chroms[pointer],]
					FPKM_fin <- rbind(FPKM_fin,FPKM_t)
					rawData_removed_temp <- raw_dat[raw_dat[,1]==all_chroms[pointer],] #creating the raw data removed data frame 
					rawData_removed <- rbind(rawData_removed,rawData_removed_temp) 
					pointer <- pointer + 1
				}
				
				ps_t <- c() #just in case no conditions are met, this way we won't get the previous ps_t (this shouldn't happen, it's just a precaution)
				FPKM_t <- c() #just in case no conditions are met, this way we won't get the previous FPKM_t (this shouldn't happen, it's just a precaution)
				#getting the indices to remove
				if(first==1)#if first is the 1st line in the file
				{
					if(second!=nrow(raw_dat)) #if second isn't the last of the file
					{			
						if(raw_dat[second,1]!=raw_dat[(second+1),1]) #if it covers all the chromosome (if second is the last in chromosome) (1)
						{
							ps_t <- c()
							FPKM_t <- c()
						}
						else if(raw_dat[second,1]==raw_dat[(second+1),1]) #if second isn't the last of the chromosome (2)
						{
							#remove anything below the index at second
							ps_t <- ps_afterCO[(ps_afterCO[,1]==ch)&(ps_afterCO[,2]>=raw_dat[(second+1),2]),]					
							FPKM_t <- FPKM_afterCO[(FPKM_afterCO[,1]==ch)&(FPKM_afterCO[,2]>=raw_dat[(second+1),2]),]
						}
					}	
					else #the line numbers cover all the data, so none of it is used (9)
					{
						ps_t <- c()	
						FPKM_t <- c()
					}				
				}
				else if(raw_dat[first,1]!=raw_dat[(first-1),1]) #if it is the first line in chromosome
				{
					if(second!=nrow(raw_dat)) #if second isn't last line
					{
						if(raw_dat[second,1]!=raw_dat[(second+1),1]) #if second is the last in chromosome (4)
						{
							ps_t <- c()
							FPKM_t <- c()
						}
						else if(raw_dat[second,1]==raw_dat[(second+1),1]) #if second isn't last in chromosome (3)
						{
							#remove anything below the index at second
							ps_t <- ps_afterCO[(ps_afterCO[,1]==ch)&(ps_afterCO[,2]>=raw_dat[(second+1),2]),]
							FPKM_t <- FPKM_afterCO[(FPKM_afterCO[,1]==ch)&(FPKM_afterCO[,2]>=raw_dat[(second+1),2]),]
						}
					}
					else #if second is last line of the file (7)
					{
						ps_t <- c()
						FPKM_t <- c()
					}	
				}
				else if(raw_dat[first,1]==raw_dat[(first-1),1])#if first is in the middle of a chromosome
				{
					if(second!=nrow(raw_dat)) #if second isn't last line
					{
						if(raw_dat[second,1]!=raw_dat[(second+1),1]) #if second is the last in chromosome (6)
						{
							#remove anything above the index at first
							ps_t <- ps_afterCO[(ps_afterCO[,1]==ch)&(ps_afterCO[,2]<=raw_dat[(first-1),2]),]
							FPKM_t <- FPKM_afterCO[(FPKM_afterCO[,1]==ch)&(FPKM_afterCO[,2]<=raw_dat[(first-1),2]),]
						}
						else if(raw_dat[second,1]==raw_dat[(second+1),1]) #if second isn't last in chromosome (5)
						{
							#remove anything above the index at first and below the index at second
							ps_t <- ps_afterCO[(ps_afterCO[,1]==ch)&((ps_afterCO[,2]<=raw_dat[(first-1),2])|(ps_afterCO[,2]>=raw_dat[(second+1),2])),]
							FPKM_t <- FPKM_afterCO[(FPKM_afterCO[,1]==ch)&((FPKM_afterCO[,2]<=raw_dat[(first-1),2])|(FPKM_afterCO[,2]>=raw_dat[(second+1),2])),]
						}
					}
					else #if second is last line of the file (8)
					{
						#remove anything above the index at first
						ps_t <- ps_afterCO[(ps_afterCO[,1]==ch)&(ps_afterCO[,2]<=raw_dat[(first-1),2]),]
						FPKM_t <- FPKM_afterCO[(FPKM_afterCO[,1]==ch)&(FPKM_afterCO[,2]<=raw_dat[(first-1),2]),]
					}	
				}
				ps_fin <- rbind(ps_fin,ps_t)	
				FPKM_fin <- rbind(FPKM_fin,FPKM_t)
				rawData_removed_temp <- raw_dat[-first:-second,] 
				rawData_removed <- rbind(rawData_removed,rawData_removed_temp[rawData_removed_temp[,1]==ch,]) 
				pointer <- pointer + 1
			}

			#adding the last chromosomes data (if there are any) that didn't get the section removal treatment
			while(pointer<=numOFchroms)
			{
				ps_t <- ps_afterCO[ps_afterCO[,1]==all_chroms[pointer],]
				ps_fin <- rbind(ps_fin,ps_t)	
				FPKM_t <- FPKM_afterCO[FPKM_afterCO[,1]==all_chroms[pointer],]
				FPKM_fin <- rbind(FPKM_fin,FPKM_t)
				rawData_removed_temp <- raw_dat[raw_dat[,1]==all_chroms[pointer],] 
				rawData_removed <- rbind(rawData_removed,rawData_removed_temp) 				
				pointer <- pointer + 1
			}

			ps_afterCO <- ps_fin
			rownames(ps_afterCO) <- seq(length=nrow(ps_afterCO))
			
			FPKM_afterCO <- FPKM_fin
			rownames(FPKM_afterCO) <- seq(length=nrow(FPKM_afterCO))
			
			rownames(rawData_removed) <- seq(length=nrow(rawData_removed)) 		
		}
		
		if(z == 1) #maybe i could remove this condition and have that every time you could enter which ever chrom you want, even different ones
		{
			choice2 <- as.integer(readline(prompt=cat("\nwhat chromosomes should the intersection be applied to? (enter the option number)\n1) all chromosomes\n2) trans\n3) cis\n4) specific chromosome\n\n")))
		}

		#getting the p-score and FPKM data by chromosome according to the request of the user, and also the chromosomes chosen
		if(choice2 == 1) #all chromosomes
		{
			ps_dat <- ps_afterCO
			FPKM_dat <- FPKM_afterCO
			chroms <- c(1:numOFchroms)
			int_chroms <- "all" #this will be entered later in 'expressionVScontacts_correlation_plots.txt'
		}
		else if(choice2 == 2) #trans
		{
			ps_dat <- ps_afterCO[ps_afterCO[,1] != cis,]
			FPKM_dat <- FPKM_afterCO[FPKM_afterCO[,1] != cis,]
			chroms <- Filter(function(x) x!=cis,1:numOFchroms) #this here excludes cis chromosome number in the sequence, another way would be - (1:numOFchroms)[-cis]
			int_chroms <- "trans" #this will be entered later in 'expressionVScontacts_correlation_plots.txt'
		}
		else if(choice2 == 3) #cis
		{
			ps_dat <- ps_afterCO[ps_afterCO[,1] == cis,]
			FPKM_dat <- FPKM_afterCO[FPKM_afterCO[,1] == cis,]			
			chroms <- cis
			int_chroms <- "cis" #this will be entered later in 'expressionVScontacts_correlation_plots.txt'
		}
		else if(choice2 == 4) #by chromosome
		{
			chr_chosen <- as.integer(readline(prompt=cat("\nchoose a chromosome you would like to use:\n\n")))
			ps_dat <- ps_afterCO[ps_afterCO[,1] == chr_chosen,]
			FPKM_dat <- FPKM_afterCO[FPKM_afterCO[,1] == chr_chosen,]			
			chroms <- chr_chosen
			int_chroms <- paste("chr",chr_chosen,sep="") #this will be entered later in 'expressionVScontacts_correlation_plots.txt'
		}

		#add an index to each p-score file
		ps_bed <- c()
		for(u in chroms)
		{
			if(nrow(ps_dat[ps_dat[,1]==u,]) != 0)
			{
				ps_chrom_temp <- ps_dat[ps_dat[,1]==u,]
				inds_temp <- ps_chrom_temp[,2]+1
				if(tail(inds_temp,1)>genome_sizes[u,2]) #if the addition of the extra index is larger than the size of the chromosome, we will remove that p-score line
				{
					inds_temp <- head(inds_temp,-1)
					ps_chrom_temp <- head(ps_chrom_temp,-1)
				}
				ps_chrom_bed <- cbind(ps_chrom_temp[,1],ps_chrom_temp[,2],inds_temp,ps_chrom_temp[,3])
				ps_bed <- rbind(ps_bed,ps_chrom_bed)
			}	
		}
		ps_bed <- data.frame(ps_bed)
		
		write.table(ps_bed,"~/Analyze4C/temp/ps_bed.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t") #creating a temp file of the p-score data chosen			
		write.table(FPKM_dat,"~/Analyze4C/temp/FPKM.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t") #creating a temp file of the FPKM data chosen			

		##!
		#getting the raw data that corresponds to the p-score file, in order to get the corresponding reads for each RE site
		#cat("\nthis is the current p-score file that is being tested:\n",ps.file,"\n")
		#repeat
		#{
		#	choice <- as.integer(readline(prompt=cat("\nenter the number of the folder which contains the raw data that the p-score file was created from (look at the name of the p-score file for clues):\n\n1) original\n2) rearranged\n3) coverage removed\n\n")))
		#	if(choice == 2)
		#	{
		#		ls_files <- system("ls ~/Analyze4C/rawData/rearranged",intern=TRUE)
		#	}
		#	else if(choice == 3)
		#	{
		#		ls_files <- system("ls ~/Analyze4C/rawData/coverage_removed",intern=TRUE)
		#	}
		#	else
		#	{
		#		ls_files <- system("ls ~/Analyze4C/rawData/original",intern=TRUE)
		#	}		
		#
		#	if(length(ls_files) != 0)
		#	{	
		#		file.name <- readline(prompt=cat("\nchoose from these raw data files (make sure the name is similar to the p-score files):\n",ls_files,"\n",sep="\n"))
		#		ind_files <- pmatch(file.name,ls_files,nomatch=0)
		#		if(ind_files == 0)
		#		{
		#			cat("no such file exists.\nplease try again.\n\n")
		#		}
		#		else
		#		{
		#			break
		#		}
		#	}
		#	else
		#	{
		#		cat("\nthis folder is empty. please choose a different folder from which you would like to use a file\n\n")
		#	}			
		#}
		#importing the data from the file
		#if(choice == 2)
		#{
		#	rawData <- read.table(paste("~/Analyze4C/rawData/rearranged/",file.name,sep=""))
		#}
		#else if(choice == 3)
		#{
		#	rawData <- read.table(paste("~/Analyze4C/rawData/coverage_removed/",file.name,sep=""))
		#}
		#else
		#{
		#	rawData <- read.table(paste("~/Analyze4C/rawData/original/",file.name,sep=""))
		#}
		##!
		
		#add an index to the rawData_removed file 
		rawData_bed <- c()
		for(u in chroms)
		{
			if(nrow(rawData_removed[rawData_removed[,1]==u,]) != 0) 
			{
				rawData_chrom_temp <- rawData_removed[rawData_removed[,1]==u,] 
				rawData_inds_temp <- rawData_chrom_temp[,2]+1
				if(tail(rawData_inds_temp,1)>genome_sizes[u,2]) #if the addition of the extra index is larger than the size of the chromosome, we will remove that p-score line
				{
					rawData_inds_temp <- head(rawData_inds_temp,-1)
					rawData_chrom_temp <- head(rawData_chrom_temp,-1)
				}
				rawData_chrom_bed <- cbind(rawData_chrom_temp[,1],rawData_chrom_temp[,2],rawData_inds_temp,rawData_chrom_temp[,3])
				rawData_bed <- rbind(rawData_bed,rawData_chrom_bed)
			}	
		}

		rawData_bed <- data.frame(rawData_bed)
		write.table(rawData_bed,"~/Analyze4C/temp/rawData_bed.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t") #creating a temp file of the rawData data chosen
		system("bedtools intersect -a ~/Analyze4C/temp/ps_bed.bed -b ~/Analyze4C/temp/rawData_bed.bed -loj > ~/Analyze4C/temp/reads.bed")
		reads <- read.table("~/Analyze4C/temp/reads.bed",header=FALSE,stringsAsFactors=FALSE,sep="\t")
		reads_bed <- as.data.frame(cbind(ps_bed[,1],ps_bed[,2],ps_bed[,3],reads[,8]))
		write.table(reads_bed,"~/Analyze4C/temp/reads_bed.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

		#intersecting	
		if(choice4 == 1) #if by individual
		{
			system("bedtools intersect -a ~/Analyze4C/temp/ps_bed.bed -b ~/Analyze4C/temp/FPKM.bed -wo > ~/Analyze4C/temp/intersected.bed")			
			intersected <- read.table("~/Analyze4C/temp/intersected.bed")	
			intersected_pairs <- as.data.frame(cbind(intersected[,4],intersected[,8]))
			
			system("bedtools intersect -a ~/Analyze4C/temp/reads_bed.bed -b ~/Analyze4C/temp/FPKM.bed -wo > ~/Analyze4C/temp/intersected_reads.bed")
			intersected_reads <- read.table("~/Analyze4C/temp/intersected_reads.bed")	
			intersected_reads_pairs <- as.data.frame(cbind(intersected_reads[,4],intersected_reads[,8]))			
		}	
		else #if by window
		{	
			system("bedtools map -a ~/Analyze4C/temp/windows.bed -b ~/Analyze4C/temp/ps_bed.bed -c 4 > ~/Analyze4C/temp/ps_windows.bed")
			ps_windows <- read.table("~/Analyze4C/temp/ps_windows.bed",header=FALSE,stringsAsFactors=FALSE,sep="\t")
			ps_windows[ps_windows[,4]==".",4] <- 0
			#write.table(ps_windows,"~/Analyze4C/temp/ps_windows.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

			system("bedtools map -a ~/Analyze4C/temp/windows.bed -b ~/Analyze4C/temp/FPKM.bed -c 4 > ~/Analyze4C/temp/FPKM_windows.bed")
			FPKM_windows <- read.table("~/Analyze4C/temp/FPKM_windows.bed",header=FALSE,stringsAsFactors=FALSE,sep="\t")
			FPKM_windows[FPKM_windows[,4]==".",4] <- 0
			#write.table(FPKM_windows,"~/Analyze4C/temp/FPKM_windows.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
			
			intersected_pairs <- as.data.frame(cbind(as.numeric(ps_windows[,4]),as.numeric(FPKM_windows[,4])))
			
			system("bedtools map -a ~/Analyze4C/temp/windows.bed -b ~/Analyze4C/temp/reads_bed.bed -c 4 > ~/Analyze4C/temp/reads_windows.bed")
			reads_windows <- read.table("~/Analyze4C/temp/reads_windows.bed",header=FALSE,stringsAsFactors=FALSE,sep="\t")
			reads_windows[reads_windows[,4]==".",4] <- 0
			#write.table(reads_windows,"~/Analyze4C/temp/reads_windows.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
			
			intersected_reads_pairs <- as.data.frame(cbind(as.numeric(reads_windows[,4]),as.numeric(FPKM_windows[,4])))			
		}
		
		#adding the name given by user to the data, and creating a column for it
		if(choice1 == 1)
		{
			ps_datName <- readline(prompt=cat("\nhow would you like this p-score data to be refered to in the plot? (no spaces)\n\n"))
			ps_datName_all <- c(ps_datName_all,ps_datName)
			intersected_pairs <- cbind(intersected_pairs,rep(ps_datName,nrow(intersected_pairs)))
			intersected_reads_pairs <- cbind(intersected_reads_pairs,rep(ps_datName,nrow(intersected_reads_pairs)))			
		}
		else
		{
			FPKM_datName <- readline(prompt=cat("\nhow would you like this FPKM data to be refered to in the plot? (no spaces)\n\n"))
			FPKM_datName_all <- c(FPKM_datName_all,FPKM_datName)
			intersected_pairs <- cbind(intersected_pairs,rep(FPKM_datName,nrow(intersected_pairs)))
			intersected_reads_pairs <- cbind(intersected_reads_pairs,rep(FPKM_datName,nrow(intersected_reads_pairs)))			
		}
		
		colnames(intersected_pairs) <- c("pScore","FPKM","names")
		#intersected_pairs <- as.data.frame(intersected_pairs)
		intersected_pairs_all <- rbind(intersected_pairs_all,intersected_pairs)

		colnames(intersected_reads_pairs) <- c("reads","FPKM","names")
		#intersected_reads_pairs <- as.data.frame(intersected_reads_pairs)
		intersected_reads_pairs_all <- rbind(intersected_reads_pairs_all,intersected_reads_pairs)
		
		#spearman
		cat("\nthe spearman correlation coefficient for the p-score vs. FPKM:\n")
		print(cor.test(intersected_pairs$pScore,intersected_pairs$FPKM,method="spearman",exact=FALSE))
		#print(cor.test(as.numeric(intersected_pairs$pScore),as.numeric(intersected_pairs$FPKM),method="spearman",exact=FALSE))
		#pearson
		cat("\nthe pearson correlation coefficient for the p-score vs. FPKM:\n")
		print(cor.test(intersected_pairs$pScore,intersected_pairs$FPKM))
		#print(cor.test(as.numeric(intersected_pairs$pScore),as.numeric(intersected_pairs$FPKM)))

		#spearman
		cat("\nthe spearman correlation coefficient for the reads vs. FPKM:\n")
		print(cor.test(intersected_reads_pairs_all$reads,intersected_reads_pairs_all$FPKM,method="spearman",exact=FALSE))
		#print(cor.test(as.numeric(intersected_reads_pairs_all$reads),as.numeric(intersected_reads_pairs_all$FPKM),method="spearman",exact=FALSE))
		#pearson
		cat("\nthe pearson correlation coefficient for the reads vs. FPKM:\n")
		print(cor.test(intersected_reads_pairs_all$reads,intersected_reads_pairs_all$FPKM))
		#print(cor.test(as.numeric(intersected_reads_pairs_all$reads),as.numeric(intersected_reads_pairs_all$FPKM)))

		if(choice1 == 1)
		{
			#asking if to choose another p-score file to correlate
			ans6 <- readline(prompt=cat("\nwould you like to choose another p-score file to intersect with the FPKM data?\ny/n\n\n"))
			if(ans6 == "y")
			{
				z <- z + 1				
				system("rm ~/Analyze4C/temp/ps_bed.bed")
				system("rm ~/Analyze4C/temp/reads.bed")
				system("rm ~/Analyze4C/temp/reads_bed.bed")
				system("rm ~/Analyze4C/temp/rawData_bed.bed")
				if(choice4 == 1)
				{
					system("rm ~/Analyze4C/temp/intersected.bed")
					system("rm ~/Analyze4C/temp/intersected_reads.bed")
				}
				else
				{
					winds <- c(winds," ; ",winSize)
					system("rm ~/Analyze4C/temp/ps_windows.bed")
					system("rm ~/Analyze4C/temp/reads_windows.bed")
					system("rm ~/Analyze4C/temp/FPKM_windows.bed")
					choice5 <- readline(prompt=cat("\nwould you like to choose a different window size?\ny/n\n\n"))
					if(choice5 == "y")
					{
						system("rm ~/Analyze4C/temp/windows.bed")
						wind_flag <- 0
					}
					else
					{
						wind_flag <- 1
					}	
				}	
			}
			else
			{
				files_flag <- 1
				num_of_pScore_files <- z
				#erasing all the files from temp
				system("rm ~/Analyze4C/temp/*")				
			}
		}
		else
		{
			#asking if to choose another FPKM file to correlate
			ans6 <- readline(prompt=cat("\nwould you like to choose another FPKM file to intersect with the FPKM data?\ny/n\n\n"))
			if(ans6 == "y")
			{
				z <- z + 1
				system("rm ~/Analyze4C/temp/FPKM.bed")
				system("rm ~/Analyze4C/temp/reads.bed")
				system("rm ~/Analyze4C/temp/reads_bed.bed")				
				system("rm ~/Analyze4C/temp/rawData_bed.bed")				
				if(choice4 == 1)
				{
					system("rm ~/Analyze4C/temp/intersected.bed")
					system("rm ~/Analyze4C/temp/intersected_reads.bed")
				}
				else
				{
					winds <- c(winds," ; ",winSize)
					system("rm ~/Analyze4C/temp/FPKM_windows.bed")
					system("rm ~/Analyze4C/temp/ps_windows.bed")
					system("rm ~/Analyze4C/temp/reads_windows.bed")
					choice5 <- readline(prompt=cat("\nwould you like to choose a different window size?\ny/n\n\n"))
					if(choice5 == "y")
					{
						system("rm ~/Analyze4C/temp/windows.bed")
						wind_flag <- 0
					}
					else
					{
						wind_flag <- 1
					}	
				}
			}
			else
			{
				files_flag <- 1
				num_of_FPKM_files <- z
				#erasing all the files from temp
				system("rm ~/Analyze4C/temp/*")
			}		
		}
	}

	#getting the time and date
	DandT1 <- toString(Sys.time())
	DandT2 <- gsub(" ","_",DandT1)
	DandT2 <- gsub(":","",DandT2)
	
	#intersected_pairs_all <- as.data.frame(intersected_pairs_all)
	cat("\ncreating the p-score VS. FPKM plot\nplease wait...\n\n")
	#print(ggplot2::ggplot(intersected_pairs_all, ggplot2::aes(x=pScore,y=FPKM,color=names)) + ggplot2::scale_x_log10() +  ggplot2::scale_y_log10() + ggplot2::geom_point() + ggplot2::geom_smooth(method=lm,se=FALSE,fullrange=TRUE))
	#print(ggplot2::ggplot(intersected_pairs_all, ggplot2::aes(x=pScore,y=FPKM,color=names)) + ggplot2::geom_point() + ggplot2::geom_smooth(method=lm,se=FALSE,fullrange=TRUE) + ggplot2::scale_y_continuous(limits = c(0,5000)) + ggplot2::ggtitle("contact VS. expression - correlation") + ggplot2::theme(plot.title = ggplot2::element_text(size=15,hjust = 0.5)))
	print(ggplot2::ggplot(intersected_pairs_all, ggplot2::aes(x=pScore,y=FPKM,color=names)) + ggplot2::geom_point() + ggplot2::geom_smooth(method=lm,se=FALSE,fullrange=TRUE) + ggplot2::ggtitle("contact VS. expression - correlation") + ggplot2::theme(plot.title = ggplot2::element_text(size=15,hjust = 0.5)))

	#saving the plot if the user wants
	ans7 <- readline(prompt=cat("\nwould you like to save this plot?\ny/n\n\n"))
	if(ans7 == "y")
	{
		#save the plot
		nm <- paste("expressionVScontacts_correlation_plots_pScores_vs_FPKM_",DandT2,".png",sep="")
		wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
		ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
		cat("\nsaving the plot\nplease wait...\n\n")
		ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
		save_flag <- 1	
	}

	#intersected_reads_pairs_all <- as.data.frame(intersected_reads_pairs_all)	
	cat("\ncreating the reads VS. FPKM plot\nplease wait...\n\n")	
	#print(ggplot2::ggplot(intersected_reads_pairs_all, ggplot2::aes(x=reads,y=FPKM,color=names)) + ggplot2::scale_x_log10() +  ggplot2::scale_y_log10() + ggplot2::geom_point() + ggplot2::geom_smooth(method=lm,se=FALSE,fullrange=TRUE))
	#print(ggplot2::ggplot(intersected_reads_pairs_all, ggplot2::aes(x=reads,y=FPKM,color=names)) + ggplot2::geom_point() + ggplot2::geom_smooth(method=lm,se=FALSE,fullrange=TRUE) + ggplot2::scale_y_continuous(limits = c(0,5000)) + ggplot2::ggtitle("contact VS. expression - reads correlation") + ggplot2::theme(plot.title = ggplot2::element_text(size=15,hjust = 0.5)))
	print(ggplot2::ggplot(intersected_reads_pairs_all, ggplot2::aes(x=reads,y=FPKM,color=names)) + ggplot2::geom_point() + ggplot2::geom_smooth(method=lm,se=FALSE,fullrange=TRUE) + ggplot2::ggtitle("contact VS. expression - reads correlation") + ggplot2::theme(plot.title = ggplot2::element_text(size=15,hjust = 0.5)))

	#saving the plot if the user wants
	ans10 <- readline(prompt=cat("\nwould you like to save this plot?\ny/n\n\n"))
	if(ans10 == "y")
	{
		#save the plot
		nm <- paste("expressionVScontacts_correlation_plots_reads_vs_FPKM_",DandT2,".png",sep="")
		wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
		ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
		cat("\nsaving the plot\nplease wait...\n\n")
		ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
		save_flag <- 1	
	}
	
	if(save_flag == 1)
	{
		#replacing the " ; " at the start of the lines with nothing
		CO_ps_type <- gsub("^ ; ","",CO_ps_type)
		CO_ps <- gsub("^ ; ","",CO_ps)
		CO_ps_chroms <- gsub("^ ; ","",CO_ps_chroms)
		CO_FPKM_type <- gsub("^ ; ","",CO_FPKM_type)
		CO_FPKM <- gsub("^ ; ","",CO_FPKM)
		CO_FPKM_chroms <- gsub("^ ; ","",CO_FPKM_chroms)
		
		if(choice4 == 2)
		{
			winds <- gsub("^ ; ","",winds)
		}
		
		#saving the parameters and details to expressionVScontacts_correlation_plots.txt
		expressionVScontacts_correlation_plots[nrow(expressionVScontacts_correlation_plots)+1,] <- c(psORreads,num_of_pScore_files,ps_Experiments,num_of_FPKM_files,FPKM_Experiments,CO_ps_type,CO_ps,CO_ps_chroms,CO_FPKM_type,CO_FPKM,CO_FPKM_chroms,int_chroms,secs_rmv,by_winds,winds,DandT1)
		#sorting the list of experiments by bait alphabetically (and sorting the row indices)
		expressionVScontacts_correlation_plots <- expressionVScontacts_correlation_plots[order(expressionVScontacts_correlation_plots$Date_and_Time),]
		rownames(expressionVScontacts_correlation_plots) <- seq(length=nrow(expressionVScontacts_correlation_plots))
		#adding the new data to the file (by erasing the old one and creating a new one)
		system("rm expressionVScontacts_correlation_plots.txt")
		write.table(expressionVScontacts_correlation_plots,"expressionVScontacts_correlation_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)	
	}
}
