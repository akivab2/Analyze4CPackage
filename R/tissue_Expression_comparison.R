#' @export

#this function compares one contact bands file with 2 different expression (FPKM) files
#it will find the difference between the comparison of each FPKM data with the contact bands data

#need to add an option of saving the intersection files in order to view in the browser - probably don't need this, since i could just take the contact bands file and the FPKM file and apply the cutoff in the browser

tissue_Expression_comparison <- function(Experiments_4C,expressionVScontacts_plots,rearranged_rawData)
{
	cat("\n**************************************************************************\n\n")
	cat("\n                        Tissue expression comparison\n\n")
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
	
	#aking user if to use already existing contact band file or create a new one	
	ans1 <- as.numeric(readline(prompt=cat("\nwould you like to:\n1) choose an already existing contact bands file\n2) create a new contact bands data file\n\n")))
	if(ans1 == 1) #taking an already existing contact bands file
	{
		ls_conts <- system("ls ~/Analyze4C/contact_bands/",intern=TRUE)
		repeat
		{
			if(length(ls_conts) != 0)
			{	
				conts_filename <- readline(prompt=cat("\nchoose the contact bands file that you would like to use:\n",ls_conts,"",sep="\n"))
				ind_files1 <- pmatch(conts_filename,ls_conts,nomatch=0)
				if(ind_files1 == 0)
				{
					cat("no such file exists.\nplease try again.\n\n")
				}
				else
				{
					conts <- read.table(paste("~/Analyze4C/contact_bands/",conts_filename,sep=""))
					#conts <- remove_translocated(conts,"y")
					
					
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
		conts_filename <- makeContactBands.mainMenu(Experiments_4C,rearranged_rawData)
		conts <- read.table(paste("~/Analyze4C/contact_bands/",conts_filename,sep=""))
	}

	#finding the specific contact bands information in Experiments_4C
	sp1 <- unlist(strsplit(conts_filename,"_"))
	out <- findIn.Experiments_4C(sp1[1],sp1[2],sp1[3],Experiments_4C)

	#getting the cis chromosome number
	cis <- as.numeric(out[2])
		
	#get two different rnaseq FPKM files:
	cat("\n\nchoosing the FPKM files:\n\nyou must choose 2 different FPKM data files\n")
	ls_FPKM <- system("ls ~/Analyze4C/RNAseq/FPKM | grep '.bed$'",intern=TRUE)
	
	if(length(ls_FPKM) >= 2)
	{	
		repeat
		{
			FPKM_filename1 <- readline(prompt=cat("\nchoose the first FPKM file that you would like to use:\n",ls_FPKM,"",sep="\n"))
			ind_files2 <- pmatch(FPKM_filename1,ls_FPKM,nomatch=0)
			if(ind_files2 == 0)
			{
				cat("no such file exists.\nplease try again.\n\n")
			}
			else
			{
				FPKM1 <- read.table(paste("~/Analyze4C/RNAseq/FPKM/",FPKM_filename1,sep=""))
				break
			}
		}
		
		repeat
		{
			FPKM_filename2 <- readline(prompt=cat("\nchoose the second FPKM file that you would like to use:\n",ls_FPKM[ls_FPKM != FPKM_filename1],"",sep="\n"))
			ind_files3 <- pmatch(FPKM_filename2,ls_FPKM[ls_FPKM != FPKM_filename1],nomatch=0)
			if(ind_files3 == 0)
			{
				cat("no such file exists.\nplease try again.\n\n")
			}
			else
			{
				FPKM2 <- read.table(paste("~/Analyze4C/RNAseq/FPKM/",FPKM_filename2,sep=""))
				break
			}
		}
		
		#asking if to not use all the FPKMs that are 0, this is important because when getting the quantiles it considers the 0 examples as well
		#and since FPKMs of 0 could be considered as if there is no expression at all in those areas, then we might be able to remove them
		#the user decided
		remZeroFPKM <- readline(prompt=cat("\nshould FPKM's of 0 be removed?\ny/n\n\n"))
		if(remZeroFPKM == "y") #removing all FPKMs that are 0 (since they are basically rna-seq reads that don't exist)
		{
			FPKM1 <- FPKM1[FPKM1[,4]>0,]
			FPKM2 <- FPKM2[FPKM2[,4]>0,]
		}		
	}
	else
	{
		cat("\nat least two different FPKM bed files are required to be in the folder\nit seems threre are less than two currently\n")
		cat("\nif you would like to execute this function, first download more RNA-seq datasets\n\n")
		return()
	}
	
	#applying a cutoff to the FPKM data
	CO_type <- as.numeric(readline(prompt=cat("\napplying a cutoff to the FPKM data:\n\nhow would you like to apply the cutoff?\n1) by percentage\n2) by a set FPKM value\n\n")))
	if(CO_type == 1)
	{
		FPKM_COper <- as.numeric(readline(prompt=cat("\nenter a cutoff percentage (between 0 and 1) that will be applied to the FPKM\nonly the FPKMs that are above or equal to the cutoff will be intersected to the contact bands.\nif you don't want any cutoff, enter 0:\n")))
		CO_type <- "by percentage"
		FPKM_cutoff <- FPKM_COper
		ans4 <- as.integer(readline(prompt=cat("\nshould the cutoff percentage be of and applied to:\n1) the whole genome\n2) all trans\n3) each chromosome separately\n\n")))
		if(ans4 == 1) #whole genome
		{
			CO1 <- quantile(FPKM1[,4],FPKM_COper)
			FPKM1_afterCO <- FPKM1[FPKM1[,4]>=CO1,]
			CO2 <- quantile(FPKM2[,4],FPKM_COper)
			FPKM2_afterCO <- FPKM2[FPKM2[,4]>=CO2,]
		}
		else if(ans4 == 2) #trans
		{
			FPKM1_afterCO <- c()
			FPKM2_afterCO <- c()
			CO1 <- quantile(FPKM1[FPKM1[,1]!=cis,4],FPKM_COper)
			CO2 <- quantile(FPKM2[FPKM2[,1]!=cis,4],FPKM_COper)
			for(k in 1:numOFchroms) 
			{
				if(k!=cis)
				{
					FPKM1_temp <- FPKM1[FPKM1[,1]==k & FPKM1[,4]>=CO1,]
					FPKM1_afterCO <- rbind(FPKM1_afterCO,FPKM1_temp)
					
					FPKM2_temp <- FPKM2[FPKM2[,1]==k & FPKM2[,4]>=CO2,]
					FPKM2_afterCO <- rbind(FPKM2_afterCO,FPKM2_temp)
				}
				else
				{
					FPKM1_temp <- FPKM1[FPKM1[,1]==k,]
					FPKM1_afterCO <- rbind(FPKM1_afterCO,FPKM1_temp)
					
					FPKM2_temp <- FPKM2[FPKM2[,1]==k,]
					FPKM2_afterCO <- rbind(FPKM2_afterCO,FPKM2_temp)
				}
			}	
		}
		else if(ans4 == 3) #specific chromosomes
		{		
			FPKM1_afterCO <- c()
			FPKM2_afterCO <- c()
			for(k in 1:numOFchroms)
			{
				CO1 <- quantile(FPKM1[FPKM1[,1]==k,4],FPKM_COper)
				FPKM1_temp <- FPKM1[FPKM1[,1]==k & FPKM1[,4]>=CO1,]
				FPKM1_afterCO <- rbind(FPKM1_afterCO,FPKM1_temp)
				
				CO2 <- quantile(FPKM2[FPKM2[,1]==k,4],FPKM_COper)
				FPKM2_temp <- FPKM2[FPKM2[,1]==k & FPKM2[,4]>=CO2,]
				FPKM2_afterCO <- rbind(FPKM2_afterCO,FPKM2_temp)
			}		
		}	
	}
	else
	{
		CO1 <- as.numeric(readline(prompt=cat("\nenter a cutoff (in FPKM) that will be applied to",FPKM_filename1,"\nonly the FPKMs that are above or equal to the cutoff will be intersected to the contact bands.\nif you don't want any cutoff, enter 0:\n")))
		CO2 <- as.numeric(readline(prompt=cat("\nenter a cutoff (in FPKM) that will be applied to",FPKM_filename2,"\nonly the FPKMs that are above or equal to the cutoff will be intersected to the contact bands.\nif you don't want any cutoff, enter 0:\n")))
		CO_type <- "by value"
		FPKM_cutoff <- c(CO1,CO2)
		ans4 <- as.integer(readline(prompt=cat("\nshould the cutoffs be applied to:\n1) the whole genome\n2) all trans\n3) each chromosome separately\n\n")))
		if(ans4 == 1) #whole genome
		{
			FPKM1_afterCO <- FPKM1[FPKM1[,4]>=CO1,]
			FPKM2_afterCO <- FPKM2[FPKM2[,4]>=CO2,]
		}
		else if(ans4 == 2) #trans
		{
			FPKM1_afterCO <- c()
			FPKM2_afterCO <- c()
			for(k in 1:numOFchroms) 
			{
				if(k!=cis)
				{
					FPKM1_temp <- FPKM1[FPKM1[,1]==k & FPKM1[,4]>=CO1,]
					FPKM1_afterCO <- rbind(FPKM1_afterCO,FPKM1_temp)
					
					FPKM2_temp <- FPKM2[FPKM2[,1]==k & FPKM2[,4]>=CO2,]
					FPKM2_afterCO <- rbind(FPKM2_afterCO,FPKM2_temp)
				}
				else
				{
					FPKM1_temp <- FPKM1[FPKM1[,1]==k,]
					FPKM1_afterCO <- rbind(FPKM1_afterCO,FPKM1_temp)
					
					FPKM2_temp <- FPKM2[FPKM2[,1]==k,]
					FPKM2_afterCO <- rbind(FPKM2_afterCO,FPKM2_temp)
				}
			}	
		}
		else if(ans4 == 3) #specific chromosomes
		{		
			FPKM1_afterCO <- c()
			FPKM2_afterCO <- c()
			for(k in 1:numOFchroms)
			{
				FPKM1_temp <- FPKM1[FPKM1[,1]==k & FPKM1[,4]>=CO1,]
				FPKM1_afterCO <- rbind(FPKM1_afterCO,FPKM1_temp)
				
				FPKM2_temp <- FPKM2[FPKM2[,1]==k & FPKM2[,4]>=CO2,]
				FPKM2_afterCO <- rbind(FPKM2_afterCO,FPKM2_temp)
			}		
		}	
	}
	
	#removing the translocated areas from the FPKM data, the same areas that were removed from the contact bands files
	if(("rearranged" %in% sp1) & !("removed" %in% sp1)) #if the conts filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
	{	
		#getting the line numbers that need to be removed
		for(d in 1:length(rearranged_rawData$Experiment))
		{
			sTemp <- unlist(strsplit(rearranged_rawData$Experiment[d],"_"))
			dtTemp <- unlist(strsplit(rearranged_rawData$Date_and_Time[d]," "))
			tTemp <- dtTemp[2]
			dtTemp[2] <- paste(unlist(strsplit(tTemp,":")),collapse="")
			if(all(sTemp %in% sp1) & !is.na(rearranged_rawData$Added_Sections_Lines[d]) & all(dtTemp %in% sp1))
			{
				lin_nums1 <- rearranged_rawData$Added_Sections_Lines[d]
				lin_nums2 <- as.integer(unlist(strsplit(lin_nums1,", ")))

				#getting the raw data file that the contact bands file was created from:
				if(!identical(grep("covRemoved",sp1[e]),integer(0)))
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
						if(!identical(grep("covRemoved",sp1[e]),integer(0)))
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
		
		#removing the sections from the FPKM data
		all_chroms <- seq(numOFchroms)
		pointer <- 1
		FPKM1_fin <- c()
		FPKM2_fin <- c()
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
			
			#taking all the chromosomes that aren't getting the removal treatment, and adding them to the FPKM data
			while(all_chroms[pointer]!=ch)
			{
				FPKM1_t <- FPKM1_afterCO[FPKM1_afterCO[,1]==all_chroms[pointer],]
				FPKM2_t <- FPKM2_afterCO[FPKM2_afterCO[,1]==all_chroms[pointer],]	
				FPKM1_fin <- rbind(FPKM1_fin,FPKM1_t)
				FPKM2_fin <- rbind(FPKM2_fin,FPKM2_t)
				pointer <- pointer + 1
			}
			
			#just in case no conditions are met, this way we won't get the previous FPKM1_t and FPKM2_t(this shouldn't happen, it's just a precaution)
			FPKM1_t <- c()
			FPKM2_t <- c()
			#getting the indices to remove
			if(first==1)#if first is the 1st line in the file
			{
				if(second!=nrow(raw_dat)) #if second isn't the last of the file
				{			
					if(raw_dat[second,1]!=raw_dat[(second+1),1]) #if it covers all the chromosome (if second is the last in chromosome) (1)
					{
						FPKM1_t <- c()
						FPKM2_t <- c()
					}
					else if(raw_dat[second,1]==raw_dat[(second+1),1]) #if second isn't the last of the chromosome (2)
					{
						#remove anything below the index at second
						FPKM1_t <- FPKM1_afterCO[(FPKM1_afterCO[,1]==ch)&(FPKM1_afterCO[,2]>=raw_dat[(second+1),2]),]
						FPKM2_t <- FPKM2_afterCO[(FPKM2_afterCO[,1]==ch)&(FPKM2_afterCO[,2]>=raw_dat[(second+1),2]),]
					}
				}	
				else #the line numbers cover all the data, so none of it is used (9)
				{
					FPKM1_t <- c()
					FPKM2_t <- c()
				}				
			}
			else if(raw_dat[first,1]!=raw_dat[(first-1),1]) #if it is the first line in chromosome
			{
				if(second!=nrow(raw_dat)) #if second isn't last line
				{
					if(raw_dat[second,1]!=raw_dat[(second+1),1]) #if second is the last in chromosome (4)
					{
						FPKM1_t <- c()
						FPKM2_t <- c()
					}
					else if(raw_dat[second,1]==raw_dat[(second+1),1]) #if second isn't last in chromosome (3)
					{
						#remove anything below the index at second
						FPKM1_t <- FPKM1_afterCO[(FPKM1_afterCO[,1]==ch)&(FPKM1_afterCO[,2]>=raw_dat[(second+1),2]),]
						FPKM2_t <- FPKM2_afterCO[(FPKM2_afterCO[,1]==ch)&(FPKM2_afterCO[,2]>=raw_dat[(second+1),2]),]
					}
				}
				else #if second is last line of the file (7)
				{
					FPKM1_t <- c()
					FPKM2_t <- c()
				}	
			}
			else if(raw_dat[first,1]==raw_dat[(first-1),1])#if first is in the middle of a chromosome
			{
				if(second!=nrow(raw_dat)) #if second isn't last line
				{
					if(raw_dat[second,1]!=raw_dat[(second+1),1]) #if second is the last in chromosome (6)
					{
						#remove anything above the index at first
						FPKM1_t <- FPKM1_afterCO[(FPKM1_afterCO[,1]==ch)&(FPKM1_afterCO[,2]<=raw_dat[(first-1),2]),]
						FPKM2_t <- FPKM2_afterCO[(FPKM2_afterCO[,1]==ch)&(FPKM2_afterCO[,2]<=raw_dat[(first-1),2]),]
					}
					else if(raw_dat[second,1]==raw_dat[(second+1),1]) #if second isn't last in chromosome (5)
					{
						#remove anything above the index at first and below the index at second
						FPKM1_t <- FPKM1_afterCO[(FPKM1_afterCO[,1]==ch)&((FPKM1_afterCO[,2]<=raw_dat[(first-1),2])|(FPKM1_afterCO[,2]>=raw_dat[(second+1),2])),]
						FPKM2_t <- FPKM2_afterCO[(FPKM2_afterCO[,1]==ch)&((FPKM2_afterCO[,2]<=raw_dat[(first-1),2])|(FPKM2_afterCO[,2]>=raw_dat[(second+1),2])),]
					}
				}
				else #if second is last line of the file (8)
				{
					#remove anything above the index at first
					FPKM1_t <- FPKM1_afterCO[(FPKM1_afterCO[,1]==ch)&(FPKM1_afterCO[,2]<=raw_dat[(first-1),2]),]
					FPKM2_t <- FPKM2_afterCO[(FPKM2_afterCO[,1]==ch)&(FPKM2_afterCO[,2]<=raw_dat[(first-1),2]),]
				}	
			}			
			FPKM1_fin <- rbind(FPKM1_fin,FPKM1_t)
			FPKM2_fin <- rbind(FPKM2_fin,FPKM2_t)
			pointer <- pointer + 1
		}
		
		#adding the last chromosomes data (if there are any) that didn't get the section removal treatment
		while(pointer<=numOFchroms)
		{
			FPKM1_t <- FPKM1_afterCO[FPKM1_afterCO[,1]==all_chroms[pointer],]
			FPKM2_t <- FPKM2_afterCO[FPKM2_afterCO[,1]==all_chroms[pointer],]	
			FPKM1_fin <- rbind(FPKM1_fin,FPKM1_t)
			FPKM2_fin <- rbind(FPKM2_fin,FPKM2_t)
			pointer <- pointer + 1
		}
		
		FPKM1_afterCO <- FPKM1_fin
		FPKM2_afterCO <- FPKM2_fin
		rownames(FPKM1_afterCO) <- seq(length=nrow(FPKM1_afterCO))
		rownames(FPKM2_afterCO) <- seq(length=nrow(FPKM2_afterCO))		
	}
		
	write.table(FPKM1_afterCO,"~/Analyze4C/temp/FPKM1_afterCO.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
	write.table(FPKM2_afterCO,"~/Analyze4C/temp/FPKM2_afterCO.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
	
	#wherever there is a contact - intersect with both FPKM files
	cat("\nintersecting the contact bands file and the expression data from the two files that you chose\nplease wait...\n\n")
	system(paste("bedtools intersect -a ~/Analyze4C/temp/FPKM1_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename," -wo > ~/Analyze4C/temp/expressionComparison1.bed",sep=""))
	system(paste("bedtools intersect -a ~/Analyze4C/temp/FPKM2_afterCO.bed -b ~/Analyze4C/contact_bands/",conts_filename," -wo > ~/Analyze4C/temp/expressionComparison2.bed",sep=""))
	expressionComparison1 <- read.table("~/Analyze4C/temp/expressionComparison1.bed")	
	expressionComparison2 <- read.table("~/Analyze4C/temp/expressionComparison2.bed")

	#getting the date and time in order to distinguish between file names of plots
	DandT1 <- toString(Sys.time())
	DandT2 <- gsub(" ","_",DandT1)
	DandT2 <- gsub(":","",DandT2)
	
	#compare the numbers of both - ask user if to compare all, trans, cis , or individual chromosomes
	repeat
	{
		cat("\ncomparing the intersections of the expressions:\n")
		ans3 <- as.numeric(readline(prompt=cat("\nchoose a chromosome that you would like to compare the intersections for:\n\n1) All\n2) trans\n3) specific chromosome (will include the cis chromosome)\n4) exit\n\n")))
		if(ans3 == 1) #all
		{
			#1) sum all FPKM and compare
			#2) get a distribution of the FPKMs
			#3) compare specific contacts and the results of FPKM, this is only when both FPKM files intersect at the same contact band
			#4) get the average of FPKM for each file and compare them
			
			#the mcnemar test
			ans.mcnemar <- readline("\nwould you like to perform a mcnemar test?\ny/n\n\n")
			if(ans.mcnemar == "y")
			{
				mcnemar_tester(paste("~/Analyze4C/contact_bands/",conts_filename,sep=""),"~/Analyze4C/temp/FPKM1_afterCO.bed","~/Analyze4C/temp/FPKM2_afterCO.bed")
			}

			ans.prop <- readline(prompt=cat("\n\nwould you like to perform a prop test on the intersected bp and the intersected FPKM?\ny/n\n\n"))
			if(ans.prop == "y")
			{	
				#BPs:

				#getting the sum of all the bps of the intersections
				intersectionBPsum1 <- sum(expressionComparison1[,8])
				intersectionBPsum2 <- sum(expressionComparison2[,8])

				#getting the sum of all the bps of the contact bands
				all_contsBPsum <- sum(conts[,3]-conts[,2])
				
				#getting the sum of all the bps that are above the FPKM cutoff
				all_fpkmBPsum1 <- sum(FPKM1_afterCO[,3]-FPKM1_afterCO[,2])
				all_fpkmBPsum2 <- sum(FPKM2_afterCO[,3]-FPKM2_afterCO[,2])
				
				#FPKMs:

				#getting the sum of all the FPKMs of the intersections
				expressionComparison1_all <- expressionComparison1[row.names(unique(expressionComparison1[,c(colnames(expressionComparison1)[1],colnames(expressionComparison1)[2])])),]
				intersectionFPKMsum1 <- sum(expressionComparison1_all[,4])
				expressionComparison2_all <- expressionComparison2[row.names(unique(expressionComparison2[,c(colnames(expressionComparison2)[1],colnames(expressionComparison2)[2])])),]
				intersectionFPKMsum2 <- sum(expressionComparison2_all[,4])
				
				#getting the sum of all the FPKMs above the FPKM cutoff
				all_fpkmSum1 <- sum(FPKM1_afterCO[,4])
				all_fpkmSum2 <- sum(FPKM2_afterCO[,4])
				
				cat("\nthe contact bands file chosen:\n",conts_filename,"\n")
				cat("\nthe number of bps of contact bands on the whole genome:",all_contsBPsum,"\n")
				cat("\nthe number of bps of FPKMs above the cutoff on the whole genome of the expression file",FPKM_filename1,"is:",all_fpkmBPsum1,"\n")
				cat("\nthe number of bps of FPKMs above the cutoff on the whole genome of the expression file",FPKM_filename2,"is:",all_fpkmBPsum2,"\n")
				cat("\nthe sum of FPKMs above the cutoff on the whole genome of the expression file",FPKM_filename1,"is:",all_fpkmSum1,"\n")
				cat("\nthe sum of FPKMs above the cutoff on the whole genome of the expression file",FPKM_filename2,"is:",all_fpkmSum2,"\n")			
				cat("\nexpression file - ",FPKM_filename1,":\nthe number of bps intersected are ",intersectionBPsum1," which are ",(intersectionBPsum1/all_contsBPsum)*100,"% of all the bps of the contact bands in the whole genome\n\n",sep="")
				cat("\nexpression file - ",FPKM_filename2,":\nthe number of bps intersected are ",intersectionBPsum2," which are ",(intersectionBPsum2/all_contsBPsum)*100,"% of all the bps of the contact bands in the whole genome\n\n",sep="")
				cat("\nexpression file - ",FPKM_filename1,":\nthe number of bps intersected are ",intersectionBPsum1," which are ",(intersectionBPsum1/all_fpkmBPsum1)*100,"% of all the bps of FPKMs above the cutoff\n\n",sep="")
				cat("\nexpression file - ",FPKM_filename2,":\nthe number of bps intersected are ",intersectionBPsum2," which are ",(intersectionBPsum2/all_fpkmBPsum2)*100,"% of all the bps of FPKMs above the cutoff\n\n",sep="")
				cat("\nexpression file - ",FPKM_filename1,":\nthe sum of FPKMs intersected are ",intersectionFPKMsum1," which are ",(intersectionFPKMsum1/all_fpkmSum1)*100,"% of all the FPKMs above the cutoff\n\n",sep="")
				cat("\nexpression file - ",FPKM_filename2,":\nthe sum of FPKMs intersected are ",intersectionFPKMsum2," which are ",(intersectionFPKMsum2/all_fpkmSum2)*100,"% of all the FPKMs above the cutoff\n\n",sep="")			

				#creating plot		
				barName1 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",FPKM_filename1,"on the plot:\n"))
				barName2 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",FPKM_filename2,"on the plot:\n"))
				expression_source <- c(barName1,barName2)
				conts_base_pairs_percentage <- c((intersectionBPsum1/all_contsBPsum)*100,(intersectionBPsum2/all_contsBPsum)*100)
				FPKM_base_pairs_percentage <- c((intersectionBPsum1/all_fpkmBPsum1)*100,(intersectionBPsum2/all_fpkmBPsum2)*100)
				FPKM_percentage <- c((intersectionFPKMsum1/all_fpkmSum1)*100,(intersectionFPKMsum2/all_fpkmSum2)*100)
				
				conts_bp_plotData <- data.frame(expression_source,conts_base_pairs_percentage)
				fpkm_bp_plotData <- data.frame(expression_source,FPKM_base_pairs_percentage)
				fpkm_plotData <- data.frame(expression_source,FPKM_percentage)
				
				#creating the percentage of bp intersected from contact bands plot:
				plo_contsBP <- ggplot2::ggplot(conts_bp_plotData, ggplot2::aes(expression_source,conts_base_pairs_percentage,fill=expression_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of bp from contact bands intersected between expression and contacts of",sp1[1],sp1[2],sp1[3],"in the whole genome")) + ggplot2::labs(x="expression source",y="intersected bp percentage of contact bands(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				if(intersectionBPsum1>0 & intersectionBPsum2>0)
				{			
					#testing the significance of difference between percentages of both data groups
					cat("\n\ntesting the significance of the difference between the percentages of bp intersected of both data groups...\n")
					#prop_contsBP <- prop.test(conts_bp_plotData$conts_base_pairs_percentage,c(all_contsBPsum,all_contsBPsum),correct=FALSE)
					cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the contact bands\nthe 2nd creates a table of all the successes (bp intersections) and fails (all the bp of the contact bands that didn't intersect) and uses them\n\n")
					prop_contsBP1 <- prop.test(c(intersectionBPsum1,intersectionBPsum2),c(all_contsBPsum,all_contsBPsum),correct=FALSE)
					print(prop_contsBP1)
					prop_contsBP2 <- prop.test(matrix(c(all_contsBPsum-intersectionBPsum1,intersectionBPsum1,all_contsBPsum-intersectionBPsum2,intersectionBPsum2),ncol=2),correct=FALSE)
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
					compute.es::propes((intersectionBPsum1/all_contsBPsum),(intersectionBPsum2/all_contsBPsum),all_contsBPsum,all_contsBPsum)
					cat("\n")
					
					#phi
					cat("\nphi:\n")
					print(psych::phi(matrix(c(intersectionBPsum1,all_contsBPsum-intersectionBPsum1,intersectionBPsum2,all_contsBPsum-intersectionBPsum2),ncol=2,byrow=TRUE)))
					cat("\n")
					
					#creating a barplot with significance stars
					cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")					
					label.df <- data.frame(expression_source = conts_bp_plotData[which.max(conts_bp_plotData$conts_base_pairs_percentage),1],conts_base_pairs_percentage = conts_bp_plotData[which.max(conts_bp_plotData$conts_base_pairs_percentage),2]+1)
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
						cat("\nthe",barName1,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum2 == 0)
					{
						cat("\nthe",barName2,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum1 == 0 & intersectionBPsum2 == 0)
					{
						cat("\nthe",barName1,"and",barName2,"data do not have any intersections\n\n")
					}			
					cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")	
					print(plo_contsBP)
				}

				ans8 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans8 == "y")
				{
					nm <- paste("ContactsBP_FPKMintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file"
					
					bpORfpkm <- "bp"
					
					percentageOf <- "Contact bands"
					
					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}
					
					Intersection_chromosomes <- "whole genome"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(nm,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
				}
				
				#creating the percentage of bp intersected from RNA-seq bands plot:
				plo_fpkmBP <- ggplot2::ggplot(fpkm_bp_plotData, ggplot2::aes(expression_source,FPKM_base_pairs_percentage,fill=expression_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of bp from RNA-seq bands intersected between expression and contacts of",sp1[1],sp1[2],sp1[3],"in the whole genome")) + ggplot2::labs(x="expression source",y="intersected bp percentage RNA-seq bands(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				if(intersectionBPsum1>0 & intersectionBPsum2>0)
				{			
					#testing the significance of difference between percentages of both data groups
					cat("\n\ntesting the significance of the difference between the percentages of bp intersected out of the bp of the FPKM bands, for both data groups...\n\n")
					#prop_fpkmBP <- prop.test(fpkm_bp_plotData$FPKM_base_pairs_percentage,c(all_fpkmBPsum1,all_fpkmBPsum2),correct=FALSE)
					cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the FPKM file\nthe 2nd creates a table of all the successes (bp of FPKMs intersections) and fails (all the bp of the FPKM bands that didn't intersect) and uses them\n\n")					
					prop_fpkmBP1 <- prop.test(c(intersectionBPsum1,intersectionBPsum2),c(all_fpkmBPsum1,all_fpkmBPsum2),correct=FALSE)
					print(prop_fpkmBP1)
					prop_fpkmBP2 <- prop.test(matrix(c(all_fpkmBPsum1-intersectionBPsum1,intersectionBPsum1,all_fpkmBPsum2-intersectionBPsum2,intersectionBPsum2),ncol=2),correct=FALSE)
					print(prop_fpkmBP2)
					prop_ans2 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
					if(prop_ans2 == 1)
					{
						prop_fpkmBP <- prop_fpkmBP1
					}
					else
					{
						prop_fpkmBP <- prop_fpkmBP2
					}
					
					#calculating effect sizes
					cat("\neffect sizes:\n")
					compute.es::propes((intersectionBPsum1/all_fpkmBPsum1),(intersectionBPsum2/all_fpkmBPsum2),all_fpkmBPsum1,all_fpkmBPsum2)
					cat("\n")
					
					#phi
					cat("\nphi:\n")
					print(psych::phi(matrix(c(intersectionBPsum1,all_fpkmBPsum1-intersectionBPsum1,intersectionBPsum2,all_fpkmBPsum2-intersectionBPsum2),ncol=2,byrow=TRUE)))
					cat("\n")
					
					#creating a barplot with significance stars	
					cat("\nplotting the FPKM bands bp percentage data\nplease wait...\n\n")					
					label.df <- data.frame(expression_source = fpkm_bp_plotData[which.max(fpkm_bp_plotData$FPKM_base_pairs_percentage),1],FPKM_base_pairs_percentage = fpkm_bp_plotData[which.max(fpkm_bp_plotData$FPKM_base_pairs_percentage),2]+1)
					if(prop_fpkmBP$p.value <= 0.0001)
					{
						print(plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "****",size=8))
					}
					else if(prop_fpkmBP$p.value <= 0.001)
					{
						print(plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "***",size=8))
					}
					else if(prop_fpkmBP$p.value <= 0.01)
					{
						print(plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "**",size=8))
					}
					else if(prop_fpkmBP$p.value <= 0.05)
					{
						print(plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "*",size=8))
					}
					else
					{
						print(plo_fpkmBP)
					}
				}
				else
				{
					if(intersectionBPsum1 == 0)
					{
						cat("\nthe",barName1,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum2 == 0)
					{
						cat("\nthe",barName2,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum1 == 0 & intersectionBPsum2 == 0)
					{
						cat("\nthe",barName1,"and",barName2,"data do not have any intersections\n\n")
					}			
					cat("\nplotting the FPKM bands bp percentage data\nplease wait...\n\n")	
					print(plo_fpkmBP)
				}

				ans5 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans5 == "y")
				{
					nm <- paste("fpkmBP_FPKMintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file"
					
					bpORfpkm <- "bp"
					
					percentageOf <- "RNA-seq fragments"

					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}
					
					Intersection_chromosomes <- "whole genome"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(nm,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
				}
				
				#creating the percentage of FPKM intersected of sum of FPKM plot:
				ploFPKM <- ggplot2::ggplot(fpkm_plotData, ggplot2::aes(expression_source,FPKM_percentage,fill=expression_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected FPKMs between expression and contacts of",sp1[1],sp1[2],sp1[3],"in the whole genome")) + ggplot2::labs(x="expression source",y="intersected FPKM percentage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				if(intersectionFPKMsum1>0 & intersectionFPKMsum2>0)
				{
					#testing the significance of difference between percentages of both data groups
					cat("\n\ntesting the significance of the difference between the percentages of FPKM intersected of both data groups...\n\n")
					#propFPKM <- prop.test(fpkm_plotData$FPKM_percentage,c(all_fpkmSum1,all_fpkmSum2),correct=FALSE)
					cat("\nthe 1st is the comparison of proportions of intersections out of all FPKMs in the FPKM bands\nthe 2nd creates a table of all the successes (FPKM intersections) and fails (all the FPKM of the FPKM bands that didn't intersect) and uses them\n\n")					
					propFPKM1 <- prop.test(c(intersectionFPKMsum1,intersectionFPKMsum2),c(all_fpkmSum1,all_fpkmSum2),correct=FALSE)
					print(propFPKM1)
					propFPKM2 <- prop.test(matrix(c(all_fpkmSum1-intersectionFPKMsum1,intersectionFPKMsum1,all_fpkmSum2-intersectionFPKMsum2,intersectionFPKMsum2),ncol=2),correct=FALSE)
					print(propFPKM2)					
					prop_ans3 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
					if(prop_ans3 == 1)
					{
						propFPKM <- propFPKM1
					}
					else
					{
						propFPKM <- propFPKM2
					}
					
					#calculating effect sizes
					cat("\neffect sizes:\n")
					compute.es::propes((intersectionFPKMsum1/all_fpkmSum1),(intersectionFPKMsum2/all_fpkmSum2),all_fpkmSum1,all_fpkmSum2)
					cat("\n")
					
					#phi
					cat("\nphi:\n")
					print(psych::phi(matrix(c(intersectionFPKMsum1,all_fpkmSum1-intersectionFPKMsum1,intersectionFPKMsum2,all_fpkmSum2-intersectionFPKMsum2),ncol=2,byrow=TRUE)))
					cat("\n")
					
					#creating a barplot with significance stars
					cat("\nplotting the FPKM percentage data\nplease wait...\n\n")					
					label.df <- data.frame(expression_source = fpkm_plotData[which.max(fpkm_plotData$FPKM_percentage),1],FPKM_percentage = fpkm_plotData[which.max(fpkm_plotData$FPKM_percentage),2]+1)	
					if(propFPKM$p.value <= 0.0001)
					{
						print(ploFPKM + ggplot2::geom_text(data = label.df, label = "****",size=8))
					}
					else if(propFPKM$p.value <= 0.001)
					{
						print(ploFPKM + ggplot2::geom_text(data = label.df, label = "***",size=8))
					}
					else if(propFPKM$p.value <= 0.01)
					{
						print(ploFPKM + ggplot2::geom_text(data = label.df, label = "**",size=8))
					}
					else if(propFPKM$p.value <= 0.05)
					{
						print(ploFPKM + ggplot2::geom_text(data = label.df, label = "*",size=8))
					}
					else
					{
						print(ploFPKM)
					}
				}	
				else
				{
					if(intersectionFPKMsum1 == 0)
					{
						cat("\nthe",barName1,"data does not have any intersections\n\n")
					}
					else if(intersectionFPKMsum2 == 0)
					{
						cat("\nthe",barName2,"data does not have any intersections\n\n")
					}
					else if(intersectionFPKMsum1 == 0 & intersectionFPKMsum2 == 0)
					{
						cat("\nthe",barName1,"and",barName2,"data do not have any intersections\n\n")
					}				
					cat("\nplotting the FPKM percentage data\nplease wait...\n\n")
					print(ploFPKM)
				}
				
				ans6 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans6 == "y")
				{
					nm <- paste("fpkm_FPKMintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file"
					
					bpORfpkm <- "fpkm"
					
					percentageOf <- "FPKM"
					
					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "whole genome"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(nm,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}
			}
			
			#venn diagrams
			if(readline(prompt=cat("\nwould you like to create venn diagrams?\ny/n\n\n")) == "y")
			{
				#creating the venn data frame for 0 percent
				venn_DF <- data.frame(rep(1,nrow(conts)),matrix(0,nrow(conts),2))##!

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
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM1_afterCO.bed -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM2_afterCO.bed -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A2.bed",sep=""))
					venn_A1 <- read.table(paste("~/Analyze4C/temp/venn_A1.bed",sep=""))
					venn_A2 <- read.table(paste("~/Analyze4C/temp/venn_A2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn_A1.bed")
					system("rm ~/Analyze4C/temp/venn_A2.bed")
					if(venn_ans2 == "y")
					{
						system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM1_afterCO.bed -c > ~/Analyze4C/temp/venn_B1.bed",sep=""))
						system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM2_afterCO.bed -c > ~/Analyze4C/temp/venn_B2.bed",sep=""))
						venn_B1 <- read.table(paste("~/Analyze4C/temp/venn_B1.bed",sep=""))
						venn_B2 <- read.table(paste("~/Analyze4C/temp/venn_B2.bed",sep=""))
						#removing files
						system("rm ~/Analyze4C/temp/venn_B1.bed")
						system("rm ~/Analyze4C/temp/venn_B2.bed")
						#if the number in the 1st is 0 but the number in the 2nd is more than ovlp_num then we will consider it to be 1
						venn_A1[venn_B1[,4]>=ovlp_num,4] <- 1
						venn_A2[venn_B2[,4]>=ovlp_num,4] <- 1
					}
					#add the correct number to each column
					venn_DF[venn_A1[,4]>=1,2] <- 1
					venn_DF[venn_A2[,4]>=1,3] <- 1
				}
				else
				{
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM1_afterCO.bed -c > ~/Analyze4C/temp/venn1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM2_afterCO.bed -c > ~/Analyze4C/temp/venn2.bed",sep=""))
					venn1 <- read.table(paste("~/Analyze4C/temp/venn1.bed",sep=""))
					venn2 <- read.table(paste("~/Analyze4C/temp/venn2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn1.bed")
					system("rm ~/Analyze4C/temp/venn2.bed")
					#add the correct number to each column
					venn_DF[venn1[,4]>=1,2] <- 1
					venn_DF[venn2[,4]>=1,3] <- 1				
				}
				
				# create venn diagram using the data frame venn_DF
				venn_flag <- 0 #if venn_flag is 1 then we save a diagram and it should be recorded in 'expressionVScontacts_sumOFintersections_plots'
				vennName <- paste("venn_psContsVS2FPKMs_",DandT2,".jpg",sep="")
				venn_flag <- venn_creator(venn_DF,2,1,vennName)			

				#recording the data into 'expressionVScontacts_sumOFintersections_plots.txt'
				if(venn_flag == 1)
				{					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file - Venn Diagram"					
					bpORfpkm <- NA
					percentageOf <- NA
					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "whole genome"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(vennName,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}	
			}					
		}
		else if(ans3 == 2) #trans
		{
			#the mcnemar test
			ans.mcnemar <- readline("\nwould you like to perform a mcnemar test?\ny/n\n\n")
			if(ans.mcnemar == "y")
			{
				#we are only taking the trans since we are looking at the proportion that intersects, which in this case is only from trans
				FPKM1_afterCO_mcnemar <- FPKM1_afterCO[FPKM1_afterCO[,1]!=cis,]
				FPKM2_afterCO_mcnemar <- FPKM2_afterCO[FPKM2_afterCO[,1]!=cis,]
				conts_mcnemar <- conts[conts[,1]!=cis,]
				write.table(FPKM1_afterCO_mcnemar,"~/Analyze4C/temp/FPKM1_afterCO_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
				write.table(FPKM2_afterCO_mcnemar,"~/Analyze4C/temp/FPKM2_afterCO_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
				write.table(conts_mcnemar,"~/Analyze4C/temp/conts_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
	
				mcnemar_tester("~/Analyze4C/temp/conts_mcnemar.bed","~/Analyze4C/temp/FPKM1_afterCO_mcnemar.bed","~/Analyze4C/temp/FPKM2_afterCO_mcnemar.bed")
				
				system("rm ~/Analyze4C/temp/FPKM1_afterCO_mcnemar.bed")
				system("rm ~/Analyze4C/temp/FPKM2_afterCO_mcnemar.bed")	
				system("rm ~/Analyze4C/temp/conts_mcnemar.bed")								
			}
		
			ans.prop <- readline(prompt=cat("\n\nwould you like to perform a prop test on the intersected bp and the intersected FPKM?\ny/n\n\n"))
			if(ans.prop == "y")
			{
				#BPs:
				
				#getting the sum of all the bps of the intersections in trans
				intersectionBPsum1 <- sum(expressionComparison1[expressionComparison1[,1]!=cis,8])
				intersectionBPsum2 <- sum(expressionComparison2[expressionComparison2[,1]!=cis,8])
				
				#getting the sum of all the bps of the contact bands in trans
				all_contsBPsum <- sum(conts[conts[,1]!=cis,3]-conts[conts[,1]!=cis,2])
				
				#getting the sum of all the bps that are above the FPKM cutoff in trans
				all_fpkmBPsum1 <- sum(FPKM1_afterCO[FPKM1_afterCO[,1]!=cis,3]-FPKM1_afterCO[FPKM1_afterCO[,1]!=cis,2])
				all_fpkmBPsum2 <- sum(FPKM2_afterCO[FPKM2_afterCO[,1]!=cis,3]-FPKM2_afterCO[FPKM2_afterCO[,1]!=cis,2])
				
				#FPKMs:

				#getting the sum of all the FPKMs of the intersections
				expressionComparison1_trans <- expressionComparison1[row.names(unique(expressionComparison1[,c(colnames(expressionComparison1)[1],colnames(expressionComparison1)[2])])),]
				intersectionFPKMsum1 <- sum(expressionComparison1_trans[expressionComparison1_trans[,1]!=cis,4])
				expressionComparison2_trans <- expressionComparison2[row.names(unique(expressionComparison2[,c(colnames(expressionComparison2)[1],colnames(expressionComparison2)[2])])),]
				intersectionFPKMsum2 <- sum(expressionComparison2_trans[expressionComparison2_trans[,1]!=cis,4])
				
				#getting the sum of all the FPKMs above the FPKM cutoff in trans
				all_fpkmSum1 <- sum(FPKM1_afterCO[FPKM1_afterCO[,1]!=cis,4])
				all_fpkmSum2 <- sum(FPKM2_afterCO[FPKM2_afterCO[,1]!=cis,4])
				
				cat("\nthe contact bands file chosen:\n",conts_filename,"\n")
				cat("\nthe number of bps of contact bands in trans:",all_contsBPsum,"\n")
				cat("\nthe number of bps of FPKMs above the cutoff in the trans of the expression file",FPKM_filename1,"is:",all_fpkmBPsum1,"\n")
				cat("\nthe number of bps of FPKMs above the cutoff in the trans of the expression file",FPKM_filename2,"is:",all_fpkmBPsum2,"\n")
				cat("\nthe sum of FPKMs above the cutoff in the trans of the expression file",FPKM_filename1,"is:",all_fpkmSum1,"\n")
				cat("\nthe sum of FPKMs above the cutoff in the trans of the expression file",FPKM_filename2,"is:",all_fpkmSum2,"\n")			
				cat("\nexpression file - ",FPKM_filename1,":\nthe number of bps intersected in trans are ",intersectionBPsum1," which are ",(intersectionBPsum1/all_contsBPsum)*100,"% of all the bps of the contact bands in trans\n\n",sep="")
				cat("\nexpression file - ",FPKM_filename2,":\nthe number of bps intersected in trans are ",intersectionBPsum2," which are ",(intersectionBPsum2/all_contsBPsum)*100,"% of all the bps of the contact bands in trans\n\n",sep="")			
				cat("\nexpression file - ",FPKM_filename1,":\nthe number of bps intersected in trans are ",intersectionBPsum1," which are ",(intersectionBPsum1/all_fpkmBPsum1)*100,"% of all the bps of FPKMs above the cutoff in trans\n\n",sep="")
				cat("\nexpression file - ",FPKM_filename2,":\nthe number of bps intersected in trans are ",intersectionBPsum2," which are ",(intersectionBPsum2/all_fpkmBPsum2)*100,"% of all the bps of FPKMs above the cutoff in trans\n\n",sep="")
				cat("\nexpression file - ",FPKM_filename1,":\nthe sum of FPKMs intersected in trans are ",intersectionFPKMsum1," which are ",(intersectionFPKMsum1/all_fpkmSum1)*100,"% of all the FPKMs above the cutoff in trans\n\n",sep="")
				cat("\nexpression file - ",FPKM_filename2,":\nthe sum of FPKMs intersected in trans are ",intersectionFPKMsum2," which are ",(intersectionFPKMsum2/all_fpkmSum2)*100,"% of all the FPKMs above the cutoff in trans\n\n",sep="")		

				#creating plot		
				barName1 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",FPKM_filename1,"on the plot:\n"))
				barName2 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",FPKM_filename2,"on the plot:\n"))
				expression_source <- c(barName1,barName2)
				conts_base_pairs_percentage <- c((intersectionBPsum1/all_contsBPsum)*100,(intersectionBPsum2/all_contsBPsum)*100)
				FPKM_base_pairs_percentage <- c((intersectionBPsum1/all_fpkmBPsum1)*100,(intersectionBPsum2/all_fpkmBPsum2)*100)
				FPKM_percentage <- c((intersectionFPKMsum1/all_fpkmSum1)*100,(intersectionFPKMsum2/all_fpkmSum2)*100)
				
				conts_bp_plotData <- data.frame(expression_source,conts_base_pairs_percentage)
				fpkm_bp_plotData <- data.frame(expression_source,FPKM_base_pairs_percentage)
				fpkm_plotData <- data.frame(expression_source,FPKM_percentage)
				
				#creating the percentage of bp intersected from contact bands plot:
				plo_contsBP <- ggplot2::ggplot(conts_bp_plotData, ggplot2::aes(expression_source,conts_base_pairs_percentage,fill=expression_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of bp from contact bands intersected between expression and contacts of",sp1[1],sp1[2],sp1[3],"in trans")) + ggplot2::labs(x="expression source",y="intersected bp percentage of contact bands(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				if(intersectionBPsum1>0 & intersectionBPsum2>0)
				{		
					#testing the significance of difference between percentages of both data groups
					cat("\n\ntesting the significance of the difference between the percentages of bp intersected of both data groups...\n")
					#prop_contsBP <- prop.test(conts_bp_plotData$conts_base_pairs_percentage,c(all_contsBPsum,all_contsBPsum),correct=FALSE)
					cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the contact bands\nthe 2nd creates a table of all the successes (bp intersections) and fails (all the bp of the contact bands that didn't intersect) and uses them\n\n")
					prop_contsBP1 <- prop.test(c(intersectionBPsum1,intersectionBPsum2),c(all_contsBPsum,all_contsBPsum),correct=FALSE)
					print(prop_contsBP1)
					prop_contsBP2 <- prop.test(matrix(c(all_contsBPsum-intersectionBPsum1,intersectionBPsum1,all_contsBPsum-intersectionBPsum2,intersectionBPsum2),ncol=2),correct=FALSE)
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
					compute.es::propes((intersectionBPsum1/all_contsBPsum),(intersectionBPsum2/all_contsBPsum),all_contsBPsum,all_contsBPsum)
					cat("\n")
					
					#phi
					cat("\nphi:\n")
					print(psych::phi(matrix(c(intersectionBPsum1,all_contsBPsum-intersectionBPsum1,intersectionBPsum2,all_contsBPsum-intersectionBPsum2),ncol=2,byrow=TRUE)))
					cat("\n")
					
					#creating a barplot with significance stars
					cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")					
					label.df <- data.frame(expression_source = conts_bp_plotData[which.max(conts_bp_plotData$conts_base_pairs_percentage),1],conts_base_pairs_percentage = conts_bp_plotData[which.max(conts_bp_plotData$conts_base_pairs_percentage),2]+1)					
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
						cat("\nthe",barName1,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum2 == 0)
					{
						cat("\nthe",barName2,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum1 == 0 & intersectionBPsum2 == 0)
					{
						cat("\nthe",barName1,"and",barName2,"data do not have any intersections\n\n")
					}			
					cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")	
					print(plo_contsBP)
				}

				ans8 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans8 == "y")
				{
					nm <- paste("ContactsBP_FPKMintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file"
					
					bpORfpkm <- "bp"
					
					percentageOf <- "Contact bands"
					
					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}
					
					Intersection_chromosomes <- "trans"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(nm,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
				}
				
				#creating the percentage of bp from RNA-seq bands plot:
				plo_fpkmBP <- ggplot2::ggplot(fpkm_bp_plotData, ggplot2::aes(expression_source,FPKM_base_pairs_percentage,fill=expression_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of bp from RNA-seq bands intersected between expression and contacts of",sp1[1],sp1[2],sp1[3],"in trans")) + ggplot2::labs(x="expression source",y="intersected bp percentage RNA-seq bands(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				if(intersectionBPsum1>0 & intersectionBPsum2>0)
				{				
					#testing the significance of difference between percentages of both data groups
					cat("\n\ntesting the significance of the difference between the percentages of bp intersected out of the bp of the FPKM bands, for both data groups...\n\n")
					#prop_fpkmBP <- prop.test(fpkm_bp_plotData$FPKM_base_pairs_percentage,c(all_fpkmBPsum1,all_fpkmBPsum2),correct=FALSE)
					cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the FPKM file\nthe 2nd creates a table of all the successes (bp of FPKMs intersections) and fails (all the bp of the FPKM bands that didn't intersect) and uses them\n\n")					
					prop_fpkmBP1 <- prop.test(c(intersectionBPsum1,intersectionBPsum2),c(all_fpkmBPsum1,all_fpkmBPsum2),correct=FALSE)
					print(prop_fpkmBP1)
					prop_fpkmBP2 <- prop.test(matrix(c(all_fpkmBPsum1-intersectionBPsum1,intersectionBPsum1,all_fpkmBPsum2-intersectionBPsum2,intersectionBPsum2),ncol=2),correct=FALSE)
					print(prop_fpkmBP2)
					prop_ans2 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
					if(prop_ans2 == 1)
					{
						prop_fpkmBP <- prop_fpkmBP1
					}
					else
					{
						prop_fpkmBP <- prop_fpkmBP2
					}
					
					#calculating effect sizes
					cat("\neffect sizes:\n")
					compute.es::propes((intersectionBPsum1/all_fpkmBPsum1),(intersectionBPsum2/all_fpkmBPsum2),all_fpkmBPsum1,all_fpkmBPsum2)
					cat("\n")
					
					#phi
					cat("\nphi:\n")
					print(psych::phi(matrix(c(intersectionBPsum1,all_fpkmBPsum1-intersectionBPsum1,intersectionBPsum2,all_fpkmBPsum2-intersectionBPsum2),ncol=2,byrow=TRUE)))
					cat("\n")
					
					#creating a barplot with significance stars
					cat("\nplotting the FPKM bands bp percentage data\nplease wait...\n\n")					
					label.df <- data.frame(expression_source = fpkm_bp_plotData[which.max(fpkm_bp_plotData$FPKM_base_pairs_percentage),1],FPKM_base_pairs_percentage = fpkm_bp_plotData[which.max(fpkm_bp_plotData$FPKM_base_pairs_percentage),2]+1)
					if(prop_fpkmBP$p.value <= 0.0001)
					{
						print(plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "****",size=8))
					}
					else if(prop_fpkmBP$p.value <= 0.001)
					{
						print(plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "***",size=8))
					}
					else if(prop_fpkmBP$p.value <= 0.01)
					{
						print(plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "**",size=8))
					}
					else if(prop_fpkmBP$p.value <= 0.05)
					{
						print(plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "*",size=8))
					}
					else
					{
						print(plo_fpkmBP)
					}
				}
				else
				{
					if(intersectionBPsum1 == 0)
					{
						cat("\nthe",barName1,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum2 == 0)
					{
						cat("\nthe",barName2,"data does not have any intersections\n\n")
					}
					else if(intersectionBPsum1 == 0 & intersectionBPsum2 == 0)
					{
						cat("\nthe",barName1,"and",barName2,"data do not have any intersections\n\n")
					}			
					cat("\nplotting the FPKM bands bp percentage data\nplease wait...\n\n")
					print(plo_fpkmBP)
				}

				ans5 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans5 == "y")
				{
					nm <- paste("fpkmBP_FPKMintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file"
					
					bpORfpkm <- "bp"
					
					percentageOf <- "RNA-seq fragments"

					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}
					
					Intersection_chromosomes <- "trans"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(nm,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
				}
				
				#creating the percentage of FPKM intersected of sum of FPKM plot:
				ploFPKM <- ggplot2::ggplot(fpkm_plotData, ggplot2::aes(expression_source,FPKM_percentage,fill=expression_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected FPKMs between expression and contacts of",sp1[1],sp1[2],sp1[3],"in trans")) + ggplot2::labs(x="expression source",y="intersected FPKM percentage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				if(intersectionFPKMsum1>0 & intersectionFPKMsum2>0)
				{
					#testing the significance of difference between percentages of both data groups
					cat("\n\ntesting the significance of the difference between the percentages of FPKM intersected of both data groups...\n\n")
					#propFPKM <- prop.test(fpkm_plotData$FPKM_percentage,c(all_fpkmSum1,all_fpkmSum2),correct=FALSE)
					cat("\nthe 1st is the comparison of proportions of intersections out of all FPKMs in the FPKM bands\nthe 2nd creates a table of all the successes (FPKM intersections) and fails (all the FPKM of the FPKM bands that didn't intersect) and uses them\n\n")					
					propFPKM1 <- prop.test(c(intersectionFPKMsum1,intersectionFPKMsum2),c(all_fpkmSum1,all_fpkmSum2),correct=FALSE)
					print(propFPKM1)
					propFPKM2 <- prop.test(matrix(c(all_fpkmSum1-intersectionFPKMsum1,intersectionFPKMsum1,all_fpkmSum2-intersectionFPKMsum2,intersectionFPKMsum2),ncol=2),correct=FALSE)
					print(propFPKM2)					
					prop_ans3 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
					if(prop_ans3 == 1)
					{
						propFPKM <- propFPKM1
					}
					else
					{
						propFPKM <- propFPKM2
					}
					
					#calculating effect sizes
					cat("\neffect sizes:\n")
					compute.es::propes((intersectionFPKMsum1/all_fpkmSum1),(intersectionFPKMsum2/all_fpkmSum2),all_fpkmSum1,all_fpkmSum2)
					cat("\n")
					
					#phi
					cat("\nphi:\n")
					print(psych::phi(matrix(c(intersectionFPKMsum1,all_fpkmSum1-intersectionFPKMsum1,intersectionFPKMsum2,all_fpkmSum2-intersectionFPKMsum2),ncol=2,byrow=TRUE)))
					cat("\n")
					
					#creating a barplot with significance stars					
					cat("\nplotting the FPKM percentage data\nplease wait...\n\n")
					label.df <- data.frame(expression_source = fpkm_plotData[which.max(fpkm_plotData$FPKM_percentage),1],FPKM_percentage = fpkm_plotData[which.max(fpkm_plotData$FPKM_percentage),2]+1)
					if(propFPKM$p.value <= 0.0001)
					{
						print(ploFPKM + ggplot2::geom_text(data = label.df, label = "****",size=8))
					}
					else if(propFPKM$p.value <= 0.001)
					{
						print(ploFPKM + ggplot2::geom_text(data = label.df, label = "***",size=8))
					}
					else if(propFPKM$p.value <= 0.01)
					{
						print(ploFPKM + ggplot2::geom_text(data = label.df, label = "**",size=8))
					}
					else if(propFPKM$p.value <= 0.05)
					{
						print(ploFPKM + ggplot2::geom_text(data = label.df, label = "*",size=8))
					}
					else
					{
						print(ploFPKM)
					}
				}	
				else
				{
					if(intersectionFPKMsum1 == 0)
					{
						cat("\nthe",barName1,"data does not have any intersections\n\n")
					}
					else if(intersectionFPKMsum2 == 0)
					{
						cat("\nthe",barName2,"data does not have any intersections\n\n")
					}
					else if(intersectionFPKMsum1 == 0 & intersectionFPKMsum2 == 0)
					{
						cat("\nthe",barName1,"and",barName2,"data do not have any intersections\n\n")
					}				
					cat("\nplotting the FPKM percentage data\nplease wait...\n\n")
					print(ploFPKM)
				}
				
				ans6 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans6 == "y")
				{
					nm <- paste("fpkm_FPKMintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file"
					
					bpORfpkm <- "fpkm"

					percentageOf <- "FPKM"
					
					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "trans"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(nm,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}
			}

			#venn diagrams
			if(readline(prompt=cat("\nwould you like to create venn diagrams?\ny/n\n\n")) == "y")
			{
				#creating the venn data frame for 0 percent
				venn_DF <- data.frame(rep(1,nrow(conts)),matrix(0,nrow(conts),2))##!

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
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM1_afterCO.bed -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM2_afterCO.bed -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A2.bed",sep=""))
					venn_A1 <- read.table(paste("~/Analyze4C/temp/venn_A1.bed",sep=""))
					venn_A2 <- read.table(paste("~/Analyze4C/temp/venn_A2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn_A1.bed")
					system("rm ~/Analyze4C/temp/venn_A2.bed")
					if(venn_ans2 == "y")
					{
						system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM1_afterCO.bed -c > ~/Analyze4C/temp/venn_B1.bed",sep=""))
						system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM2_afterCO.bed -c > ~/Analyze4C/temp/venn_B2.bed",sep=""))
						venn_B1 <- read.table(paste("~/Analyze4C/temp/venn_B1.bed",sep=""))
						venn_B2 <- read.table(paste("~/Analyze4C/temp/venn_B2.bed",sep=""))
						#removing files
						system("rm ~/Analyze4C/temp/venn_B1.bed")
						system("rm ~/Analyze4C/temp/venn_B2.bed")
						#if the number in the 1st is 0 but the number in the 2nd is more than ovlp_num then we will consider it to be 1
						venn_A1[venn_B1[,4]>=ovlp_num,4] <- 1
						venn_A2[venn_B2[,4]>=ovlp_num,4] <- 1
					}
					#add the correct number to each column
					venn_DF[venn_A1[,4]>=1 & venn_A1[,1]!=cis,2] <- 1
					venn_DF[venn_A2[,4]>=1 & venn_A2[,1]!=cis,3] <- 1
				}
				else
				{
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM1_afterCO.bed -c > ~/Analyze4C/temp/venn1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM2_afterCO.bed -c > ~/Analyze4C/temp/venn2.bed",sep=""))
					venn1 <- read.table(paste("~/Analyze4C/temp/venn1.bed",sep=""))
					venn2 <- read.table(paste("~/Analyze4C/temp/venn2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn1.bed")
					system("rm ~/Analyze4C/temp/venn2.bed")
					#add the correct number to each column
					venn_DF[venn1[,4]>=1 & venn1[,1]!=cis,2] <- 1
					venn_DF[venn2[,4]>=1 & venn2[,1]!=cis,3] <- 1				
				}
				
				# create venn diagram using the data frame venn_DF
				venn_flag <- 0 #if venn_flag is 1 then we save a diagram and it should be recorded in 'expressionVScontacts_sumOFintersections_plots'
				vennName <- paste("venn_psContsVS2FPKMs_",DandT2,".jpg",sep="")
				venn_flag <- venn_creator(venn_DF,2,1,vennName)			

				#recording the data into 'expressionVScontacts_sumOFintersections_plots.txt'
				if(venn_flag == 1)
				{					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file - Venn Diagram"					
					bpORfpkm <- NA
					percentageOf <- NA
					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "trans"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(vennName,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}	
			}							
		}	
		else if(ans3 == 3) #chromosomes separately
		{
			#the mcnemar test
			ans.mcnemar <- readline("\nwould you like to perform a mcnemar test?\ny/n\n\n")
			if(ans.mcnemar == "y")
			{
				for(u in 1:numOFchroms)
				{	
					cat("\nmcnemar test - chromosome ",u,":\n\n",sep="")
					FPKM1_afterCO_mcnemar <- FPKM1_afterCO[FPKM1_afterCO[,1]==u,]
					FPKM2_afterCO_mcnemar <- FPKM2_afterCO[FPKM2_afterCO[,1]==u,]
					conts_mcnemar <- conts[conts[,1]==u,]
					write.table(FPKM1_afterCO_mcnemar,"~/Analyze4C/temp/FPKM1_afterCO_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
					write.table(FPKM2_afterCO_mcnemar,"~/Analyze4C/temp/FPKM2_afterCO_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
					write.table(conts_mcnemar,"~/Analyze4C/temp/conts_mcnemar.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
		
					mcnemar_tester("~/Analyze4C/temp/conts_mcnemar.bed","~/Analyze4C/temp/FPKM1_afterCO_mcnemar.bed","~/Analyze4C/temp/FPKM2_afterCO_mcnemar.bed")
					
					system("rm ~/Analyze4C/temp/FPKM1_afterCO_mcnemar.bed")
					system("rm ~/Analyze4C/temp/FPKM2_afterCO_mcnemar.bed")	
					system("rm ~/Analyze4C/temp/conts_mcnemar.bed")
				}	
			}
		
			ans.prop <- readline(prompt=cat("\n\nwould you like to perform a prop test on the intersected bp and the intersected FPKM?\ny/n\n\n"))
			if(ans.prop == "y")
			{
				#intitiating these vectors with the size of the number of chromosomes, in order to put the sums per chromosome in them
				intersectionBPsum1 <- rep(0,numOFchroms)
				intersectionBPsum2 <- rep(0,numOFchroms)
				all_contsBPsum <- rep(0,numOFchroms)	
				all_fpkmBPsum1 <- rep(0,numOFchroms)
				all_fpkmBPsum2 <- rep(0,numOFchroms)
				intersectionFPKMsum1 <- rep(0,numOFchroms)
				intersectionFPKMsum2 <- rep(0,numOFchroms)
				all_fpkmSum1 <- rep(0,numOFchroms)
				all_fpkmSum2 <- rep(0,numOFchroms)
				contsBP_intersectedPercent1 <- rep(0,numOFchroms)
				contsBP_intersectedPercent2 <- rep(0,numOFchroms)
				fpkmBP_intersectedPercent1 <- rep(0,numOFchroms)
				fpkmBP_intersectedPercent2 <- rep(0,numOFchroms)			
				fpkm_intersectedPercent1 <- rep(0,numOFchroms)
				fpkm_intersectedPercent2 <- rep(0,numOFchroms)
				
				for(w in 1:numOFchroms)
				{
					#BPs:
					
					#getting the sum of all the bps of the intersections per chromosome
					intersectionBPsum1[w] <- sum(expressionComparison1[expressionComparison1[,1]==w,8])
					intersectionBPsum2[w] <- sum(expressionComparison2[expressionComparison2[,1]==w,8])
					
					#getting the sum of all the bps of the contact bands per chromosome
					all_contsBPsum[w] <- sum(conts[conts[,1]==w,3]-conts[conts[,1]==w,2])
					
					#getting the sum of all the bps that are above the FPKM cutoff per chromosome
					all_fpkmBPsum1[w] <- sum(FPKM1_afterCO[FPKM1_afterCO[,1]==w,3]-FPKM1_afterCO[FPKM1_afterCO[,1]==w,2])
					all_fpkmBPsum2[w] <- sum(FPKM2_afterCO[FPKM2_afterCO[,1]==w,3]-FPKM2_afterCO[FPKM2_afterCO[,1]==w,2])
					
					#FPKMs:

					#getting the sum of all the FPKMs of the intersections per chromosome
					expressionComparison1_chr <- expressionComparison1[row.names(unique(expressionComparison1[,c(colnames(expressionComparison1)[1],colnames(expressionComparison1)[2])])),]
					intersectionFPKMsum1[w] <- sum(expressionComparison1_chr[expressionComparison1_chr[,1]==w,4])
					expressionComparison2_chr <- expressionComparison2[row.names(unique(expressionComparison2[,c(colnames(expressionComparison2)[1],colnames(expressionComparison2)[2])])),]
					intersectionFPKMsum2[w] <- sum(expressionComparison2_chr[expressionComparison2_chr[,1]==w,4])
					
					#getting the sum of all the FPKMs above the FPKM cutoff per chromosome
					all_fpkmSum1[w] <- sum(FPKM1_afterCO[FPKM1_afterCO[,1]==w,4])
					all_fpkmSum2[w] <- sum(FPKM2_afterCO[FPKM2_afterCO[,1]==w,4])

					contsBP_intersectedPercent1[w] <- (intersectionBPsum1[w]/all_contsBPsum[w])*100
					contsBP_intersectedPercent2[w] <- (intersectionBPsum2[w]/all_contsBPsum[w])*100				
					fpkmBP_intersectedPercent1[w] <- (intersectionBPsum1[w]/all_fpkmBPsum1[w])*100
					fpkmBP_intersectedPercent2[w] <- (intersectionBPsum2[w]/all_fpkmBPsum2[w])*100
					fpkm_intersectedPercent1[w] <- (intersectionFPKMsum1[w]/all_fpkmSum1[w])*100
					fpkm_intersectedPercent2[w] <- (intersectionFPKMsum2[w]/all_fpkmSum2[w])*100
				}
				
				#printing the results:
				
				cat("\nthe contact bands file chosen:\n",conts_filename,"\n")
				for(s in 1:numOFchroms)
				{
					cat("\nthe number of bps of contact bands in chromosome",s,"is:",all_contsBPsum[s],"\n")
					cat("\nthe number of bps of FPKMs above the cutoff in the expression file",FPKM_filename1,"in chromosome",s,"is:",all_fpkmBPsum1[s],"\n")
					cat("\nthe number of bps of FPKMs above the cutoff in the expression file",FPKM_filename2,"in chromosome",s,"is:",all_fpkmBPsum2[s],"\n")
					cat("\nthe sum of FPKMs above the cutoff in the expression file",FPKM_filename1,"in chromosome",s,"is:",all_fpkmSum1[s],"\n")
					cat("\nthe sum of FPKMs above the cutoff in the expression file",FPKM_filename2,"in chromosome",s,"is:",all_fpkmSum2[s],"\n")					
				}
				
				#percentage of all bps in contact bands
				cat("\nexpression file - ",FPKM_filename1,":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis)
					{
						cat("chr",q," (cis): the number of bps intersected are ",intersectionBPsum1[q]," which are ",contsBP_intersectedPercent1[q],"% of all the bps of contact bands in the chromosome\n",sep="")							
					}
					else
					{
						cat("chr",q,": the number of bps intersected are ",intersectionBPsum1[q]," which are ",contsBP_intersectedPercent1[q],"% of all the bps contact bands in the chromosome\n",sep="")				
					}
				}
				cat("\n")
				
				cat("\nexpression file - ",FPKM_filename2,":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis)
					{
						cat("chr",q," (cis): the number of bps intersected are ",intersectionBPsum2[q]," which are ",contsBP_intersectedPercent2[q],"% of all the bps of contact bands in the chromosome\n",sep="")							
					}
					else
					{
						cat("chr",q,": the number of bps intersected are ",intersectionBPsum2[q]," which are ",contsBP_intersectedPercent2[q],"% of all the bps of contact bands in the chromosome\n",sep="")				
					}
				}
				cat("\n")

				#perecentage of all bps of FPKM fragments (RNA-seq)
				cat("\nexpression file - ",FPKM_filename1,":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis)
					{
						cat("chr",q," (cis): the number of bps intersected are ",intersectionBPsum1[q]," which are ",fpkmBP_intersectedPercent1[q],"% of all the bps of FPKMs above the cutoff in the chromosome\n",sep="")							
					}
					else
					{
						cat("chr",q,": the number of bps intersected are ",intersectionBPsum1[q]," which are ",fpkmBP_intersectedPercent1[q],"% of all the bps of FPKMs above the cutoff in the chromosome\n",sep="")				
					}
				}
				cat("\n")
				
				cat("\nexpression file - ",FPKM_filename2,":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis)
					{
						cat("chr",q," (cis): the number of bps intersected are ",intersectionBPsum2[q]," which are ",fpkmBP_intersectedPercent2[q],"% of all the bps of FPKMs above the cutoff in the chromosome\n",sep="")							
					}
					else
					{
						cat("chr",q,": the number of bps intersected are ",intersectionBPsum2[q]," which are ",fpkmBP_intersectedPercent2[q],"% of all the bps of FPKMs above the cutoff in the chromosome\n",sep="")				
					}
				}
				cat("\n")
				
				#percentage of sum of FPKMs
				cat("\nexpression file - ",FPKM_filename1,":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis)
					{
						cat("chr",q," (cis): the number of FPKMs intersected are ",intersectionFPKMsum1[q]," which are ",fpkm_intersectedPercent1[q],"% of all the FPKMs above the cutoff in the chromosome\n",sep="")								
					}
					else
					{
						cat("chr",q,": the number of FPKMs intersected are ",intersectionFPKMsum1[q]," which are ",fpkm_intersectedPercent1[q],"% of all the FPKMs above the cutoff in the chromosome\n",sep="")				
					}	
				}
				cat("\n")	

				cat("\nexpression file - ",FPKM_filename2,":\n")
				for(q in 1:numOFchroms)
				{
					if(q == cis)
					{
						cat("chr",q," (cis): the number of FPKMs intersected are ",intersectionFPKMsum2[q]," which are ",fpkm_intersectedPercent2[q],"% of all the FPKMs above the cutoff in the chromosome\n",sep="")								
					}
					else
					{
						cat("chr",q,": the number of FPKMs intersected are ",intersectionFPKMsum2[q]," which are ",fpkm_intersectedPercent2[q],"% of all the FPKMs above the cutoff in the chromosome\n",sep="")				
					}	
				}
				cat("\n")

				#creating plot		
				barName1 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",FPKM_filename1,"on the plot:\n"))
				barName2 <- readline(prompt=cat("\nenter the name that you would like to give the bars that represent the data from",FPKM_filename2,"on the plot:\n"))
				expression_source <- c(replicate(numOFchroms,barName1),replicate(numOFchroms,barName2))
				chromosomes <- c(chr_names,chr_names)
				conts_base_pairs_percentage <- c(contsBP_intersectedPercent1,contsBP_intersectedPercent2)
				FPKM_base_pairs_percentage <- c(fpkmBP_intersectedPercent1,fpkmBP_intersectedPercent2)
				FPKM_percentage <- c(fpkm_intersectedPercent1,fpkm_intersectedPercent2)
				
				conts_bp_plotData <- data.frame(expression_source,chromosomes,conts_base_pairs_percentage)
				fpkm_bp_plotData <- data.frame(expression_source,chromosomes,FPKM_base_pairs_percentage)
				fpkm_plotData <- data.frame(expression_source,chromosomes,FPKM_percentage)
				
				#creating the percentage of bp intersected from contact bands plot:
				
				#cat("\nplotting the contact bands bp percentage data\nplease wait...\n\n")
				#print(ggplot2::ggplot(conts_bp_plotData, ggplot2::aes(expression_source,conts_base_pairs_percentage,fill=chromosomes)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected bp from contact bands between expression and contacts of",sp1[1],sp1[2],sp1[3],"for each chromosome separately")) + ggplot2::labs(x="expression source",y="intersected bp percentage(%)")) + ggplot2::theme(plot.title = ggplot2::element_text(size=13))

				plo_contsBP <- ggplot2::ggplot(conts_bp_plotData, ggplot2::aes(chromosomes,conts_base_pairs_percentage,fill=expression_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected bp from contact bands between expression and contacts of",sp1[1],sp1[2],sp1[3],"for each chromosome separately")) + ggplot2::labs(x="expression source",y="intersected bp percentage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				cat("\nperforming prop tests\nthe 1st is the comparison of proportions of intersections out of all bps in the contact bands\nthe 2nd creates a table of all the successes (bp intersections) and fails (all the bp of the contact bands that didn't intersect) and uses them\n\n")
				for(r in 1:length(unique(chromosomes)))
				{
					#prop_contsBP <- prop.test(conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],]$conts_base_pairs_percentage,c(all_contsBPsum[r],all_contsBPsum[r]),correct=FALSE)
					if(intersectionBPsum1[r]>0 & intersectionBPsum2[r]>0)
					{
						cat("\nchromosome ",r,":\n",sep="")
						prop_contsBP1 <- prop.test(c(intersectionBPsum1[r],intersectionBPsum2[r]),c(all_contsBPsum[r],all_contsBPsum[r]),correct=FALSE)
						print(prop_contsBP1)
						prop_contsBP2 <- prop.test(matrix(c(all_contsBPsum[r]-intersectionBPsum1[r],intersectionBPsum1[r],all_contsBPsum[r]-intersectionBPsum2[r],intersectionBPsum2[r]),ncol=2),correct=FALSE)
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
						compute.es::propes((intersectionBPsum1[r]/all_contsBPsum[r]),(intersectionBPsum2[r]/all_contsBPsum[r]),all_contsBPsum[r],all_contsBPsum[r])
						cat("\n")
						
						#phi
						cat("\nphi:\n")
						print(psych::phi(matrix(c(intersectionBPsum1[r],all_contsBPsum[r]-intersectionBPsum1[r],intersectionBPsum2[r],all_contsBPsum[r]-intersectionBPsum2[r]),ncol=2,byrow=TRUE)))
						cat("\n")
						
						#creating a barplot with significance stars	
						if(prop_contsBP$p.value <= 0.05)
						{
							exps <- conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],]$conts_base_pairs_percentage),]$expression_source		
							chroms <- conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],]$conts_base_pairs_percentage),]$chromosomes
							perces <- conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(conts_bp_plotData[conts_bp_plotData$chromosomes==unique(chromosomes)[r],]$conts_base_pairs_percentage),]$conts_base_pairs_percentage+1
							label.df <- data.frame(expression_source=exps,chromosomes=chroms,conts_base_pairs_percentage=perces)

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
					nm <- paste("contactsBP_FPKMintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file"
					
					bpORfpkm <- "bp"
					
					percentageOf <- "Contact bands"

					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}
					
					Intersection_chromosomes <- "each chromosome separately"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(nm,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
				}
		
				#creating the percentage of bp from RNA-seq bands plot:
				
				#cat("\nplotting the FPKM bands bp percentage data\nplease wait...\n\n")
				#print(ggplot2::ggplot(fpkm_bp_plotData, ggplot2::aes(expression_source,FPKM_base_pairs_percentage,fill=chromosomes)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected bp from RNA-seq bands between expression and contacts of",sp1[1],sp1[2],sp1[3],"for each chromosome separately")) + ggplot2::labs(x="expression source",y="intersected bp percentage(%)")) + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				
				plo_fpkmBP <- ggplot2::ggplot(fpkm_bp_plotData, ggplot2::aes(chromosomes,FPKM_base_pairs_percentage,fill=expression_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected bp from RNA-seq bands between expression and contacts of",sp1[1],sp1[2],sp1[3],"for each chromosome separately")) + ggplot2::labs(x="expression source",y="intersected bp percentage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				cat("\nthe 1st is the comparison of proportions of intersections out of all bps in the FPKM bands\nthe 2nd creates a table of all the successes (bp of FPKMs intersections) and fails (all the bp of the FPKM bands that didn't intersect) and uses them\n\n")
				for(r in 1:length(unique(chromosomes)))
				{
					#prop_fpkmBP <- prop.test(fpkm_bp_plotData[fpkm_bp_plotData$chromosomes==unique(chromosomes)[r],]$FPKM_base_pairs_percentage,c(all_fpkmBPsum1[r],all_fpkmBPsum2[r]),correct=FALSE)
					if(intersectionBPsum1[r]>0 & intersectionBPsum2[r]>0)
					{
						cat("\nchromosome ",r,":\n",sep="")
						prop_fpkmBP1 <- prop.test(c(intersectionBPsum1[r],intersectionBPsum2[r]),c(all_fpkmBPsum1[r],all_fpkmBPsum2[r]),correct=FALSE)
						print(prop_fpkmBP1)
						prop_fpkmBP2 <- prop.test(matrix(c(all_fpkmBPsum1[r]-intersectionBPsum1[r],intersectionBPsum1[r],all_fpkmBPsum2[r]-intersectionBPsum2[r],intersectionBPsum2[r]),ncol=2),correct=FALSE)
						print(prop_fpkmBP2)
						prop_ans2 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
						if(prop_ans2 == 1)
						{
							prop_fpkmBP <- prop_fpkmBP1
						}
						else
						{
							prop_fpkmBP <- prop_fpkmBP2
						}
						
						#calculating effect sizes
						cat("\neffect sizes:\n")
						compute.es::propes((intersectionBPsum1[r]/all_fpkmBPsum1[r]),(intersectionBPsum2[r]/all_fpkmBPsum2[r]),all_fpkmBPsum1[r],all_fpkmBPsum2[r])
						cat("\n")
						
						#phi
						cat("\nphi:\n")
						print(psych::phi(matrix(c(intersectionBPsum1[r],all_fpkmBPsum1[r]-intersectionBPsum1[r],intersectionBPsum2[r],all_fpkmBPsum2[r]-intersectionBPsum2[r]),ncol=2,byrow=TRUE)))
						cat("\n")
						
						#creating a barplot with significance stars	
						if(prop_fpkmBP$p.value <= 0.05)
						{
							exps <- fpkm_bp_plotData[fpkm_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(fpkm_bp_plotData[fpkm_bp_plotData$chromosomes==unique(chromosomes)[r],]$FPKM_base_pairs_percentage),]$expression_source		
							chroms <- fpkm_bp_plotData[fpkm_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(fpkm_bp_plotData[fpkm_bp_plotData$chromosomes==unique(chromosomes)[r],]$FPKM_base_pairs_percentage),]$chromosomes
							perces <- fpkm_bp_plotData[fpkm_bp_plotData$chromosomes==unique(chromosomes)[r],][which.max(fpkm_bp_plotData[fpkm_bp_plotData$chromosomes==unique(chromosomes)[r],]$FPKM_base_pairs_percentage),]$FPKM_base_pairs_percentage+1
							label.df <- data.frame(expression_source=exps,chromosomes=chroms,FPKM_base_pairs_percentage=perces)

							if(prop_fpkmBP$p.value <= 0.0001)
							{
								plo_fpkmBP <- plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "****",size=8)
							}
							else if(prop_fpkmBP$p.value <= 0.001)
							{
								plo_fpkmBP <- plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "***",size=8)
							}
							else if(prop_fpkmBP$p.value <= 0.01)
							{
								plo_fpkmBP <- plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "**",size=8)
							}
							else if(prop_fpkmBP$p.value <= 0.05)
							{
								plo_fpkmBP <- plo_fpkmBP + ggplot2::geom_text(data = label.df, label = "*",size=8)
							}
						}
					}	
				}
				cat("\nplotting the FPKM bands bp percentage data\nplease wait...\n\n")
				print(plo_fpkmBP)
							
				ans5 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans5 == "y")
				{
					nm <- paste("fpkmBP_FPKMintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file"
					
					bpORfpkm <- "bp"
					
					percentageOf <- "RNA-seq fragments"

					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}
					
					Intersection_chromosomes <- "each chromosome separately"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(nm,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)				
				}
		
				#creating the percentage of FPKM intersected of sum of FPKM plot:
				
				#cat("\nplotting the FPKM percentage data\nplease wait...\n\n")
				#print(ggplot2::ggplot(fpkm_plotData, ggplot2::aes(expression_source,FPKM_percentage,fill=chromosomes)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected FPKMs between expression and contacts of",sp1[1],sp1[2],sp1[3],"for each chromosome separately")) + ggplot2::labs(x="expression source",y="intersected FPKM percentage(%)")) + ggplot2::theme(plot.title = ggplot2::element_text(size=13))

				ploFPKM <- ggplot2::ggplot(fpkm_plotData, ggplot2::aes(chromosomes,FPKM_percentage,fill=expression_source)) + ggplot2::geom_bar(stat="identity", position = "dodge") + ggplot2::ggtitle(paste("percentage of intersected FPKMs between expression and contacts of",sp1[1],sp1[2],sp1[3],"for each chromosome separately")) + ggplot2::labs(x="expression source",y="intersected FPKM percentage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=13))
				cat("\nthe 1st is the comparison of proportions of intersections out of all FPKMs in the FPKM bands\nthe 2nd creates a table of all the successes (FPKM intersections) and fails (all the FPKM of the FPKM bands that didn't intersect) and uses them\n\n")					
				for(r in 1:length(unique(chromosomes)))
				{
					if(intersectionFPKMsum1[r]>0 & intersectionFPKMsum2[r]>0)
					{
						cat("\nchromosome ",r,":\n",sep="")
						#propFPKM <- prop.test(fpkm_plotData[fpkm_plotData$chromosomes==unique(chromosomes)[r],]$FPKM_percentage,c(all_fpkmSum1[r],all_fpkmSum2[r]),correct=FALSE)
						propFPKM1 <- prop.test(c(intersectionFPKMsum1[r],intersectionFPKMsum2[r]),c(all_fpkmSum1[r],all_fpkmSum2[r]),correct=FALSE)
						print(propFPKM1)
						propFPKM2 <- prop.test(matrix(c(all_fpkmSum1[r]-intersectionFPKMsum1[r],intersectionFPKMsum1[r],all_fpkmSum2[r]-intersectionFPKMsum2[r],intersectionFPKMsum2[r]),ncol=2),correct=FALSE)
						print(propFPKM2)					
						prop_ans3 <- as.integer(readline(prompt=cat("\nwhich of the results would you like to use?\n1) first\n2) second\n\n")))
						if(prop_ans3 == 1)
						{
							propFPKM <- propFPKM1
						}
						else
						{
							propFPKM <- propFPKM2
						}
						
						#calculating effect sizes
						cat("\neffect sizes:\n")
						compute.es::propes((intersectionFPKMsum1[r]/all_fpkmSum1[r]),(intersectionFPKMsum2[r]/all_fpkmSum2[r]),all_fpkmSum1[r],all_fpkmSum2[r])
						cat("\n")
						
						#phi
						cat("\nphi:\n")
						print(psych::phi(matrix(c(intersectionFPKMsum1[r],all_fpkmSum1[r]-intersectionFPKMsum1[r],intersectionFPKMsum2[r],all_fpkmSum2[r]-intersectionFPKMsum2[r]),ncol=2,byrow=TRUE)))
						cat("\n")
						
						#creating a barplot with significance stars	
						if(propFPKM$p.value <= 0.05)
						{
							exps <- fpkm_plotData[fpkm_plotData$chromosomes==unique(chromosomes)[r],][which.max(fpkm_plotData[fpkm_plotData$chromosomes==unique(chromosomes)[r],]$FPKM_percentage),]$expression_source		
							chroms <- fpkm_plotData[fpkm_plotData$chromosomes==unique(chromosomes)[r],][which.max(fpkm_plotData[fpkm_plotData$chromosomes==unique(chromosomes)[r],]$FPKM_percentage),]$chromosomes
							perces <- fpkm_plotData[fpkm_plotData$chromosomes==unique(chromosomes)[r],][which.max(fpkm_plotData[fpkm_plotData$chromosomes==unique(chromosomes)[r],]$FPKM_percentage),]$FPKM_percentage+1
							label.df <- data.frame(expression_source=exps,chromosomes=chroms,FPKM_percentage=perces)

							if(propFPKM$p.value <= 0.0001)
							{
								ploFPKM <- ploFPKM + ggplot2::geom_text(data = label.df, label = "****",size=8)
							}
							else if(propFPKM$p.value <= 0.001)
							{
								ploFPKM <- ploFPKM + ggplot2::geom_text(data = label.df, label = "***",size=8)
							}
							else if(propFPKM$p.value <= 0.01)
							{
								ploFPKM <- ploFPKM + ggplot2::geom_text(data = label.df, label = "**",size=8)
							}
							else if(propFPKM$p.value <= 0.05)
							{
								ploFPKM <- ploFPKM + ggplot2::geom_text(data = label.df, label = "*",size=8)
							}
						}
					}	
				}
				cat("\nplotting the FPKM percentage data\nplease wait...\n\n")
				print(ploFPKM)
				
				ans6 <- readline(prompt=cat("\n\nsave the plot?\ny/n\n"))
				if(ans6 == "y")
				{
					#getting the date and time in order to distinguish between file names of plots
					nm <- paste("fpkm_FPKMintersectedContacts_",DandT2,".png",sep="")
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file"
					
					bpORfpkm <- "fpkm"
					
					percentageOf <- "FPKM"
					
					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "each chromosome separately"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(nm,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}
			}

			#venn diagrams
			if(readline(prompt=cat("\nwould you like to create venn diagrams?\ny/n\n\n")) == "y")
			{
				#creating the venn data frame for 0 percent
				venn_DF <- data.frame(rep(1,nrow(conts)),matrix(0,nrow(conts),2))##!

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
				
				venn_flag <- 0 #if venn_flag is 1 then we save a diagram and it should be recorded in 'expressionVScontacts_sumOFintersections_plots'
				#filling in the data for the venn diagrams
				if(venn_ans1 == "y")
				{
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM1_afterCO.bed -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM2_afterCO.bed -f ",ovlp_cov," -c > ~/Analyze4C/temp/venn_A2.bed",sep=""))
					venn_A1 <- read.table(paste("~/Analyze4C/temp/venn_A1.bed",sep=""))
					venn_A2 <- read.table(paste("~/Analyze4C/temp/venn_A2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn_A1.bed")
					system("rm ~/Analyze4C/temp/venn_A2.bed")
					if(venn_ans2 == "y")
					{
						system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM1_afterCO.bed -c > ~/Analyze4C/temp/venn_B1.bed",sep=""))
						system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM2_afterCO.bed -c > ~/Analyze4C/temp/venn_B2.bed",sep=""))
						venn_B1 <- read.table(paste("~/Analyze4C/temp/venn_B1.bed",sep=""))
						venn_B2 <- read.table(paste("~/Analyze4C/temp/venn_B2.bed",sep=""))
						#removing files
						system("rm ~/Analyze4C/temp/venn_B1.bed")
						system("rm ~/Analyze4C/temp/venn_B2.bed")
						#if the number in the 1st is 0 but the number in the 2nd is more than ovlp_num then we will consider it to be 1
						venn_A1[venn_B1[,4]>=ovlp_num,4] <- 1
						venn_A2[venn_B2[,4]>=ovlp_num,4] <- 1
					}
					#add the correct number to each column and create venn diagram using the data frame venn_DF
					for(w in 1:numOFchroms)
					{
						cat("\nchromosome ",w,":\n\n",sep="")
						venn_DF[venn_A1[,4]>=1 & venn_A1[,1]==w,2] <- 1
						venn_DF[venn_A2[,4]>=1 & venn_A2[,1]==w,3] <- 1
						vennName <- paste("venn_psContsVS2FPKMs_chr",w,"_",DandT2,".jpg",sep="")
						venn_flag <- venn_creator(venn_DF,2,1,vennName)						
					}	
				}
				else
				{
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM1_afterCO.bed -c > ~/Analyze4C/temp/venn1.bed",sep=""))
					system(paste("bedtools intersect -a ~/Analyze4C/contact_bands/",conts_filename," -b ~/Analyze4C/temp/FPKM2_afterCO.bed -c > ~/Analyze4C/temp/venn2.bed",sep=""))
					venn1 <- read.table(paste("~/Analyze4C/temp/venn1.bed",sep=""))
					venn2 <- read.table(paste("~/Analyze4C/temp/venn2.bed",sep=""))
					system("rm ~/Analyze4C/temp/venn1.bed")
					system("rm ~/Analyze4C/temp/venn2.bed")
					#add the correct number to each column and create venn diagram using the data frame venn_DF
					for(w in 1:numOFchroms)
					{
						cat("\nchromosome ",w,":\n\n",sep="")
						venn_DF[venn1[,4]>=1 & venn1[,1]==w,2] <- 1
						venn_DF[venn2[,4]>=1 & venn2[,1]==w,3] <- 1
						vennName <- paste("venn_psContsVS2FPKMs_chr",w,"_",DandT2,".jpg",sep="")
						venn_flag <- venn_creator(venn_DF,2,1,vennName)	
					}	
				}
								
				#recording the data into 'expressionVScontacts_sumOFintersections_plots.txt'
				if(venn_flag == 1)
				{					
					#getting parameters for expressionVScontacts_plots.txt:
					struct <- "two expression files intersected with a contact bands file - Venn Diagram"					
					bpORfpkm <- NA
					percentageOf <- NA
					if(ans4 == 1)
					{
						FPKM_CO_appliedTo <- "whole genome"
					}
					else if(ans4 == 2)
					{
						FPKM_CO_appliedTo <- "trans"
					}
					else if(ans4 == 3)
					{
						FPKM_CO_appliedTo <- "each chromosome separately"
					}

					Intersection_chromosomes <- "each chromosome separately"
					
					ans7 <- readline(prompt=cat("\nwould you like to include any additional notes about the plot (this all will be written in the expressionVScontacts_plots.txt file)?\ny/n\n"))
					if(ans7 == "y")
					{
						notes <- readline(prompt=cat("\nenter the notes:\n\n"))
					}
					else
					{
						notes <- NA
					}
					
					#saving the parameters and details to expressionVScontacts_plots.txt
					expressionVScontacts_plots[nrow(expressionVScontacts_plots)+1,] <- c(vennName,struct,conts_filename,FPKM_filename1,FPKM_filename2,bpORfpkm,percentageOf,DandT1,CO_type,FPKM_cutoff,FPKM_CO_appliedTo,Intersection_chromosomes,notes)
					#sorting the list of experiments by bait alphabetically (and sorting the row indices)
					expressionVScontacts_plots <- expressionVScontacts_plots[order(expressionVScontacts_plots$Plotfile_name),]
					rownames(expressionVScontacts_plots) <- seq(length=nrow(expressionVScontacts_plots))
					#adding the new data to the file (by erasing the old one and creating a new one)
					system("rm expressionVScontacts_plots.txt")
					write.table(expressionVScontacts_plots,"expressionVScontacts_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
				}	
			}							
		}		
		else if(ans3 == 4) #exit
		{
			#removing the files in temp
			system("rm ~/Analyze4C/temp/FPKM1_afterCO.bed")
			system("rm ~/Analyze4C/temp/FPKM2_afterCO.bed")
			system("rm ~/Analyze4C/temp/expressionComparison1.bed")
			system("rm ~/Analyze4C/temp/expressionComparison2.bed")
			break
		}				
	}
}
