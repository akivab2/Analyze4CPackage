#' @export

#the main purpose of the function is to take a raw data file, remove different amounts of coverage and see how that affects the similarities between this data and any other files data (after calculating contact bands)
#the measurements are precision, recall, and f-measure
#the order of things is first getting a raw data file, then calculating p-scores, then creating contacts. Then intersecting the contacts with itself. Each step here is done with parameters decided by user.
#after this, coverage is removed from the raw data, p-score, and contact band making.
#The contacts are then intersected with the static raw datas contact bands. This will give use precision, recall, and f-measures.
#these measures how different or the same the static contact bands are from the new ones (after coverage removal)
#The user can choose how much coverage to remove and how many times.
#These options include limiting the removal by coverage percent, limiting by minimum f-measure reached, or limiting by actual minimum coverage
#the user chooses a step for each (or number of steps in actual minimum coverage limiting)
#the user also chooses from what chromosomes the coverage is removed, and for which the precision, recall, and f-measure are calculated
#the data is then plotted. x axis is the coverage or coverage removal percentages, the y axis is for the precision, recall and f-measure.
#Each one of these are plotted on a separate line, to show their change as the coverage changes
#if a plot is saved, its details are saved in coverageVSfmeasure_plots.txt

coverage_change_tester2 <- function(Experiments_4C,coverageVSfmeasure_plots,rearranged_rawData)
{
#description of the function:
{
	cat("\ncoverage/coverage removal comparison to F-measure with self:\n\n")
	cat("this function tests in multiple ways how the coverage change changes the contacts of any single raw data file\n")
	cat("the way this is done is by reducing coverages for a specific raw data file, getting contacts from before and after coverage removal, and calculating the f-measure.\n")
	cat("this data is plotted and we are able to see how the coverage change influences the f-measure\n\n")
}

#checking that the 'temp' folder is empty, if not then empty it out
{
	ls_sgr <- system("ls ~/Analyze4C/temp/",intern=TRUE)
	if(identical(grep(".sgr",ls_sgr),integer(0)) == "FALSE")
	{
		system("rm ~/Analyze4C/temp/*.sgr")
	}
	
	ls_bed <- system("ls ~/Analyze4C/temp/",intern=TRUE)
	if(identical(grep(".bed",ls_bed),integer(0)) == "FALSE")
	{
		system("rm ~/Analyze4C/temp/*.bed")
	}
}

#choosing the 1st raw data file to work with  (the one which will have coverage removed from):
{
	#getting the name of the raw data file the user wants to use
	repeat
	{
		choice2_rmv <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the 1st raw data file (the one that will have coverage removed from it) from:\n\n1) original\n2) rearranged\n3) coverage removed\n\n")))
		if(choice2_rmv == 2)
		{
			ls_files_rmv <- system("ls ~/Analyze4C/rawData/rearranged",intern=TRUE)
		}
		else if(choice2_rmv == 3)
		{
			ls_files_rmv <- system("ls ~/Analyze4C/rawData/coverage_removed",intern=TRUE)
		}
		else
		{
			ls_files_rmv <- system("ls ~/Analyze4C/rawData/original",intern=TRUE)
		}		
		file.name_rmv <- readline(prompt=cat("These are the raw data files available\nwhich would you like to use?\n",ls_files_rmv,"\n",sep="\n"))
		ind_files_rmv <- pmatch(file.name_rmv,ls_files_rmv,nomatch=0)
		if(ind_files_rmv == 0)
		{
			cat("no such file exists.\nplease try again.\n\n")
		}
		else
		{
			break
		}
	}
	#importing the data from the file
	if(choice2_rmv == 2)
	{
		rawData_rmv <- read.table(paste("~/Analyze4C/rawData/rearranged/",file.name_rmv,sep=""))
	}
	else if(choice2_rmv == 3)
	{
		rawData_rmv <- read.table(paste("~/Analyze4C/rawData/coverage_removed/",file.name_rmv,sep=""))
	}
	else
	{
		rawData_rmv <- read.table(paste("~/Analyze4C/rawData/original/",file.name_rmv,sep=""))
	}
	
	#finding the specific raw data files information in Experiments_4C
	sp1_rmv <- unlist(strsplit(file.name_rmv,"[.]"))
	sp2_rmv <- strsplit(sp1_rmv[[1]][1],"_")
	out_rmv <- findIn.Experiments_4C(sp2_rmv[[1]][1],sp2_rmv[[1]][2],sp2_rmv[[1]][3],Experiments_4C)
	
	#getting the cis_rmv chromosome number
	cis_rmv <- as.numeric(out_rmv[2])
	
	#getting the bait index
	vp.pos_rmv <- as.numeric(out_rmv[3])
	
	#getting the amount of bp around the bait that will be erased
	erase_rmv <- 1e6 #for now this will be the default. later i could have the user asked how many bp they want to be removed
}

#choosing the 2nd raw data file to work with:
{
	#getting the name of the raw data file the user wants to use
	repeat
	{
		choice2_static <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the 2nd raw data file from:\n\n1) original\n2) rearranged\n3) coverage removed\n\n")))
		if(choice2_static == 2)
		{
			ls_files_static <- system("ls ~/Analyze4C/rawData/rearranged",intern=TRUE)
		}
		else if(choice2_static == 3)
		{
			ls_files_static <- system("ls ~/Analyze4C/rawData/coverage_removed",intern=TRUE)
		}
		else
		{
			ls_files_static <- system("ls ~/Analyze4C/rawData/original",intern=TRUE)
		}		
		file.name_static <- readline(prompt=cat("These are the raw data files available\nwhich would you like to use?\n",ls_files_static,"\n",sep="\n"))
		ind_files_static <- pmatch(file.name_static,ls_files_static,nomatch=0)
		if(ind_files_static == 0)
		{
			cat("no such file exists.\nplease try again.\n\n")
		}
		else
		{
			break
		}
	}
	#importing the data from the file
	if(choice2_static == 2)
	{
		rawData_static <- read.table(paste("~/Analyze4C/rawData/rearranged/",file.name_static,sep=""))
	}
	else if(choice2_static == 3)
	{
		rawData_static <- read.table(paste("~/Analyze4C/rawData/coverage_removed/",file.name_static,sep=""))
	}
	else
	{
		rawData_static <- read.table(paste("~/Analyze4C/rawData/original/",file.name_static,sep=""))
	}
	
	#finding the specific raw data files information in Experiments_4C
	sp1_static <- unlist(strsplit(file.name_static,"[.]"))
	sp2_static <- strsplit(sp1_static[[1]][1],"_")
	out_static <- findIn.Experiments_4C(sp2_static[[1]][1],sp2_static[[1]][2],sp2_static[[1]][3],Experiments_4C)
	
	#getting the cis_static chromosome number
	cis_static <- as.numeric(out_static[2])
	
	#getting the bait index
	vp.pos_static <- as.numeric(out_static[3])
	
	#getting the amount of bp around the bait that will be erased
	erase_static <- 1e6 #for now this will be the default. later i could have the user asked how many bp they want to be removed
}

#choosing a method for function to use and getting parameters:	
{
	#getting the current coverage (either of trans or all)
	cat("\ndata stats for ",file.name_rmv,":\n",sep="")
	rawData.stats(Experiments_4C,rawData_rmv,file.name_rmv)

	cat("\ndata stats for ",file.name_static,":\n",sep="")
	rawData.stats(Experiments_4C,rawData_static,file.name_static)

	choice0 <- as.integer(readline(prompt=cat("\nwhat method would you like to use?:\n\n1) choosing a limit of coverage percentage removal\n2) choosing an f-measure limit\n3) choosing of a minimum coverage limit\n4) removing minimum number of reads until a choice of specific number of reads\n\n")))
	if(choice0 == 1)
	{
		covPer_lim <- as.numeric(readline(prompt=cat("\n\nenter the most percentage (the limit) of coverage you would like removed from the original raw data (between 0 and 100):\n\n")))
		stp <- as.numeric(readline(prompt=cat("\nenter the step of coverage removal (between 0 and 100):\n\n"))) #getting the step size of coverage removal
		
		#getting the number of cycles of the same coverage removal percentage
		cyc <- as.integer(readline(prompt=cat("\nthe same coverage removal percentage can be done multiple times. The f-measure that will eventually be used is the median these.\nHow many times would you like the same coverage removal percentage be tested?\n\n")))				
	}
	else if(choice0 == 2)
	{
		fmeasure_lim <- as.numeric(readline(prompt=cat("\nenter the lowest f-measure you would like to be reached:\n\n")))
		stp <- as.numeric(readline(prompt=cat("\nenter the step of coverage removal (between 0 and 100):\n\n"))) #getting the step size of coverage removal

		#getting the number of cycles of the same coverage removal percentage
		cyc <- as.integer(readline(prompt=cat("\nthe same coverage removal percentage can be done multiple times. The f-measure that will eventually be used is the median these.\nHow many times would you like the same coverage removal percentage be tested?\n\n")))		
	}
	else if(choice0 == 3)
	{
		cov_lim <- as.numeric(readline(prompt=cat("\nenter the lowest coverage you would like be reached (between 0 and 1):\n\n")))
		numberOfSteps <- as.integer(readline(prompt=cat("\nin how many steps would you like to get there?\n\n")))
		
		#getting the number of cycles of the same coverage removal percentage
		cyc <- as.integer(readline(prompt=cat("\nthe same coverage removal percentage can be done multiple times. The f-measure that will eventually be used is the median these.\nHow many times would you like the same coverage removal percentage be tested?\n\n")))		
	}
	else if(choice0 == 4)
	{
		first_min <- as.numeric(readline(prompt=cat("\nenter the first number of minimum RE sites you would like be tested:\n\n")))
		last_min_reads <- as.numeric(readline(prompt=cat("\nenter the last number of minimum RE sites you would like be tested:\n\n")))
		min_read_stp <- as.integer(readline(prompt=cat("\nhow many RE sites should be added per step?\n\n")))
	}
			
	#asking for what chromosomes will there be coverage removal
	choice1 <- as.integer(readline(prompt=cat("\nwhat chromosomes should the coverage removal be applied to? (enter the option number)\n1) all chromosomes\n2) trans\n3) cis\n4) specific chromosome\n\n")))
	if(choice1 == 1) #getting the starting coverage of the whole genome
	{
		rem_chr <- 0 #if no chromosome is chosen, rem_chr will get 0 for coverageChanger_chrChooser
		cur_cov <- sum(rawData_rmv[,3]>0)/sum(rawData_rmv[,3]>=0)
		cat("\nthe starting coverage of chromosome the whole genome is",cur_cov,"\n\n")		
	}
	else if(choice1 == 2) #getting the starting coverage for trans
	{
		rem_chr <- 0 #if no chromosome is chosen, rem_chr will get 0 for coverageChanger_chrChooser
		cur_cov <- sum(rawData_rmv[rawData_rmv[,1]!=cis_rmv,3]>0)/sum(rawData_rmv[rawData_rmv[,1]!=cis_rmv,3]>=0)
		cat("\nthe starting coverage of chromosome trans is",cur_cov,"\n\n")		
	}
	else if(choice1 == 3) #getting the starting coverage for cis
	{
		rem_chr <- 0 #if no chromosome is chosen, rem_chr will get 0 for coverageChanger_chrChooser
		cur_cov <- sum(rawData_rmv[rawData_rmv[,1]==cis_rmv,3]>0)/sum(rawData_rmv[rawData_rmv[,1]==cis_rmv,3]>=0)
		cat("\nthe starting coverage of chromosome cis is",cur_cov,"\n\n")		
	}	
	else if(choice1 == 4) #getting the chromosome the user would like coverage removal from, and the starting coverage of the chromosome
	{
		rem_chr <- as.integer(readline(prompt=cat("\nchoose the chromosome you would like to remove coverage from:\n\n")))
		cur_cov <- sum(rawData_rmv[rawData_rmv[,1]!=rem_chr,3]>0)/sum(rawData_rmv[rawData_rmv[,1]!=rem_chr,3]>=0)
		cat("\nthe starting coverage of chromosome",rem_chr,"is",cur_cov,"\n\n")
	}
}
		
#getting all the parameters that will be used for all the coverage removals and contact band making:
{	
	#p-score parameters:
	cat("\n*******************************************************************\n\n")
	cat("\n                       p-score parameters:\n\n")
	cat("\nfor file ",file.name_rmv,":\n",sep="")
	#get the window size in RE site
	REs_rmv <- as.integer(readline(prompt=cat("\nplease enter the size of window by RE sites (should be an even number):\n\n")))
	
	#getting the bp window size if not a rearranged data
	if(choice2_rmv != 2)
	{
		#get the window size in bp
		wind_rmv <- as.integer(readline(prompt=cat("\nplease enter the size of window by bp (note that the size should probably be proportionate to the RE site window size):\n\n")))
	}
		
	cat("\nfor file ",file.name_static,":\n",sep="")
	#get the window size in RE site
	REs_static <- as.integer(readline(prompt=cat("\nplease enter the size of window by RE sites (should be an even number):\n\n")))
	
	#getting the bp window size if not a rearranged data
	if(choice2_static != 2)
	{
		#get the window size in bp
		wind_static <- as.integer(readline(prompt=cat("\nplease enter the size of window by bp (note that the size should probably be proportionate to the RE site window size):\n\n")))
	}

	#contact band parameters:
	cat("\n*******************************************************************\n\n")
	cat("\n                      conact bands parameters:\n\n")
	
	ls_genomes <- system("ls ~/Analyze4C/genomes/",intern=TRUE)
	genome_file.name <- readline(prompt=cat("\nchoose the file with genome sizes that you would like to use:\n",ls_genomes,"\n",sep="\n"))
	genome_sizes <- read.table(paste("~/Analyze4C/genomes/",genome_file.name,sep=""))
	
	#getting the different paramaters and inputs from user
	cat("\nfor file ",file.name_rmv,":\n",sep="")
	RE_gap_rmv <- as.numeric(readline(prompt=cat("\nplease choose up to how many negative RE sites in between positive RE sites it is not considered a break in the contact band:\n\n")))
	bp_gap_rmv <- as.numeric(readline(prompt=cat("\nplease choose the max number of bp that are allowed to be in between RE sites and not be considered a gap:\n\n")))

	#asking user if to create contact bands using a p-score cutoff for all chromosomes or a percent cutoff
	if(choice2_rmv != 2) #if the data chosen is not rearranged
	{
		choice3_rmv <- as.integer(readline(prompt=cat("\nenter the number of what type of contact band making you would like to do:\n\n1) all chromosomes with the same cutoff (in p-score)\n2) all chromosomes with the same cutoff percentage\n\n")))
		if(choice3_rmv == 1)
		{
			co_rmv <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))
		}
		else
		{
			co_per_rmv <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff percent (between 0 and 1) that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))
		}
	}
	else #if the data chosen is rearranged, then we could only use a percentage cutoff (used later in contactBands_byChrom)
	{
		choice3_rmv <- 2
		co_per_rmv <- as.numeric(readline(prompt=cat("\nthe method of cutoff will be to apply a cutoff percentage\nplease enter a p-score cutoff percent (between 0 and 1) that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))	
	}
	
	#getting the different paramaters and inputs from user
	cat("\nfor file ",file.name_static,":\n",sep="")
	RE_gap_static <- as.numeric(readline(prompt=cat("\nplease choose up to how many negative RE sites in between positive RE sites it is not considered a break in the contact band:\n\n")))
	bp_gap_static <- as.numeric(readline(prompt=cat("\nplease choose the max number of bp that are allowed to be in between RE sites and not be considered a gap:\n\n")))
	
	#asking user if to create contact bands using a p-score cutoff for all chromosomes or a percent cutoff
	if(choice2_static != 2) #if the data chosen is not rearranged
	{
		choice3_static <- as.integer(readline(prompt=cat("\nenter the number of what type of contact band making you would like to do:\n\n1) all chromosomes with the same cutoff (in p-score)\n2) all chromosomes with the same cutoff percentage\n\n")))
		if(choice3_static == 1)
		{
			co_static <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))
		}
		else
		{
			co_per_static <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff percent (between 0 and 1) that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))
		}
	}
	else #if the data chosen is rearranged, then we could only use a percentage cutoff (used later in contactBands_byChrom)
	{
		choice3_static <- 2
		co_per_static <- as.numeric(readline(prompt=cat("\nthe method of cutoff will be to apply a cutoff percentage\nplease enter a p-score cutoff percent (between 0 and 1) that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))	
	}
		
	#f-measure parameters:
	cat("\n*******************************************************************\n\n")
	cat("\n                      f-measure parameters:\n\n")
	
	#deciding what chromosomes to calculate contact bands and f-measure for
	choice4 <- as.integer(readline(prompt=cat("\nwhat chromosomes will you like f-measure be calculated for?(choose the number of the option you want):\n\n1) all chromosomes\n2) trans\n3) cis\n4) specific chromosome\n\n")))
	if(choice4 == 4)
	{
		f_chr <- as.integer(readline(prompt=cat("\nchoose the chromosome you would like to remove coverage from:\n\n")))
	}
	
	#deciding what beta will be. The default is 1
	beta <- 1
	ans_beta <- readline(prompt=cat("\nthe default beta for the f-measure is 1.\nwould you like to change it?\ny/n\n\n"))
	if(ans_beta == "y")
	{
		beta <- -1
		while(beta <= 0)
		{
			beta <- as.numeric(readline(prompt=cat("\na beta between 0 and 1 will give more weight on precision. above 1 puts more weight on recall.\nplease enter a beta:\n\n")))
			if(beta <= 0 )
			{
				cat("\nthe beta entered is incompatible. The beta must be over 0.\n\n")
			}
		}
	}
}

#calculating p-scores for the static data (if a translocated data is chosen this must be with no kb and before translocation areas are removed):
{	
	cat("\n\ncalculating p-scores for",file.name_static,"\nplease wait...\n\n")
	if(choice2_static == 2)
	{
		#using the parameters above we calculate the p score with RE site window sizes only
		ps_static <- pScore_nokb(rawData_static,cis_static,vp.pos_static,erase_static,REs_static/2)
	}
	else
	{
		#using the parameters above we calculate the p score with RE site and bp window sizes
		ps_static <- pScore(rawData_static,cis_static,vp.pos_static,erase_static,wind_static,REs_static/2)
	}
	
	#write p score to file
	#write.table(ps_static,"~/Analyze4C/temp/pScore_static.sgr",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}

#ps_forF_static will get the p-score data that will get the contact and f-measure treatment
ps_forF_static <- ps_static #the default is all (choice4 == 1).
if(choice4 == 2) #trans only
{
	ps_forF_static <- ps_static[ps_static[,1]!=cis_static,]
}
else if(choice4 == 3) #cis only
{
	ps_forF_static <- ps_static[ps_static[,1]==cis_static,]
}
else if(choice4 == 4) #specific chromosome
{
	ps_forF_static <- ps_static[ps_static[,1]==f_chr,]
}

#applying cutoff and getting contact bands for the static data(use contactBands_byChrom_translocatedRemoval if translocation, this will first apply cutoff then remove translocations and then get contact bands):
{
	cat("\n\ncreating contact bands for the static data\nplease wait...\n\n")
	#creating the contact bands:
	if(choice2_static == 2) #if the data is rearranged
	{
		#checking if we are dealing with a file that has a translocation that isn't removed
		name_spl <- unlist(strsplit(file.name_static,"_")) #splits the raw data file name
		if(("rearranged" %in% name_spl) & !("removed" %in% name_spl)) #if the filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
		{
			conts_static <- contactBands_byChrom(ps_forF_static,RE_gap_static,bp_gap_static,co_per_static,genome_sizes,rearranged_rawData,flag=1,file.name=file.name_static) #the function 'contactBands_byChrom' with the input 1 will remove the added sections after applying the cutoff and then get contact bands
			conts_static <- conts_static[[1]]
		}
		else
		{
			conts_static <- contactBands_byChrom(ps_forF_static,RE_gap_static,bp_gap_static,co_per_static,genome_sizes,rearranged_rawData) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
			conts_static <- conts_static[[1]]
		}
	}	
	else if(choice3_static == 1)
	{
		conts_static <- contactBands(ps_forF_static,RE_gap_static,bp_gap_static,cis_static,co_static,genome_sizes)
	}
	else if(choice3_static == 2)
	{
		conts_static <- contactBands_byChrom(ps_forF_static,RE_gap_static,bp_gap_static,co_per_static,genome_sizes,rearranged_rawData)
		conts_static <- conts_static[[1]]
	}
	write.table(conts_static,paste("~/Analyze4C/temp/conts_static.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)	
}

#getting the self intersection for later use:
{
	system("bedtools intersect -a ~/Analyze4C/temp/conts_static.bed -b ~/Analyze4C/temp/conts_static.bed -wo > ~/Analyze4C/temp/static_self.bed")
}

#getting the contacts and intersections between the static data and the removeable data before removing coverage (0% cov removed)
{
	prf_all <- c() #will contain all the prf values in a way that each row is a different coverage and the columns are precision, recall, and f-measure (respectively)

	#calculating p-scores for the 1st raw data before removal of coverage (if a translocated data is chosen this must be with no kb and before translocation areas are removed):	
	cat("\n\ncalculating p-scores for",file.name_rmv,"\nplease wait...\n\n")
	if(choice2_rmv == 2)
	{
		#using the parameters above we calculate the p score with RE site window sizes only
		ps_temp <- pScore_nokb(rawData_rmv,cis_rmv,vp.pos_rmv,erase_rmv,REs_rmv/2)
	}
	else
	{
		#using the parameters above we calculate the p score with RE site and bp window sizes
		ps_temp <- pScore(rawData_rmv,cis_rmv,vp.pos_rmv,erase_rmv,wind_rmv,REs_rmv/2)
	}
	
	#ps_forF will get the p-score data that will get the contact and f-measure treatment
	ps_forF <- ps_temp #the default is all (choice4 == 1).
	if(choice4 == 2) #trans only
	{
		ps_forF <- ps_temp[ps_temp[,1]!=cis_rmv,]
	}
	else if(choice4 == 3) #cis_rmv only
	{
		ps_forF <- ps_temp[ps_temp[,1]==cis_rmv,]
	}
	else if(choice4 == 4) #specific chromosome
	{
		ps_forF <- ps_temp[ps_temp[,1]==f_chr,]
	}

	#applying cutoff and getting contact bands for 1st raw data(use contactBands_byChrom_translocatedRemoval if translocation, this will first apply cutoff then remove translocations and then get contact bands):
	if(choice4 == 1) #whole genome
	{
		cat("\n\ncreating contact bands for whole genome\nplease wait...\n\n",sep="")
	}
	else if(choice4 == 2) #trans only
	{
		cat("\n\ncreating contact bands for trans\nplease wait...\n\n",sep="")
	}
	else if(choice4 == 3) #cis_rmv only
	{
		cat("\n\ncreating contact bands for cis\nplease wait...\n\n",sep="")
	}
	else if(choice4 == 4) #specific chromosome
	{
		cat("\n\ncreating contact bands for chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
	}

	if(choice2_rmv == 2) #if the data is rearranged
	{
		#checking if we are dealing with a file that has a translocation that isn't removed
		name_spl <- unlist(strsplit(file.name_rmv,"_")) #splits the raw data file name
		if(("rearranged" %in% name_spl) & !("removed" %in% name_spl)) #if the filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
		{
			conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData,flag=1,file.name=file.name_rmv) #the function 'contactBands_byChrom' with the input 1 will remove the added sections after applying the cutoff and then get contact bands
			conts_rmv <- conts_rmv[[1]]
		}
		else
		{
			conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
			conts_rmv <- conts_rmv[[1]]
		}
	}	
	else if(choice3_rmv == 1)
	{
		conts_rmv <- contactBands(ps_forF,RE_gap_rmv,bp_gap_rmv,cis_rmv,co_rmv,genome_sizes)
	}
	else if(choice3_rmv == 2)
	{
		conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData)
		conts_rmv <- conts_rmv[[1]]
	}
	write.table(conts_rmv,paste("~/Analyze4C/temp/conts_rmv.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)	

	#getting the self intersection
	system("bedtools intersect -a ~/Analyze4C/temp/conts_rmv.bed -b ~/Analyze4C/temp/conts_rmv.bed -wo > ~/Analyze4C/temp/rmv_int_self.bed")

	#get the intersections with the static data
	system("bedtools intersect -a ~/Analyze4C/temp/conts_static.bed -b ~/Analyze4C/temp/conts_rmv.bed -wo > ~/Analyze4C/temp/rmv_int.bed")

	#calculate f-measure
	results <- system(paste("perl ~/Analyze4C/proxy/precisionRecallFmeasure.pl ~/Analyze4C/temp/static_self.bed ~/Analyze4C/temp/rmv_int_self.bed ~/Analyze4C/temp/rmv_int.bed",beta),intern=TRUE)

	#record the f-measure and the iteration
	prf <- rep(0,3) #the precision, recall, and f-measure will eventually be recorded into here (respectively)
	counter2 <- 1 #counts the index of prf
	for(z in 1:(length(results))) #iterates over the lines that were recieved into results, from the function precisionRecallFmeasure.pl
	{
		spl_results <- strsplit(results[z],":") #splits the line where there is a ":", the number result should be what comes after that (possibly with a \n will be before the number)
		if((identical(spl_results[[1]],character(0)) == "FALSE")) #testing to see that the we don't get a 'character(0)' when doing strplit, this happens when the line is just \n
		{
			#putting each result number in prf
			prf[counter2] <- as.numeric(spl_results[[1]][2])
			counter2 <-  counter2 + 1
		}
	}

	prf_all <- rbind(prf_all,prf) #adding the prf to prf_all

	#removing the contact and intersection files in case they are created multiple times
	system("rm ~/Analyze4C/temp/rmv_int_self.bed")
	system("rm ~/Analyze4C/temp/rmv_int.bed")
	system("rm ~/Analyze4C/temp/conts_rmv.bed")
}

#change coverage and print new coverage:
{	
	#cur_rmv <- stp	#"cur_rmv" is the percentage of coverage removed
	counter1 <- 1 #counts the number of iterations
	#prf_all <- c() #will contain all the prf values in a way that each row is a different coverage and the columns are precision, recall, and f-measure (respectively)
	
	if(choice0 == 1) #limited by coverage removal percentage
	{
		cur_rmv <- stp	#"cur_rmv" is the percentage of coverage removed
		while(cur_rmv <= covPer_lim)
		{	
			cat("\n*****************************************************************\n\n")
			prf_cyc <- c() #gets the prf for each coverage removal percentage, for all its cycles
			for(p in 1:cyc)
			{
				cat("\ncycle ",p," for ",cur_rmv,"% coverage removal:\n\n",sep="")
				#coverage removal
				rmv_fraction <- cur_rmv/100 #changing the percentage into fraction for input to coverageChanger_chrChooser
				out2 <- coverageChanger_chrChooser(rawData_rmv,rmv_fraction,cis_rmv,choice1,rem_chr)
				newData <- out2[[1]]
				
				#making file of removed coverage from raw data
				write.table(newData,paste("~/Analyze4C/temp/raw_",cur_rmv,"perRemoved.sgr",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
				
				#printing message according to chromosome and coverage amount removed
				if(choice1 == 1)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from the whole genome\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 2)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 3)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 4)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
				}
				
				#calculate p-scores for the rmv data
				if(choice2_rmv == 2)
				{
					#using the parameters above we calculate the p score with RE site window sizes only
					ps_temp <- pScore_nokb(newData,cis_rmv,vp.pos_rmv,erase_rmv,REs_rmv/2)
				}
				else
				{
					#using the parameters above we calculate the p score with RE site and bp window sizes
					ps_temp <- pScore(newData,cis_rmv,vp.pos_rmv,erase_rmv,wind_rmv,REs_rmv/2)
				}
				
				#write p score to file
				#write.table(ps_temp,paste("~/Analyze4C/temp/pscore_",cur_rmv,"perRemoved.sgr",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
				
				#ps_forF will get the p-score data that will get the contact and f-measure treatment
				ps_forF <- ps_temp #the default is all (choice4 == 1).
				if(choice4 == 2) #trans only
				{
					ps_forF <- ps_temp[ps_temp[,1]!=cis_rmv,]
				}
				else if(choice4 == 3) #cis_rmv only
				{
					ps_forF <- ps_temp[ps_temp[,1]==cis_rmv,]
				}
				else if(choice4 == 4) #specific chromosome
				{
					ps_forF <- ps_temp[ps_temp[,1]==f_chr,]
				}

				#printing message according to chromosome and coverage amount removed
				if(choice1 == 1) #whole genome
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis_rmv only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 2) #trans only
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 3) #cis only
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 4) #specific chromosome
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
				}
							
				#creating the contact bands:
				if(choice2_rmv == 2) #if the data is rearranged
				{
					#checking if we are dealing with a file that has a translocation that isn't removed
					name_spl <- unlist(strsplit(file.name_rmv,"_")) #splits the raw data file name
					if(("rearranged" %in% name_spl) & !("removed" %in% name_spl)) #if the filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
					{
						conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData,flag=1,file.name=file.name_rmv) #the function 'contactBands_byChrom' with the input 1 will remove the added sections after applying the cutoff and then get contact bands
						conts_rmv <- conts_rmv[[1]]
					}
					else
					{
						conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
						conts_rmv <- conts_rmv[[1]]
					}
				}
				else if(choice3_rmv == 1) #if the cutoff option chosen is by a common cutoff
				{
					conts_rmv <- contactBands(ps_forF,RE_gap_rmv,bp_gap_rmv,cis_rmv,co_rmv,genome_sizes)
				}
				else if(choice3_rmv == 2) #if the cutoff option chosen is a percentage cutoff
				{
					conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData)
					conts_rmv <- conts_rmv[[1]]
				}	
				write.table(conts_rmv,paste("~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)	
				
				#getting the intersections of contacts of removed coverage with self
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed -b ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed -wo > ~/Analyze4C/temp/",cur_rmv,"int_self.bed",sep=""))
				
				#get the intersections with the static data
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_static.bed -b ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed -wo > ~/Analyze4C/temp/",cur_rmv,"int.bed",sep=""))
				
				#calculate f-measure
				results <- system(paste("perl ~/Analyze4C/proxy/precisionRecallFmeasure.pl ~/Analyze4C/temp/static_self.bed ~/Analyze4C/temp/",cur_rmv,"int_self.bed ~/Analyze4C/temp/",cur_rmv,"int.bed ",beta,sep=""),intern=TRUE)
				
				#record the f-measure and the iteration
				prf <- rep(0,3) #the precision, recall, and f-measure will eventually be recorded into here (respectively)
				counter2 <- 1 #counts the index of prf
				for(z in 1:(length(results))) #iterates over the lines that were recieved into results, from the function precisionRecallFmeasure.pl
				{
					spl_results <- strsplit(results[z],":") #splits the line where there is a ":", the number result should be what comes after that (possibly with a \n will be before the number)
					if((identical(spl_results[[1]],character(0)) == "FALSE")) #testing to see that the we don't get a 'character(0)' when doing strplit, this happens when the line is just \n
					{
						#putting each result number in prf
						prf[counter2] <- as.numeric(spl_results[[1]][2])
						counter2 <-  counter2 + 1
					}
				}
				prf_cyc <- rbind(prf_cyc,prf)
				
				#removing the contact and intersection files in case they are created multiple times
				system(paste("rm ~/Analyze4C/temp/",cur_rmv,"int_self.bed",sep=""))
				system(paste("rm ~/Analyze4C/temp/",cur_rmv,"int.bed",sep=""))
				system(paste("rm ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed",sep=""))
			}
			
			prec_med <- median(prf_cyc[,1])
			rec_med <- median(prf_cyc[,2])
			fm_med <- median(prf_cyc[,3])
			prf_med <- c(prec_med,rec_med,fm_med)

			prf_all <- rbind(prf_all,prf_med) #adding the prf_med to prf_all
	
			#changing cur_rmv and getting the new value	
			counter1 <- counter1 + 1
			cur_rmv <- stp * counter1	
		}
		
		#creating a data frame of the recall, precision, and f-measure data
		rownames(prf_all) <- seq(1,nrow(prf_all),1)
		covs <- c(0,seq(stp,covPer_lim,stp))
		covs <- paste(covs,"%",sep="")
		len <- length(covs)
		colnames(prf_all) <- c("precision","recall","f.measure")
		prf_prec <- data.frame(rep("precision",len),subset(prf_all,select="precision"),covs)
		names(prf_prec) <- c("type","result","coverage.removed")
		prf_rec <- data.frame(rep("recall",len),subset(prf_all,select="recall"),covs)
		names(prf_rec) <- c("type","result","coverage.removed")
		prf_fm <- data.frame(rep("f.measure",len),subset(prf_all,select="f.measure"),covs)
		names(prf_fm) <- c("type","result","coverage.removed")
		prf_final <- rbind(prf_prec,prf_rec,prf_fm)
		rownames(prf_final) <- seq(1,nrow(prf_final),1)

		print(ggplot2::ggplot(data=prf_final,ggplot2::aes(x=coverage.removed,y=result,group=type)) + ggplot2::geom_line(ggplot2::aes(color=type)) + ggplot2::geom_point())		
	}
	else if(choice0 == 2)#limited by f-measure
	{
		cur_rmv <- stp	#"cur_rmv" is the percentage of coverage removed
		fmeasure <- 1000 #entering and starting off with a ridicously high fmeasure, this reasures that it will be above the fmeasure_lim so we can continue the code
		while(fmeasure >= fmeasure_lim)
		{	
			cat("\n*****************************************************************\n\n")
			prf_cyc <- c() #gets the prf for each coverage removal percentage, for all its cycles
			for(p in 1:cyc)
			{
				cat("\ncycle ",p," for ",cur_rmv,"% coverage removal:\n\n",sep="")
				#coverage removal
				rmv_fraction <- cur_rmv/100 #changing the percentage into fraction for input to coverageChanger_chrChooser
				out2 <- coverageChanger_chrChooser(rawData_rmv,rmv_fraction,cis_rmv,choice1,rem_chr)
				newData <- out2[[1]]
				
				#making file of removed coverage from raw data
				write.table(newData,paste("~/Analyze4C/temp/raw_",cur_rmv,"perRemoved.sgr",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

				#printing message according to chromosome and coverage amount removed
				if(choice1 == 1)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from the whole genome\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 2)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 3)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 4)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
				}
				
				#calculate p-scores for the rmv data
				if(choice2_rmv == 2)
				{
					#using the parameters above we calculate the p score with RE site window sizes only
					ps_temp <- pScore_nokb(newData,cis_rmv,vp.pos_rmv,erase_rmv,REs_rmv/2)
				}
				else
				{
					#using the parameters above we calculate the p score with RE site and bp window sizes
					ps_temp <- pScore(newData,cis_rmv,vp.pos_rmv,erase_rmv,wind_rmv,REs_rmv/2)
				}
				
				#write p score to file
				#write.table(ps_temp,paste("~/Analyze4C/temp/pscore_",cur_rmv,"perRemoved.sgr",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)			

				#ps_forF will get the p-score data that will get the contact and f-measure treatment
				ps_forF <- ps_temp #the default is all (choice4 == 1).
				if(choice4 == 2) #trans only
				{
					ps_forF <- ps_temp[ps_temp[,1]!=cis_rmv,]
				}
				else if(choice4 == 3) #cis only
				{
					ps_forF <- ps_temp[ps_temp[,1]==cis_rmv,]
				}
				else if(choice4 == 4) #specific chromosome
				{
					ps_forF <- ps_temp[ps_temp[,1]==f_chr,]
				}

				#printing message according to chromosome and coverage amount removed
				if(choice1 == 1) #whole genome
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 2) #trans only
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 3) #cis only
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 4) #specific chromosome
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
				}
				
				#creating the contact bands:
				if(choice2_rmv == 2) #if the data is rearranged
				{
					#checking if we are dealing with a file that has a translocation that isn't removed
					name_spl <- unlist(strsplit(file.name_rmv,"_")) #splits the raw data file name
					if(("rearranged" %in% name_spl) & !("removed" %in% name_spl)) #if the filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
					{
						conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData,flag=1,file.name=file.name_rmv) #the function 'contactBands_byChrom' with the input 1 will remove the added sections after applying the cutoff and then get contact bands
						conts_rmv <- conts_rmv[[1]]
					}
					else
					{
						conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
						conts_rmv <- conts_rmv[[1]]
					}
				}
				else if(choice3_rmv == 1) #if the cutoff option chosen is by a common cutoff
				{
					conts_rmv <- contactBands(ps_forF,RE_gap_rmv,bp_gap_rmv,cis_rmv,co_rmv,genome_sizes)
				}
				else if(choice3_rmv == 2) #if the cutoff option chosen is a percentage cutoff
				{
					conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData)
					conts_rmv <- conts_rmv[[1]]
				}	
				
				#note: i found this with ps_temp instead of ps_forF, this is probably a mistake, just keeping this here to make sure
				#creating the contact bands:
				#if(choice2_rmv == 2) #if the data is rearranged
				#{
				#	#checking if we are dealing with a file that has a translocation that isn't removed
				#	name_spl <- unlist(strsplit(file.name_rmv,"_")) #splits the raw data file name
				#	if(("rearranged" %in% name_spl) & !("removed" %in% name_spl)) #if the filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
				#	{
				#		conts_rmv <- contactBands_byChrom(ps_temp,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData,flag=1,file.name=file.name_rmv) #the function 'contactBands_byChrom' with the input 1 will remove the added sections after applying the cutoff and then get contact bands
				#		conts_rmv <- conts_rmv[[1]]
				#	}
				#	else
				#	{
				#		conts_rmv <- contactBands_byChrom(ps_temp,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
				#		conts_rmv <- conts_rmv[[1]]
				#	}
				#}
				#else if(choice3_rmv == 1) #if the cutoff option chosen is by a common cutoff
				#{
				#	conts_rmv <- contactBands(ps_temp,RE_gap_rmv,bp_gap_rmv,cis_rmv,co_rmv,genome_sizes)
				#}
				#else if(choice3_rmv == 2) #if the cutoff option chosen is a percentage cutoff
				#{
				#	conts_rmv <- contactBands_byChrom(ps_temp,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData)
				#	conts_rmv <- conts_rmv[[1]]
				#}
				write.table(conts_rmv,paste("~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)	
				
				#getting the intersections of contacts of removed coverage with self
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed -b ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed -wo > ~/Analyze4C/temp/",cur_rmv,"int_self.bed",sep=""))
				
				#get the intersections with the static data
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_static.bed -b ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed -wo > ~/Analyze4C/temp/",cur_rmv,"int.bed",sep=""))
				
				#calculate f-measure
				results <- system(paste("perl ~/Analyze4C/proxy/precisionRecallFmeasure.pl ~/Analyze4C/temp/static_self.bed ~/Analyze4C/temp/",cur_rmv,"int_self.bed ~/Analyze4C/temp/",cur_rmv,"int.bed ",beta,sep=""),intern=TRUE)

				#record the f-measure and the iteration
				prf <- rep(0,3) #the precision, recall, and f-measure will eventually be recorded into here (respectively)
				counter2 <- 1 #counts the index of prf
				for(z in 1:(length(results))) #iterates over the lines that were recieved into results, from the function precisionRecallFmeasure.pl
				{
					spl_results <- strsplit(results[z],":") #splits the line where there is a ":", the number result should be what comes after that (possibly with a \n will be before the number)
					if((identical(spl_results[[1]],character(0)) == "FALSE")) #testing to see that the we don't get a 'character(0)' when doing strplit, this happens when the line is just \n
					{
						#putting each result number in prf
						prf[counter2] <- as.numeric(spl_results[[1]][2])
						counter2 <-  counter2 + 1
					}
				}
				prf_cyc <- rbind(prf_cyc,prf)
				
				#removing the contact and intersection files in case they are created multiple times
				system(paste("rm ~/Analyze4C/temp/",cur_rmv,"int_self.bed",sep=""))
				system(paste("rm ~/Analyze4C/temp/",cur_rmv,"int.bed",sep=""))
				system(paste("rm ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed",sep=""))
			}

			prec_med <- median(prf_cyc[,1])
			rec_med <- median(prf_cyc[,2])
			fm_med <- median(prf_cyc[,3])
			prf_med <- c(prec_med,rec_med,fm_med)
	
			fmeasure <- fm_med #recording the f-measure
			prf_all <- rbind(prf_all,prf_med) #adding the prf_med to prf_all
			
			#changing cur_rmv and getting the new value	
			counter1 <- counter1 + 1
			covPer_lim <- cur_rmv
			cur_rmv <- stp * counter1
			if(cur_rmv >= 100)
			{
				cat("\nthe percentage of coverage removal surpasses 100% and hasn't yet reached the f-measure limit of",fmeasure_lim,"\n\n")
				break
			}
		}
		
		#creating a data frame of the recall, precision, and f-measure data
		rownames(prf_all) <- seq(1,nrow(prf_all),1)
		covs <- c(0,seq(stp,covPer_lim,stp))
		covs <- paste(covs,"%",sep="")
		len <- length(covs)
		colnames(prf_all) <- c("precision","recall","f.measure")
		prf_prec <- data.frame(rep("precision",len),subset(prf_all,select="precision"),covs)
		names(prf_prec) <- c("type","result","coverage.removed")
		prf_rec <- data.frame(rep("recall",len),subset(prf_all,select="recall"),covs)
		names(prf_rec) <- c("type","result","coverage.removed")
		prf_fm <- data.frame(rep("f.measure",len),subset(prf_all,select="f.measure"),covs)
		names(prf_fm) <- c("type","result","coverage.removed")
		prf_final <- rbind(prf_prec,prf_rec,prf_fm)
		rownames(prf_final) <- seq(1,nrow(prf_final),1)	

		print(ggplot2::ggplot(data=prf_final,ggplot2::aes(x=coverage.removed,y=result,group=type)) + ggplot2::geom_line(ggplot2::aes(color=type)) + ggplot2::geom_point())				
	}
	else if(choice0 == 3) #limited by coverage
	{
		while(cov_lim >= cur_cov)
		{
			cat("\nerror: the coverage limit chosen is larger than the actual coverage.Please choose a coverage limit lower than the actual coverage.\n")
			cat("\nthe current actual coverage is ",cur_cov,"\n")
			cov_lim <- as.numeric(readline(prompt=cat("\nenter the lowest coverage you would like be reached (between 0 and 1):\n\n")))
			numberOfSteps <- as.integer(readline(prompt=cat("\nin how many steps would you like to get there?\n\n")))
		}
		cov_lim_stp <- (((cur_cov-cov_lim)/cur_cov)/numberOfSteps)*100 #calculating the percent of coverage removal for each step. the number here is a percentage (between 1 and 100)
		covs <- cur_cov #this will have the coverages after each coverage removal, and it starts with the current coverage (before removal)
		for(s in 1:numberOfSteps)
		{	
			cat("\n*****************************************************************\n\n")
			cur_rmv <- cov_lim_stp * s #every step we get cur_rmv which multiplies percentage of each step and the step number
			prf_cyc <- c() #gets the prf for each coverage removal percentage, for all its cycles
			for(p in 1:cyc)
			{
				cat("\ncycle ",p," for ",cur_rmv,"% coverage removal:\n\n",sep="")
				#coverage removal
				rmv_fraction <- cur_rmv/100 #changing the percentage into fraction for input to coverageChanger_chrChooser				
				out2 <- coverageChanger_chrChooser(rawData_rmv,rmv_fraction,cis_rmv,choice1,rem_chr)
				newData <- out2[[1]]
				
				#making file of removed coverage from raw data
				write.table(newData,paste("~/Analyze4C/temp/raw_",cur_rmv,"perRemoved.sgr",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
				
				#printing message according to chromosome and coverage amount removed
				if(choice1 == 1)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from the whole genome\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 2)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 3)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 4)
				{
					cat("\n\ncalculating p-scores for ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
				}
				
				#calculate p-scores for the rmv data
				if(choice2_rmv == 2)
				{
					#using the parameters above we calculate the p score with RE site window sizes only
					ps_temp <- pScore_nokb(newData,cis_rmv,vp.pos_rmv,erase_rmv,REs_rmv/2)
				}
				else
				{
					#using the parameters above we calculate the p score with RE site and bp window sizes
					ps_temp <- pScore(newData,cis_rmv,vp.pos_rmv,erase_rmv,wind_rmv,REs_rmv/2)
				}
				
				#write p score to file
				#write.table(ps_temp,paste("~/Analyze4C/temp/pscore_",cur_rmv,"perRemoved.sgr",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
				
				#ps_forF will get the p-score data that will get the contact and f-measure treatment
				ps_forF <- ps_temp #the default is all (choice4 == 1).
				if(choice4 == 2) #trans only
				{
					ps_forF <- ps_temp[ps_temp[,1]!=cis_rmv,]
				}
				else if(choice4 == 3) #cis only
				{
					ps_forF <- ps_temp[ps_temp[,1]==cis_rmv,]
				}
				else if(choice4 == 4) #specific chromosome
				{
					ps_forF <- ps_temp[ps_temp[,1]==f_chr,]
				}

				#printing message according to chromosome and coverage amount removed
				if(choice1 == 1) #whole genome
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from whole genome\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 2) #trans only
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from trans\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 3) #cis only
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from cis\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 4) #specific chromosome
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, from ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, from ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, from ",cur_rmv,"% removed coverage chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", from ",cur_rmv,"% removed coverage from chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
				}
							
				#creating the contact bands:
				if(choice2_rmv == 2) #if the data is rearranged
				{
					#checking if we are dealing with a file that has a translocation that isn't removed
					name_spl <- unlist(strsplit(file.name_rmv,"_")) #splits the raw data file name
					if(("rearranged" %in% name_spl) & !("removed" %in% name_spl)) #if the filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
					{
						conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData,flag=1,file.name=file.name_rmv) #the function 'contactBands_byChrom' with the input 1 will remove the added sections after applying the cutoff and then get contact bands
						conts_rmv <- conts_rmv[[1]]
					}
					else
					{
						conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
						conts_rmv <- conts_rmv[[1]]
					}
				}
				else if(choice3_rmv == 1) #if the cutoff option chosen is by a common cutoff
				{
					conts_rmv <- contactBands(ps_forF,RE_gap_rmv,bp_gap_rmv,cis_rmv,co_rmv,genome_sizes)
				}
				else if(choice3_rmv == 2) #if the cutoff option chosen is a percentage cutoff
				{
					conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData)
					conts_rmv <- conts_rmv[[1]]
				}	
				write.table(conts_rmv,paste("~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)	
				
				#getting the intersections of contacts of removed coverage with self
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed -b ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed -wo > ~/Analyze4C/temp/",cur_rmv,"int_self.bed",sep=""))
				
				#get the intersections with the static data
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_static.bed -b ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed -wo > ~/Analyze4C/temp/",cur_rmv,"int.bed",sep=""))
				
				#calculate f-measure
				results <- system(paste("perl ~/Analyze4C/proxy/precisionRecallFmeasure.pl ~/Analyze4C/temp/static_self.bed ~/Analyze4C/temp/",cur_rmv,"int_self.bed ~/Analyze4C/temp/",cur_rmv,"int.bed ",beta,sep=""),intern=TRUE)
				
				#record the f-measure and the iteration
				prf <- rep(0,3) #the precision, recall, and f-measure will eventually be recorded into here (respectively)
				counter2 <- 1 #counts the index of prf
				for(z in 1:(length(results))) #iterates over the lines that were recieved into results, from the function precisionRecallFmeasure.pl
				{
					spl_results <- strsplit(results[z],":") #splits the line where there is a ":", the number result should be what comes after that (possibly with a \n will be before the number)
					if((identical(spl_results[[1]],character(0)) == "FALSE")) #testing to see that the we don't get a 'character(0)' when doing strplit, this happens when the line is just \n
					{
						#putting each result number in prf
						prf[counter2] <- as.numeric(spl_results[[1]][2])
						counter2 <-  counter2 + 1
					}
				}
				prf_cyc <- rbind(prf_cyc,prf)
				
				#removing the contact and intersection files in case they are created multiple times
				system(paste("rm ~/Analyze4C/temp/",cur_rmv,"int_self.bed",sep=""))
				system(paste("rm ~/Analyze4C/temp/",cur_rmv,"int.bed",sep=""))
				system(paste("rm ~/Analyze4C/temp/conts_",cur_rmv,"perRemoved.bed",sep=""))
			}
			
			prec_med <- median(prf_cyc[,1])
			rec_med <- median(prf_cyc[,2])
			fm_med <- median(prf_cyc[,3])
			prf_med <- c(prec_med,rec_med,fm_med)

			prf_all <- rbind(prf_all,prf_med) #adding the prf_med to prf_all
			
			#getting the new coverage after removal and adding it to covs
			new_cov <- cur_cov*((100-cur_rmv)/100)
			covs <- c(covs,new_cov)
		}

		#creating a data frame of the recall, precision, and f-measure data
		rownames(prf_all) <- seq(1,nrow(prf_all),1)
		len <- length(covs)
		colnames(prf_all) <- c("precision","recall","f.measure")
		prf_prec <- data.frame(rep("precision",len),subset(prf_all,select="precision"),covs)
		names(prf_prec) <- c("type","result","coverage")
		prf_rec <- data.frame(rep("recall",len),subset(prf_all,select="recall"),covs)
		names(prf_rec) <- c("type","result","coverage")
		prf_fm <- data.frame(rep("f.measure",len),subset(prf_all,select="f.measure"),covs)
		names(prf_fm) <- c("type","result","coverage")
		prf_final <- rbind(prf_prec,prf_rec,prf_fm)
		rownames(prf_final) <- seq(1,nrow(prf_final),1)		
		
		#creating a plot of the recall, precision, and f-measure data. the x axis is backwards (the high values are on the left and lower on the right)
		print(ggplot2::ggplot(data=prf_final,ggplot2::aes(x=coverage,y=result,group=type)) + ggplot2::scale_x_reverse() + ggplot2::geom_line(ggplot2::aes(color=type)) + ggplot2::geom_point())		
	}
	else if(choice0 == 4) #limited by coverage removal percentage
	{
		for(cur_min_reads in seq(first_min,last_min_reads,by=min_read_stp))
		{	
			if(cur_min_reads > 1) #since we already got the basic intersections without removing coverage, it is the same like having the minimum of 1 reads. so there is no need doing it again
			{
				cat("\n*****************************************************************\n\n")
				#cat("\n",cur_min_reads,"maximum reads removal:\n\n",sep="")
				
				#reads removal

				out2 <- coverageChanger_byMinReads(rawData_rmv,cur_min_reads,cis_rmv,choice1,rem_chr)
				newData <- out2[[1]]
					
				#making file of removed coverage from raw data
				#write.table(newData,paste("~/Analyze4C/temp/raw_",cur_min_reads,"minReads.sgr",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
				
				#printing message according to chromosome and minimum reads
				if(choice1 == 1)
				{
					cat("\n\ncalculating p-scores for ",cur_min_reads," minimum reads on the whole genome\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 2)
				{
					cat("\n\ncalculating p-scores for ",cur_min_reads," minimum reads on trans\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 3)
				{
					cat("\n\ncalculating p-scores for ",cur_min_reads," minimum reads on cis\nplease wait...\n\n",sep="")
				}
				else if(choice1 == 4)
				{
					cat("\n\ncalculating p-scores for ",cur_min_reads," minimum reads on chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
				}
				
				#calculate p-scores for the rmv data
				if(choice2_rmv == 2)
				{
					#using the parameters above we calculate the p score with RE site window sizes only
					ps_temp <- pScore_nokb(newData,cis_rmv,vp.pos_rmv,erase_rmv,REs_rmv/2)
				}
				else
				{
					#using the parameters above we calculate the p score with RE site and bp window sizes
					ps_temp <- pScore(newData,cis_rmv,vp.pos_rmv,erase_rmv,wind_rmv,REs_rmv/2)
				}
				
				#write p score to file
				#write.table(ps_temp,paste("~/Analyze4C/temp/pscore_",cur_min_reads,"minReads.sgr",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
				
				#ps_forF will get the p-score data that will get the contact and f-measure treatment
				ps_forF <- ps_temp #the default is all (choice4 == 1).
				if(choice4 == 2) #trans only
				{
					ps_forF <- ps_temp[ps_temp[,1]!=cis_rmv,]
				}
				else if(choice4 == 3) #cis_rmv only
				{
					ps_forF <- ps_temp[ps_temp[,1]==cis_rmv,]
				}
				else if(choice4 == 4) #specific chromosome
				{
					ps_forF <- ps_temp[ps_temp[,1]==f_chr,]
				}

				#printing message according to chromosome and coverage amount removed
				if(choice1 == 1) #whole genome
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, with ",cur_min_reads," minimum reads on the whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, with ",cur_min_reads," minimum reads on the whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, with ",cur_min_reads," minimum reads on the whole genome\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", with ",cur_min_reads," minimum reads on the whole genome\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 2) #trans only
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, with ",cur_min_reads," minimum reads on trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, with ",cur_min_reads," minimum reads on trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, with ",cur_min_reads," minimum reads on trans\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", with ",cur_min_reads," minimum reads on trans\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 3) #cis only
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, with ",cur_min_reads," minimum reads on cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, with ",cur_min_reads," minimum reads on cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, with ",cur_min_reads," minimum reads on cis\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", with ",cur_min_reads," minimum reads on cis\nplease wait...\n\n",sep="")
					}
				}
				else if(choice1 == 4) #specific chromosome
				{
					if(choice4 == 1) #whole genome
					{
						cat("\n\ncreating contact bands for whole genome, with ",cur_min_reads," minimum reads on chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 2) #trans only
					{
						cat("\n\ncreating contact bands for trans, with ",cur_min_reads," minimum reads on chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 3) #cis only
					{
						cat("\n\ncreating contact bands for cis, with ",cur_min_reads," minimum reads on chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
					else if(choice4 == 4) #specific chromosome
					{
						cat("\n\ncreating contact bands for chromosome ",rem_chr,", with ",cur_min_reads," minimum reads on chromosome ",rem_chr,"\nplease wait...\n\n",sep="")
					}
				}
				
				#creating the contact bands:
				if(choice2_rmv == 2) #if the data is rearranged
				{
					#checking if we are dealing with a file that has a translocation that isn't removed
					name_spl <- unlist(strsplit(file.name_rmv,"_")) #splits the raw data file name
					if(("rearranged" %in% name_spl) & !("removed" %in% name_spl)) #if the filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
					{
						conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData,flag=1,file.name=file.name_rmv) #the function 'contactBands_byChrom' with the input 1 will remove the added sections after applying the cutoff and then get contact bands
						conts_rmv <- conts_rmv[[1]]
					}
					else
					{
						conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
						conts_rmv <- conts_rmv[[1]]
					}
				}
				else if(choice3_rmv == 1) #if the cutoff option chosen is by a common cutoff
				{
					conts_rmv <- contactBands(ps_forF,RE_gap_rmv,bp_gap_rmv,cis_rmv,co_rmv,genome_sizes)
				}
				else if(choice3_rmv == 2) #if the cutoff option chosen is a percentage cutoff
				{
					conts_rmv <- contactBands_byChrom(ps_forF,RE_gap_rmv,bp_gap_rmv,co_per_rmv,genome_sizes,rearranged_rawData)
					conts_rmv <- conts_rmv[[1]]
				}	
				write.table(conts_rmv,paste("~/Analyze4C/temp/conts_",cur_min_reads,"minReads.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)	
				
				#getting the intersections of contacts of removed coverage with self
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_",cur_min_reads,"minReads.bed -b ~/Analyze4C/temp/conts_",cur_min_reads,"minReads.bed -wo > ~/Analyze4C/temp/",cur_min_reads,"int_self.bed",sep=""))
				
				#get the intersections with the static data
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_static.bed -b ~/Analyze4C/temp/conts_",cur_min_reads,"minReads.bed -wo > ~/Analyze4C/temp/",cur_min_reads,"int.bed",sep=""))
				
				#calculate f-measure
				results <- system(paste("perl ~/Analyze4C/proxy/precisionRecallFmeasure.pl ~/Analyze4C/temp/static_self.bed ~/Analyze4C/temp/",cur_min_reads,"int_self.bed ~/Analyze4C/temp/",cur_min_reads,"int.bed ",beta,sep=""),intern=TRUE)
				
				#record the f-measure
				prf <- rep(0,3) #the precission, recall, and f-measure will eventually be recorded into here (respectively)
				counter2 <- 1 #counts the index of prf
				for(z in 1:(length(results))) #iterates over the lines that were recieved into results, from the function precisionRecallFmeasure.pl
				{
					spl_results <- strsplit(results[z],":") #splits the line where there is a ":", the number result should be what comes after that (possibly with a \n will be before the number)
					if((identical(spl_results[[1]],character(0)) == "FALSE")) #testing to see that the we don't get a 'character(0)' when doing strplit, this happens when the line is just \n
					{
						#putting each result number in prf
						prf[counter2] <- as.numeric(spl_results[[1]][2])
						counter2 <-  counter2 + 1
					}
				}
				
				#removing the contact and intersection files in case they are created multiple times
				system(paste("rm ~/Analyze4C/temp/",cur_min_reads,"int_self.bed",sep=""))
				system(paste("rm ~/Analyze4C/temp/",cur_min_reads,"int.bed",sep=""))
				system(paste("rm ~/Analyze4C/temp/conts_",cur_min_reads,"minReads.bed",sep=""))
				
				prf_all <- rbind(prf_all,prf) #adding the prf to prf_all
			}	
		}

		#creating a data frame of the recall, precision, and f-measure data
		rownames(prf_all) <- seq(1,nrow(prf_all),1)
		mins <- seq(first_min,last_min_reads,by=min_read_stp)
		len <- length(mins)
		colnames(prf_all) <- c("precision","recall","f.measure")
		prf_prec <- data.frame(rep("precision",len),subset(prf_all,select="precision"),mins)
		names(prf_prec) <- c("type","result","minimum.reads")
		prf_rec <- data.frame(rep("recall",len),subset(prf_all,select="recall"),mins)
		names(prf_rec) <- c("type","result","minimum.reads")
		prf_fm <- data.frame(rep("f.measure",len),subset(prf_all,select="f.measure"),mins)
		names(prf_fm) <- c("type","result","minimum.reads")
		prf_final <- rbind(prf_prec,prf_rec,prf_fm)
		rownames(prf_final) <- seq(1,nrow(prf_final),1)

		print(ggplot2::ggplot(data=prf_final,ggplot2::aes(x=minimum.reads,y=result,group=type)) + ggplot2::geom_line(ggplot2::aes(color=type)) + ggplot2::geom_point() + ggplot2::scale_x_continuous(breaks = seq(first_min,last_min_reads,by=min_read_stp)))		
	}

	#asking if to save the plot
	ans <- readline(prompt=cat("\n\nwould you like to save the plot?\ny/n\n\n"))
	if(ans == "y")
	{
		#getting the date and time in order to distinguish between file names of plots
		DandT1 <- toString(Sys.time())
		DandT2 <- gsub(" ","_",DandT1)
		DandT2 <- gsub(":","",DandT2)
		nm <- paste("coverage_removal_VS_Fmeasure_",DandT2,".png",sep="")
		wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
		ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
		ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)

		#getting parameters for coverageVSfmeasure_plots.txt
		#choice0 will give me lim_meth
		#covPer_lim, fmeasure_lim,cov_lim,last_min_reads will give the limit. it will go into limit
		#choice1 will give the removal location, into rem_loc
		#cov_lim_stp,stp, and min_read_stp will go into step
		#if pScore is calculated without bp, then wind will get NA
		#choice3 will say what cutoff method it is, will go into CO_method
		#co or co_per will go into CO
		#choice4 will say what chromosomes f-measure is calculated for, will go into Fmeasure_chromosomes
		#whatever doesn't have a value gets NA
		{
			#getting the type of coverage removal, the limit, and step
			if(choice0 == 1)
			{
				lim_meth <- "minimum coverage removal percentage"
				limit <- covPer_lim
				step <- stp
			}
			else if(choice0 == 2)
			{
				lim_meth <- "minimum f-measure"
				limit <- fmeasure_lim
				step <- stp
			}
			else if(choice0 == 3)
			{
				lim_meth <- "minimum coverage"
				limit <- cov_lim
				step <- cov_lim_stp
			}
			else if(choice0 == 4)
			{
				lim_meth <- "minimum number of reads per RE site"
				limit <- last_min_reads
				step <- min_read_stp
				cyc <- NA
			}
			
			#getting the chromosome the user would like coverage removal from
			if(choice1 == 1) 
			{
				rem_loc <- "whole genome"
			}
			else if(choice1 == 2)
			{
				rem_loc <- "trans"
			}
			else if(choice1 == 3)
			{
				rem_loc <- "cis"
			}	
			else if(choice1 == 4)
			{
				rem_loc <- toString(rem_chr)
			}	
			

			if(choice2_rmv == 2)#if pScore_nokb is activated for rmv data
			{
				wind_rmv <- NA
			}
			
			if(choice2_static == 2)#if pScore_nokb is activated for static data
			{
				wind_static <- NA
			}
			
			#cutoff method and cutoff for rmv and static
			if(choice3_rmv == 1)
			{
				CO_method_rmv <- "same CO for all"
				CO_rmv <- co_rmv
			}
			else if(choice3_rmv == 2)
			{
				CO_method_rmv <- "chromosome independent CO"
				CO_rmv <- co_per_rmv
			}
			
			#cutoff method and cutoff for static
			if(choice3_static == 1)
			{
				CO_method_static <- "same CO for all"
				CO_static <- co_static
			}
			else if(choice3_static == 2)
			{
				CO_method_static <- "chromosome independent CO"
				CO_static <- co_per_static
			}
			
			if(choice4 == 1)
			{
				Fmeasure_chromosomes <- "whole genome"
			}
			else if(choice4 == 2)
			{
				Fmeasure_chromosomes <- "trans"
			}	
			else if(choice4 == 3)
			{
				Fmeasure_chromosomes <- "cis"
			}
			else if(choice4 == 4)
			{
				Fmeasure_chromosomes <- toString(f_chr)
			}			
		}	
		
		#saving the parameters and details to coverageVSfmeasure_plots.txt
		coverageVSfmeasure_plots[nrow(coverageVSfmeasure_plots)+1,] <- c(sp1_rmv[[1]][1],REs_rmv,wind_rmv,RE_gap_rmv,bp_gap_rmv,CO_method_rmv,CO_rmv,sp1_static[[1]][1],REs_static,wind_static,RE_gap_static,bp_gap_static,CO_method_static,CO_static,lim_meth,limit,cyc,rem_loc,step,Fmeasure_chromosomes,beta,DandT1)
		#sorting the list of experiments by bait alphabetically (and sorting the row indices)
		coverageVSfmeasure_plots <- coverageVSfmeasure_plots[order(coverageVSfmeasure_plots$Coverage_Removal_Experiment),]
		rownames(coverageVSfmeasure_plots) <- seq(length=nrow(coverageVSfmeasure_plots))
		#adding the new data to the file (by erasing the old one and creating a new one)
		system("rm coverageVSfmeasure_plots.txt")
		write.table(coverageVSfmeasure_plots,"coverageVSfmeasure_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
	}			
}

#remove all the files from temp:
{
	system("rm ~/Analyze4C/temp/*")
}
}
