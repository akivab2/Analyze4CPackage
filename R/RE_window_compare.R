#' @export

#the point of this function is not to test how similar 2 files are compared to others, but rather compared to themselves using different window sizes
#the idea is to find the best window size (in regards to a specific cutoff and contact bands parameters) that will maximize the similarities between two files
#(i might need to think if i need to compare to others and how to do that, i will need to normalize the data. this might be needed to show that two of the same bait, different lanes, are more similar than same bait different cells)

#compares between 2 raw data files to find the most intersections by getting a range of window sizes and range of parameters for contact bands
#3 for loops -  first for range of windows for 1st raw data file, second for range of windows for 2nd raw data file, and 3rd for contact bands range
#1) get both raw data files
#2) get parameters of window sizes and step
#3) get parameters of contact bands and step
#4) create p-scores and contac bands (in temp)
#5) for all, trans, and each chromosome separately - create a file of contact bands and intersect
#6) create a table with all the data
#7) erase from temp
#8) plot data - x axis - 1st datas window, y axis - 2nd datas window, z axis - intersection number. have for each chromosome, all and trans

#to do:
#1) add a plane to the plots
#2) write an explanation on the top here
#3) change the size of the font of the y and x axises
#the f-measure might matter which file is first, so try both if needed (the first is the base, the second is compared to it)
RE_window_compare <- function(Experiments_4C,REwindowIntersections_plots,rearranged_rawData)
{
	#description of the function:
	cat("\nthis function compares 2 different files (they need to be with the same number of chromosomes and same cis chromosome)")
	cat("\nthe parameters being compared are precision, recall, f-measure, and the sum of the intersections of contact bands")
	cat("\nthe user will enter window sizes to be compared, the function will calculate p-scores per window size and create contact bands which will be compared by looking at the intersections between the datas")
	cat("\nall the data will be 3D plotted")
	cat("\nother paramters entered do not change throughout the function - meaning that cutoffs, f-measure beta, and contact bands parameters should be determined beforehand for the whole function")
	cat("\nthis function is advised to use when wanting to compare two of biological repeats (two lanes) to observe how similar they are\n\n")
	
	#checking that the 'temp' folder is empty, if not then empty it out
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

	cat("\nnote: the files you choose must have the same chromosomes\n\n")

	#choosing the 1st raw data file to work with:
	#getting the name of the raw data file the user wants to use
	repeat
	{
		choice2_1 <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the 1st raw data file from:\n\n1) original\n2) rearranged\n3) coverage removed\n\n")))
		if(choice2_1 == 2)
		{
			ls_files_1 <- system("ls ~/Analyze4C/rawData/rearranged",intern=TRUE)
		}
		else if(choice2_1 == 3)
		{
			ls_files_1 <- system("ls ~/Analyze4C/rawData/coverage_removed",intern=TRUE)
		}
		else
		{
			ls_files_1 <- system("ls ~/Analyze4C/rawData/original",intern=TRUE)
		}		
		file.name_1 <- readline(prompt=cat("These are the raw data files available\nwhich would you like to use?\n",ls_files_1,"\n",sep="\n"))
		ind_files_1 <- pmatch(file.name_1,ls_files_1,nomatch=0)
		if(ind_files_1 == 0)
		{
			cat("no such file exists.\nplease try again.\n\n")
		}
		else
		{
			spl1_1 <- unlist(strsplit(file.name_1,"_"))
			break
		}
	}
	#importing the data from the file
	if(choice2_1 == 2)
	{
		rawData_1 <- read.table(paste("~/Analyze4C/rawData/rearranged/",file.name_1,sep=""))
	}
	else if(choice2_1 == 3)
	{
		rawData_1 <- read.table(paste("~/Analyze4C/rawData/coverage_removed/",file.name_1,sep=""))
	}
	else
	{
		rawData_1 <- read.table(paste("~/Analyze4C/rawData/original/",file.name_1,sep=""))
	}
	
	#finding the specific raw data files information in Experiments_4C
	sp1_1 <- unlist(strsplit(file.name_1,"[.]"))
	sp2_1 <- strsplit(sp1_1[[1]][1],"_")
	out_1 <- findIn.Experiments_4C(sp2_1[[1]][1],sp2_1[[1]][2],sp2_1[[1]][3],Experiments_4C)
	
	#getting the cis_1 chromosome number
	cis_1 <- as.numeric(out_1[2])
	
	#getting the bait index
	vp.pos_1 <- as.numeric(out_1[3])
	
	#getting the amount of bp around the bait that will be erased
	erase_1 <- 1e6 #for now this will be the default. later i could have the user asked how many bp they want to be removed
	
	num_of_chroms <- length(unique(rawData_1[,1]))
	
	#choosing the 2nd raw data file to work with:
	#getting the name of the raw data file the user wants to use
	repeat
	{
		choice2_2 <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the 2nd raw data file from:\n\n1) original\n2) rearranged\n3) coverage removed\n\n")))
		if(choice2_2 == 2)
		{
			ls_files_2 <- system("ls ~/Analyze4C/rawData/rearranged",intern=TRUE)
		}
		else if(choice2_2 == 3)
		{
			ls_files_2 <- system("ls ~/Analyze4C/rawData/coverage_removed",intern=TRUE)
		}
		else
		{
			ls_files_2 <- system("ls ~/Analyze4C/rawData/original",intern=TRUE)
		}		
		file.name_2 <- readline(prompt=cat("These are the raw data files available\nwhich would you like to use?\n",ls_files_2,"\n",sep="\n"))
		ind_files_2 <- pmatch(file.name_2,ls_files_2,nomatch=0)
		if(ind_files_2 == 0)
		{
			cat("no such file exists.\nplease try again.\n\n")
		}
		else
		{
			spl1_2 <- unlist(strsplit(file.name_2,"_"))
			break
		}
	}
	#importing the data from the file
	if(choice2_2 == 2)
	{
		rawData_2 <- read.table(paste("~/Analyze4C/rawData/rearranged/",file.name_2,sep=""))
	}
	else if(choice2_2 == 3)
	{
		rawData_2 <- read.table(paste("~/Analyze4C/rawData/coverage_removed/",file.name_2,sep=""))
	}
	else
	{
		rawData_2 <- read.table(paste("~/Analyze4C/rawData/original/",file.name_2,sep=""))
	}
	
	#finding the specific raw data files information in Experiments_4C
	sp1_2 <- unlist(strsplit(file.name_2,"[.]"))
	sp2_2 <- strsplit(sp1_2[[1]][1],"_")
	out_2 <- findIn.Experiments_4C(sp2_2[[1]][1],sp2_2[[1]][2],sp2_2[[1]][3],Experiments_4C)
	
	#getting the cis_2 chromosome number
	cis_2 <- as.numeric(out_2[2])
	
	#getting the bait index
	vp.pos_2 <- as.numeric(out_2[3])
	
	#getting the amount of bp around the bait that will be erased
	erase_2 <- 1e6 #for now this will be the default. later i could have the user asked how many bp they want to be removed
		
	if(num_of_chroms != (length(unique(rawData_2[,1]))))
	{
		cat("\nerror: both raw data files don't have the same number of chromosomes!\n\n")
		return()
	}
	
	#getting all the parameters that will be used for all the contact band making:
	#p-score parameters:
	cat("\n*******************************************************************\n\n")
	cat("\n                       p-score parameters:\n\n")
	cat("\nfor file ",file.name_1,":\n",sep="")
	#get the window size in RE site
	min_RE1 <- as.integer(readline(prompt=cat("\nplease enter the size of the minimum window by RE sites (should be an even number):\n\n")))
	max_RE1 <- as.integer(readline(prompt=cat("\nplease enter the size of the maximum window by RE sites (should be an even number):\n\n")))
	step_RE1 <- as.integer(readline(prompt=cat("\nplease enter the size of the window step by RE sites (should be an even number):\n\n")))
	
	#getting the bp window size if not a rearranged data
	if(!("rearranged" %in% spl1_1) | ("removed" %in% spl1_1)) #if the raw data filename doesn't contain the word 'rearranged' or contains 'removed', meaning that it wasn't rearranged or that the rearranged sections that were rearranged sections were removed
	{	
		#get the proportion of RE per bp
		prp_1 <- as.integer(readline(prompt=cat("\nplease enter the proportion of REs per bp:\n\n")))
	}
	else
	{
		#if there are only windows of RE
		prp_1 <- "NA"
	}
		
	cat("\nfor file ",file.name_2,":\n",sep="")
	#get the window size in RE site
	min_RE2 <- as.integer(readline(prompt=cat("\nplease enter the size of the minimum window by RE sites (should be an even number):\n\n")))
	max_RE2 <- as.integer(readline(prompt=cat("\nplease enter the size of the maximum window by RE sites (should be an even number):\n\n")))
	step_RE2 <- as.integer(readline(prompt=cat("\nplease enter the size of the window step by RE sites (should be an even number):\n\n")))
	
	#getting the bp window size if not a rearranged data
	if(!("rearranged" %in% spl1_2) | ("removed" %in% spl1_2)) #if the raw data filename doesn't contain the word 'rearranged' or contains 'removed', meaning that it wasn't rearranged or that the rearranged sections that were rearranged sections were removed
	{	
		#get the proportion of RE per bp
		prp_2 <- as.integer(readline(prompt=cat("\nplease enter the proportion of REs per bp:\n\n")))
	}
	else
	{
		#if there are only windows of RE
		prp_2 <- "NA"	
	}

	#contact band parameters:
	cat("\n*******************************************************************\n\n")
	cat("\n                      conact bands parameters:\n\n")
	
	ls_genomes <- system("ls ~/Analyze4C/genomes/",intern=TRUE)
	genome_file.name <- readline(prompt=cat("\nchoose the file with genome sizes that you would like to use:\n",ls_genomes,"\n",sep="\n"))
	genome_sizes <- read.table(paste("~/Analyze4C/genomes/",genome_file.name,sep=""))
	
	#getting the different paramaters and inputs from user
	cat("\nfor file ",file.name_1,":\n",sep="")
	RE_gap_1 <- as.numeric(readline(prompt=cat("\nplease choose up to how many negative RE sites in between positive RE sites it is not considered a break in the contact band:\n\n")))
	bp_gap_1 <- as.numeric(readline(prompt=cat("\nplease choose the max number of bp that are allowed to be in between RE sites and not be considered a gap:\n\n")))

	#asking user if to create contact bands using a p-score cutoff for all chromosomes or a percent cutoff
	if(!("rearranged" %in% spl1_1) | ("removed" %in% spl1_1)) #if the raw data filename doesn't contain the word 'rearranged' or contains 'removed', meaning that it wasn't rearranged or that the rearranged sections that were rearranged sections were removed
	{	
		choice3_1 <- as.integer(readline(prompt=cat("\nenter the number of what type of contact band making you would like to do:\n\n1) all chromosomes with the same cutoff (in p-score)\n2) all chromosomes with the same cutoff percentage\n\n")))
		if(choice3_1 == 1)
		{
			co_1 <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))
			co_type1 <- "fixed"
			co1 <- co_1
		}
		else
		{
			co_per_1 <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff percent that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))
			co_type1 <- "top percent per chromosome"
			co1 <- co_per_1
		}
	}
	else #if the data chosen is rearranged, then we could only use a percentage cutoff (used later in contactBands_byChrom)
	{
		choice3_1 <- 2
		co_per_1 <- as.numeric(readline(prompt=cat("\nthe method of cutoff will be to apply a cutoff percentage\nplease enter a p-score cutoff percent that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))	
		co_type1 <- "top percent per chromosome"
		co1 <- co_per_1
	}
	
	#getting the different paramaters and inputs from user
	cat("\nfor file ",file.name_2,":\n",sep="")
	RE_gap_2 <- as.numeric(readline(prompt=cat("\nplease choose up to how many negative RE sites in between positive RE sites it is not considered a break in the contact band:\n\n")))
	bp_gap_2 <- as.numeric(readline(prompt=cat("\nplease choose the max number of bp that are allowed to be in between RE sites and not be considered a gap:\n\n")))
	
	#asking user if to create contact bands using a p-score cutoff for all chromosomes or a percent cutoff
	if(!("rearranged" %in% spl1_2) | ("removed" %in% spl1_2)) #if the raw data filename doesn't contain the word 'rearranged' or contains 'removed', meaning that it wasn't rearranged or that the rearranged sections that were rearranged sections were removed
	{	
		choice3_2 <- as.integer(readline(prompt=cat("\nenter the number of what type of contact band making you would like to do:\n\n1) all chromosomes with the same cutoff (in p-score)\n2) all chromosomes with the same cutoff percentage\n\n")))
		if(choice3_2 == 1)
		{
			co_2 <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))
			co_type2 <- "fixed"
			co2 <- co_2
		}
		else
		{
			co_per_2 <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff percent that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))
			co_type2 <- "top percent per chromosome"
			co2 <- co_per_2
		}
	}
	else #if the data chosen is rearranged, then we could only use a percentage cutoff (used later in contactBands_byChrom)
	{
		choice3_2 <- 2
		co_per_2 <- as.numeric(readline(prompt=cat("\nthe method of cutoff will be to apply a cutoff percentage\nplease enter a p-score cutoff percent that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))	
		co_type2 <- "top percent per chromosome"
		co2 <- co_per_2
	}

	#f-measure parameters:
	cat("\n*******************************************************************\n\n")
	cat("\n                       f-measure parameters:\n\n")
	beta <- -1
	while(beta <= 0)
	{
		beta <- as.numeric(readline(prompt=cat("\na beta between 0 and 1 will give more weight on precision. above 1 puts more weight on recall.\nplease enter a beta:\n\n")))
		if(beta <= 0)
		{
			cat("\nthe beta entered is incompatible. The beta must be over 0.\n\n")
		}
	}

	len <- (length(seq(min_RE1,max_RE1,by=step_RE1)))*(length(seq(min_RE2,max_RE2,by=step_RE2))) #getting the size of each table
	#calculating table sizes and creating empty tables for the precisions
	nms <- c("REs1","REs2","precision","chromosome")
	prec_all <- matrix(0,len,4)
	prec_all[,4] <- "all"
	colnames(prec_all) <- nms
	prec_all <- as.data.frame(prec_all,stringsAsFactors=FALSE)
	prf_all <- rep(0,3) #the precision, recall, and f-measure will eventually be recorded into here (respectively)
	prec_trans <- matrix(0,len,4)
	prec_trans[,4] <- "trans"
	colnames(prec_trans) <- nms
	prec_trans <- as.data.frame(prec_trans,stringsAsFactors=FALSE)
	prf_trans <- rep(0,3) #the precision, recall, and f-measure will eventually be recorded into here (respectively)	
	prec_chroms <- list()
	prf_chroms <- list()
	for(l in 1:num_of_chroms)
	{
		prec_chroms[[l]] <- matrix(0,len,4)
		prec_chroms[[l]][,4] <- l
		colnames(prec_chroms[[l]]) <- nms
		prec_chroms[[l]] <- as.data.frame(prec_chroms[[l]],stringsAsFactors=FALSE)
		prf_chroms[[l]] <- rep(0,3) #the precision, recall, and f-measure will eventually be recorded into here (respectively)			
	}

	#calculating table sizes and creating empty tables for the recalls
	nms <- c("REs1","REs2","recall","chromosome")
	rcll_all <- matrix(0,len,4)
	rcll_all[,4] <- "all"
	colnames(rcll_all) <- nms
	rcll_all <- as.data.frame(rcll_all,stringsAsFactors=FALSE)
	rcll_trans <- matrix(0,len,4)
	rcll_trans[,4] <- "trans"
	colnames(rcll_trans) <- nms
	rcll_trans <- as.data.frame(rcll_trans,stringsAsFactors=FALSE)
	rcll_chroms <- list()
	for(l in 1:num_of_chroms)
	{
		rcll_chroms[[l]] <- matrix(0,len,4)
		rcll_chroms[[l]][,4] <- l
		colnames(rcll_chroms[[l]]) <- nms
		rcll_chroms[[l]] <- as.data.frame(rcll_chroms[[l]],stringsAsFactors=FALSE)
	}

	#calculating table sizes and creating empty tables for the f-measures
	nms <- c("REs1","REs2","Fmeasure","chromosome")
	fms_all <- matrix(0,len,4)
	fms_all[,4] <- "all"
	colnames(fms_all) <- nms
	fms_all <- as.data.frame(fms_all,stringsAsFactors=FALSE)
	fms_trans <- matrix(0,len,4)
	fms_trans[,4] <- "trans"
	colnames(fms_trans) <- nms
	fms_trans <- as.data.frame(fms_trans,stringsAsFactors=FALSE)
	fms_chroms <- list()
	for(l in 1:num_of_chroms)
	{
		fms_chroms[[l]] <- matrix(0,len,4)
		fms_chroms[[l]][,4] <- l
		colnames(fms_chroms[[l]]) <- nms
		fms_chroms[[l]] <- as.data.frame(fms_chroms[[l]],stringsAsFactors=FALSE)
	}
	
	#calculating table sizes and creating empty tables for the intersections
	nms <- c("REs1","REs2","intersections","chromosome")
	sum_ints_all <- matrix(0,len,4)
	sum_ints_all[,4] <- "all"
	colnames(sum_ints_all) <- nms
	sum_ints_all <- as.data.frame(sum_ints_all,stringsAsFactors=FALSE)
	sum_ints_trans <- matrix(0,len,4)
	sum_ints_trans[,4] <- "trans"
	colnames(sum_ints_trans) <- nms
	sum_ints_trans <- as.data.frame(sum_ints_trans,stringsAsFactors=FALSE)
	sum_ints_chroms <- list()
	for(l in 1:num_of_chroms)
	{
		sum_ints_chroms[[l]] <- matrix(0,len,4)
		sum_ints_chroms[[l]][,4] <- l
		colnames(sum_ints_chroms[[l]]) <- nms
		sum_ints_chroms[[l]] <- as.data.frame(sum_ints_chroms[[l]],stringsAsFactors=FALSE)
	}

	flag <- 0 #the first time we do a loop of all windows for the 2nd raw data flag is 0, the second time we don't need to calculate them since we already did, so flag is 1
	flag_rmv1 <- 0
	flag_rmv2 <- 0
	ind <- 1 #the index for the sum_ints tables
	in_lns1 <- 0
	in_lns2 <- 0
	for(i in seq(min_RE1,max_RE1,by=step_RE1))
	{
		if(flag_rmv1 == 1)
		{
			cat("\nthese are the lines that were removed previously from the intersection output (if no line numbers appear under 'lines:' then there aren't any lines that need to be removed in that chromosome):\n")
			for(yy in 1:length(conts_out1[[2]]))
			{
				cat("\nchromosome ",conts_out1[[2]][[yy]][[1]],":\n",sep="")
				cat("\nlines:\n")
				cat(conts_out1[[2]][[yy]][[2]])
				cat("\n")
			}
			
			cat("\nwould you like to use these data lines for removing repeating intersections from the out.bed file, for any intersections done with the 1st file?\n")
			ans5 <- readline(prompt=cat("\ny/n\n\n"))
			if(ans5 == "y")
			{
				flag_rmv1 <- 2
			}
		}
		
		#calculating p-scores for the 1st raw data (if a translocated data is chosen this must be with no kb and before translocation areas are removed):
		{	
			cat("\n\ncalculating p-scores for",file.name_1,"with RE windows of",i,"\nplease wait...\n\n")
			if(("rearranged" %in% spl1_1) & !("removed" %in% spl1_1)) #if the raw data filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
			{	
				#using the parameters above we calculate the p score with RE site window sizes only
				ps_1 <- pScore_nokb(rawData_1,cis_1,vp.pos_1,erase_1,i/2)
			}
			else
			{
				#using the parameters above we calculate the p score with RE site and bp window sizes
				ps_1 <- pScore(rawData_1,cis_1,vp.pos_1,erase_1,(prp_1*i),(i/2))
			}
		}
		
		#creating contact bands files:
		#applying cutoff and getting contact bands for the static data(use contactBands_byChrom_translocatedRemoval if translocation, this will first apply cutoff then remove translocations and then get contact bands):
		{
			cat("\n\ncreating contact bands for",file.name_1,"with RE windows of",i,"\nplease wait...\n\n")
			#creating the contact bands by asking for cutoffs first
			if(("rearranged" %in% spl1_1)) #if the raw data filename contains the word 'rearranged'
			{	
				#checking if we are dealing with a file that has a translocation that isn't removed
				if(!("removed" %in% spl1_1)) #if the filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
				{
					conts_out1 <- contactBands_byChrom(ps_1,RE_gap_1,bp_gap_1,co_per_1,genome_sizes,rearranged_rawData,flag_rmv1,in_lns1,1,file.name_1) #the function 'contactBands_byChrom' with the input 1 will remove the added sections after applying the cutoff and then get contact bands
					conts_1 <- conts_out1[[1]]
					in_lns1 <- conts_out1[[2]]
					if(flag_rmv1 == 0)
					{
						flag_rmv1 <- 1
					}	
				}
				else
				{
					#conts_out1 <- contactBands_byChrom(ps_1,RE_gap_1,bp_gap_1,co_per_1,genome_sizes,rearranged_rawData,flag_rmv1,in_lns1) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
					conts_out1 <- contactBands_byChrom(ps_1,RE_gap_1,bp_gap_1,co_per_1,genome_sizes,rearranged_rawData) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
					conts_1 <- conts_out1[[1]]
					in_lns1 <- conts_out1[[2]]
					#if(flag_rmv1 == 0)
					#{
					#	flag_rmv1 <- 1
					#}	
				}
			}
			else if(choice3_1 == 1) #if the cutoff option chosen is by a common cutoff
			{
				conts_1 <- contactBands(ps_1,RE_gap_1,bp_gap_1,cis_1,co_1,genome_sizes)
			}
			else if(choice3_1 == 2) #if the cutoff option chosen is a percentage cutoff
			{
				#conts_out1 <- contactBands_byChrom(ps_1,RE_gap_1,bp_gap_1,co_per_1,genome_sizes,rearranged_rawData,flag_rmv1,in_lns1)
				conts_out1 <- contactBands_byChrom(ps_1,RE_gap_1,bp_gap_1,co_per_1,genome_sizes,rearranged_rawData)				
				conts_1 <- conts_out1[[1]]
				in_lns1 <- conts_out1[[2]]
				#if(flag_rmv1 == 0)
				#{
				#	flag_rmv1 <- 1
				#}	
			}
			write.table(conts_1,paste("~/Analyze4C/temp/conts_1_",i,"_all.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
			write.table(conts_1[conts_1[,1]!=cis_1,],paste("~/Analyze4C/temp/conts_1_",i,"_trans.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
			
			for(n in 1:num_of_chroms)
			{
				write.table(conts_1[conts_1[,1]==n,],paste("~/Analyze4C/temp/conts_1_",i,"_",n,".bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)	
			}	
		}
		
		for(j in seq(min_RE2,max_RE2,by=step_RE2))
		{
					
			if(flag_rmv2 == 1)
			{
				cat("\nthese are the lines that were removed previously from the intersection output (if no line numbers appear under 'lines:' then there aren't any lines that need to be removed in that chromosome):\n")
				for(xx in 1:length(conts_out2[[2]]))
				{
					cat("\nchromosome ",conts_out2[[2]][[xx]][[1]],":\n",sep="")
					cat("\nlines:\n")
					cat(conts_out2[[2]][[xx]][[2]])
					cat("\n")
				}
				cat("\nwould you like to use these data lines for removing repeating intersections from the out.bed file, for any intersections done with the 2nd file?\n")
				ans6 <- readline(prompt=cat("\ny/n\n\n"))
				if(ans6 == "y")
				{
					flag_rmv2 <- 2
				}
			}

			#create the p-score and contact bands files only on the first round of for loop, erase them only at the end of everything:
			if(flag == 0)
			{
				#calculating p-scores for the 2nd raw data (if a translocated data is chosen this must be with no kb and before translocation areas are removed):
				{	
					cat("\n\ncalculating p-scores for",file.name_2,"with RE windows of",j,"\nplease wait...\n\n")
					if(("rearranged" %in% spl1_2) & !("removed" %in% spl1_2)) #if the raw data filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
					{	
						#using the parameters above we calculate the p score with RE site window sizes only
						ps_2 <- pScore_nokb(rawData_2,cis_2,vp.pos_2,erase_2,j/2)
					}
					else
					{
						#using the parameters above we calculate the p score with RE site and bp window sizes
						ps_2 <- pScore(rawData_2,cis_2,vp.pos_2,erase_2,(prp_2*j),(j/2))
					}
				}
			
				#creating contact bands files:
				#applying cutoff and getting contact bands for the static data(use contactBands_byChrom_translocatedRemoval if translocation, this will first apply cutoff then remove translocations and then get contact bands):
				{
					cat("\n\ncreating contact bands for",file.name_2,"with RE windows of",j,"\nplease wait...\n\n")
					#creating the contact bands by asking for cutoffs first
					if(("rearranged" %in% spl1_2)) #if the raw data filename contains the word 'rearranged'
					{	
						#checking if we are dealing with a file that has a translocation that isn't removed
						if(!("removed" %in% spl1_2)) #if the filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
						{
							conts_out2 <- contactBands_byChrom(ps_2,RE_gap_2,bp_gap_2,co_per_2,genome_sizes,rearranged_rawData,flag_rmv2,in_lns2,1,file.name_2) #the function 'contactBands_byChrom' with the input 1 will remove the added sections after applying the cutoff and then get contact bands
							conts_2 <- conts_out2[[1]]
							in_lns2 <- conts_out2[[2]]
							if(flag_rmv2 == 0)
							{
								flag_rmv2 <- 1
							}	
						}
						else
						{
							#conts_out2 <- contactBands_byChrom(ps_2,RE_gap_2,bp_gap_2,co_per_2,genome_sizes,rearranged_rawData,flag_rmv2,in_lns2) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
							conts_out2 <- contactBands_byChrom(ps_2,RE_gap_2,bp_gap_2,co_per_2,genome_sizes,rearranged_rawData) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
							conts_2 <- conts_out2[[1]]
							in_lns2 <- conts_out2[[2]]
							#if(flag_rmv2 == 0)
							#{	
							#	flag_rmv2 <- 1
							#}	
						}
					}
					else if(choice3_2 == 1) #if the cutoff option chosen is by a common cutoff
					{
						conts_2 <- contactBands(ps_2,RE_gap_2,bp_gap_2,cis_2,co_2,genome_sizes)
					}
					else if(choice3_2 == 2) #if the cutoff option chosen is a percentage cutoff
					{
						#conts_out2 <- contactBands_byChrom(ps_2,RE_gap_2,bp_gap_2,co_per_2,genome_sizes,rearranged_rawData,flag_rmv2,in_lns2)
						conts_out2 <- contactBands_byChrom(ps_2,RE_gap_2,bp_gap_2,co_per_2,genome_sizes,rearranged_rawData)						
						conts_2 <- conts_out2[[1]]
						in_lns2 <- conts_out2[[2]]
						#if(flag_rmv2 == 0)
						#{
						#	flag_rmv2 <- 1
						#}	
					}	
					write.table(conts_2,paste("~/Analyze4C/temp/conts_2_",j,"_all.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
					write.table(conts_2[conts_2[,1]!=cis_2,],paste("~/Analyze4C/temp/conts_2_",j,"_trans.bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
					
					for(n in 1:num_of_chroms)
					{
						write.table(conts_2[conts_2[,1]==n,],paste("~/Analyze4C/temp/conts_2_",j,"_",n,".bed",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)	
					}	
				}
			}
			
			#getting precision, recall, f-measure, and sum of intersections for all:
			#intersecting files and recording in table
			system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_1_",i,"_all.bed -b ~/Analyze4C/temp/conts_1_",i,"_all.bed -wo > ~/Analyze4C/temp/conts_1_",i,"_all_self.bed",sep="")) #self for f-measure
			system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_2_",j,"_all.bed -b ~/Analyze4C/temp/conts_2_",j,"_all.bed -wo > ~/Analyze4C/temp/conts_2_",j,"_all_self.bed",sep="")) #self for f-measure
			system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_1_",i,"_all.bed -b ~/Analyze4C/temp/conts_2_",j,"_all.bed -wo > ~/Analyze4C/temp/ints_",i,"_",j,"_all.bed",sep=""))
			#in order to perform f-measure, the cis for both files must be the same
			if(cis_1 == cis_2)
			{
				fmResults_all <- system(paste("perl ~/Analyze4C/proxy/precisionRecallFmeasure.pl ~/Analyze4C/temp/conts_1_",i,"_all_self.bed ~/Analyze4C/temp/conts_2_",j,"_all_self.bed ~/Analyze4C/temp/ints_",i,"_",j,"_all.bed ",beta,sep=""),intern=TRUE)

				#record the f-measure and the iteration
				counter2 <- 1 #counts the index of prf_all
				for(z in 1:(length(fmResults_all))) #iterates over the lines that were recieved into fmResults_all, from the function precisionRecallFmeasure.pl
				{
					spl_results <- strsplit(fmResults_all[z],":") #splits the line where there is a ":", the number result should be what comes after that (possibly with a \n will be before the number)
					if((identical(spl_results[[1]],character(0)) == "FALSE")) #testing to see that the we don't get a 'character(0)' when doing strplit, this happens when the line is just \n
					{
						#putting each result number in prf_all
						prf_all[counter2] <- as.numeric(spl_results[[1]][2])
						counter2 <-  counter2 + 1
					}
				}
			}
			#putting the recall, precision, f-measure, and sum of intersections data in the tables for all
			prec_all$REs1[ind] <- i
			prec_all$REs2[ind] <- j
			prec_all$precision[ind] <- prf_all[1]
			rcll_all$REs1[ind] <- i
			rcll_all$REs2[ind] <- j
			rcll_all$recall[ind] <- prf_all[2]
			fms_all$REs1[ind] <- i
			fms_all$REs2[ind] <- j
			fms_all$Fmeasure[ind] <- prf_all[3]	

			if(prf_all[3] == 0)
			{
				sum_all <- 0
			}
			else
			{
				all <- read.delim(paste("~/Analyze4C/temp/ints_",i,"_",j,"_all.bed",sep=""),header=FALSE,quote="",stringsAsFactors=FALSE)
				sum_all <- sum(all[,7])
			}			
			sum_ints_all$REs1[ind] <- i
			sum_ints_all$REs2[ind] <- j
			sum_ints_all$intersections[ind] <- sum_all
			
			#getting precision, recall, f-measure, and sum of intersections for trans:
			#intersecting files and recording in table			
			system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_1_",i,"_trans.bed -b ~/Analyze4C/temp/conts_1_",i,"_trans.bed -wo > ~/Analyze4C/temp/conts_1_",i,"_trans_self.bed",sep="")) #self for f-measure
			system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_2_",j,"_trans.bed -b ~/Analyze4C/temp/conts_2_",j,"_trans.bed -wo > ~/Analyze4C/temp/conts_2_",j,"_trans_self.bed",sep="")) #self for f-measure				
			system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_1_",i,"_trans.bed -b ~/Analyze4C/temp/conts_2_",j,"_trans.bed -wo > ~/Analyze4C/temp/ints_",i,"_",j,"_trans.bed",sep=""))
			#in order to perform f-measure, the cis for both files must be the same
			if(cis_1 == cis_2)
			{
				fmResults_trans <- system(paste("perl ~/Analyze4C/proxy/precisionRecallFmeasure.pl ~/Analyze4C/temp/conts_1_",i,"_trans_self.bed ~/Analyze4C/temp/conts_2_",j,"_trans_self.bed ~/Analyze4C/temp/ints_",i,"_",j,"_trans.bed ",beta,sep=""),intern=TRUE)

				#record the f-measure and the iteration
				counter2 <- 1 #counts the index of prf_all
				for(z in 1:(length(fmResults_trans))) #iterates over the lines that were recieved into fmResults_trans, from the function precisionRecallFmeasure.pl
				{
					spl_results <- strsplit(fmResults_trans[z],":") #splits the line where there is a ":", the number result should be what comes after that (possibly with a \n will be before the number)
					if((identical(spl_results[[1]],character(0)) == "FALSE")) #testing to see that the we don't get a 'character(0)' when doing strplit, this happens when the line is just \n
					{
						#putting each result number in prf_all
						prf_trans[counter2] <- as.numeric(spl_results[[1]][2])
						counter2 <-  counter2 + 1
					}
				}				
			}
			
			#putting the recall, precision, f-measure, and sum of intersections data in the tables for trans
			prec_trans$REs1[ind] <- i
			prec_trans$REs2[ind] <- j
			prec_trans$precision[ind] <- prf_trans[1]
			rcll_trans$REs1[ind] <- i
			rcll_trans$REs2[ind] <- j
			rcll_trans$recall[ind] <- prf_trans[2]
			fms_trans$REs1[ind] <- i
			fms_trans$REs2[ind] <- j
			fms_trans$Fmeasure[ind] <- prf_trans[3]	

			if(prf_trans[3] == 0)
			{
				sum_trans <- 0
			}
			else
			{
				trans <- read.delim(paste("~/Analyze4C/temp/ints_",i,"_",j,"_trans.bed",sep=""),header=FALSE,quote="",stringsAsFactors=FALSE)
				sum_trans <- sum(trans[,7])
			}			
			sum_ints_trans$REs1[ind] <- i
			sum_ints_trans$REs2[ind] <- j
			sum_ints_trans$intersections[ind] <- sum_trans

			#getting precision, recall, f-measure, and sum of intersections for chroms:	
			chroms <- list()
			sum_chroms <- rep(0,num_of_chroms)
			for(m in 1:num_of_chroms)
			{
				#intersecting files and recording in table
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_1_",i,"_",m,".bed -b ~/Analyze4C/temp/conts_1_",i,"_",m,".bed -wo > ~/Analyze4C/temp/conts_1_",i,"_",m,"_self.bed",sep="")) #self for f-measure
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_2_",j,"_",m,".bed -b ~/Analyze4C/temp/conts_2_",j,"_",m,".bed -wo > ~/Analyze4C/temp/conts_2_",j,"_",m,"_self.bed",sep="")) #self for f-measure
				system(paste("bedtools intersect -a ~/Analyze4C/temp/conts_1_",i,"_",m,".bed -b ~/Analyze4C/temp/conts_2_",j,"_",m,".bed -wo > ~/Analyze4C/temp/ints_",i,"_",j,"_",m,".bed",sep=""))
				#in order to perform f-measure, the cis for both files must be the same
				if(cis_1 == cis_2)
				{
					fmResults_chroms <- system(paste("perl ~/Analyze4C/proxy/precisionRecallFmeasure.pl ~/Analyze4C/temp/conts_1_",i,"_",m,"_self.bed ~/Analyze4C/temp/conts_2_",j,"_",m,"_self.bed ~/Analyze4C/temp/ints_",i,"_",j,"_",m,".bed ",cis_1," ",beta,sep=""),intern=TRUE)

					#record the f-measure and the iteration
					counter2 <- 1 #counts the index of prf_all
					for(z in 1:(length(fmResults_chroms))) #iterates over the lines that were recieved into fmResults_chroms, from the function precisionRecallFmeasure.pl
					{
						spl_results <- strsplit(fmResults_chroms[z],":") #splits the line where there is a ":", the number result should be what comes after that (possibly with a \n will be before the number)
						if((identical(spl_results[[1]],character(0)) == "FALSE")) #testing to see that the we don't get a 'character(0)' when doing strplit, this happens when the line is just \n
						{
							#putting each result number in prf_all
							prf_chroms[[m]][counter2] <- as.numeric(spl_results[[1]][2])
							counter2 <-  counter2 + 1
						}
					}
				}				
				#putting the recall, precision, f-measure, and sum of intersections data in the tables for all
				prec_chroms[[m]]$REs1[ind] <- i
				prec_chroms[[m]]$REs2[ind] <- j
				prec_chroms[[m]]$precision[ind] <- prf_chroms[[m]][1]
				rcll_chroms[[m]]$REs1[ind] <- i
				rcll_chroms[[m]]$REs2[ind] <- j
				rcll_chroms[[m]]$recall[ind] <- prf_chroms[[m]][2]
				fms_chroms[[m]]$REs1[ind] <- i
				fms_chroms[[m]]$REs2[ind] <- j
				fms_chroms[[m]]$Fmeasure[ind] <- prf_chroms[[m]][3]	
				
				if(prf_chroms[[m]][3] == 0)
				{
					sum_chroms[m] <- 0
				}
				else
				{
					chroms[[m]] <- read.delim(paste("~/Analyze4C/temp/ints_",i,"_",j,"_",m,".bed",sep=""),header=FALSE,quote="",stringsAsFactors=FALSE)
					sum_chroms[m] <- sum(chroms[[m]][,7])
				}	
				sum_ints_chroms[[m]]$REs1[ind] <- i
				sum_ints_chroms[[m]]$REs2[ind] <- j
				sum_ints_chroms[[m]]$intersections[ind] <- sum_chroms[m]			
			}
						
			ind <- ind + 1
		}
		flag <- 1
	}
	
	#getting the date and time in order to distinguish between file names of plots
	DandT1 <- toString(Sys.time())
	DandT2 <- gsub(" ","_",DandT1)
	DandT2 <- gsub(":","",DandT2)
	
	#3d scatterplots:
	
	#creating the titles of the axis
	if(choice3_1 == 1)
	{
		xname <- paste("window sizes - ",co1," p-score cutoff\n",file.name_1,sep="")
	}
	else
	{
		xname <- paste("window sizes - ",co1*100,"% cutoff\n",file.name_1,sep="")
	}

	if(choice3_2 == 1)
	{
		yname <- paste("window sizes - ",co2," p-score cutoff\n",file.name_2,sep="")
	}
	else
	{
		yname <- paste("window sizes - ",co2*100,"% cutoff\n",file.name_2,sep="")
	}
	
	while(dev.cur() != 1)
	{
		dev.off()
	}
	save_flag <- 0 #if at least one of the plots is saved, then we will save the details of this run in the file REwindowIntersections_plots.txt, save_flag will be 1 if so

	ans1 <- readline(prompt=cat("\nwould you like to plot all?\ny/n\n\n"))
	if(ans1 == "y")
	{
		#plotting and saving all

		cat("\ncreating a plot for the precision data for 'all'\nplease wait...\n\n")
		scatterplot3d::scatterplot3d(prec_all[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. precision - all",xlab=xname,ylab=yname,zlab="precision")
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_precision_all_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_precision_all_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 25)
			scatterplot3d::scatterplot3d(prec_all[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. precision - all",xlab=xname,ylab=yname,zlab="precision")
			dev.off()
			save_flag <- 1
		}
		cat("\ncreating a plot for the recall data for 'all'\nplease wait...\n\n")
		while(dev.cur() != 1)
		{
			dev.off()
		}	
		scatterplot3d::scatterplot3d(rcll_all[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. recall - all",xlab=xname,ylab=yname,zlab="recall")
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_recall_all_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_recall_all_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 25)		
			scatterplot3d::scatterplot3d(rcll_all[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. recall - all",xlab=xname,ylab=yname,zlab="recall")
			dev.off()
			save_flag <- 1
		}
		cat("\ncreating a plot for the f-measure data for 'all'\nplease wait...\n\n")
		while(dev.cur() != 1)
		{
			dev.off()
		}	
		scatterplot3d::scatterplot3d(fms_all[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. f-measure - all",xlab=xname,ylab=yname,zlab="f-measure")
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_Fmeasure_all_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_Fmeasure_all_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 20)		
			scatterplot3d::scatterplot3d(fms_all[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. f-measure - all",xlab=xname,ylab=yname,zlab="f-measure")	
			dev.off()
			save_flag <- 1
		}
		cat("\ncreating a plot for the sum of intersections data for 'all'\nplease wait...\n\n")
		while(dev.cur() != 1)
		{
			dev.off()
		}	
		scatterplot3d::scatterplot3d(sum_ints_all[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. sumOfIntersections - all",xlab=xname,ylab=yname,zlab="sum of intersections (bp)")
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_sumOfIntersections_all_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_sumOfIntersections_all_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 20)		
			scatterplot3d::scatterplot3d(sum_ints_all[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. sumOfIntersections - all",xlab=xname,ylab=yname,zlab="sum of intersections (bp)")
			dev.off()
			save_flag <- 1
		}		
	}	
	#plotting and saving trans:
	
	ans2 <- readline(prompt=cat("\nwould you like to plot trans?\ny/n\n\n"))
	if(ans2 == "y")
	{
		cat("\ncreating a plot for the precision data for 'trans'\nplease wait...\n\n")
		while(dev.cur() != 1)
		{
			dev.off()
		}	
		scatterplot3d::scatterplot3d(prec_trans[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. precision - trans",xlab=xname,ylab=yname,zlab="precision")
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_precision_trans_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_precision_trans_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 25)		
			scatterplot3d::scatterplot3d(prec_trans[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. precision - trans",xlab=xname,ylab=yname,zlab="precision")
			dev.off()
			save_flag <- 1
		}
		cat("\ncreating a plot for the recall data for 'trans'\nplease wait...\n\n")
		while(dev.cur() != 1)
		{
			dev.off()
		}	
		scatterplot3d::scatterplot3d(rcll_trans[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. recall - trans",xlab=xname,ylab=yname,zlab="recall")
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_recall_trans_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_recall_trans_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 25)		
			scatterplot3d::scatterplot3d(rcll_trans[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. recall - trans",xlab=xname,ylab=yname,zlab="recall")
			dev.off()
			save_flag <- 1
		}
		cat("\ncreating a plot for the f-measure data for 'trans'\nplease wait...\n\n")
		while(dev.cur() != 1)
		{
			dev.off()
		}	
		scatterplot3d::scatterplot3d(fms_trans[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. f-measure - trans",xlab=xname,ylab=yname,zlab="f-measure")
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")	
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_Fmeasure_trans_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_Fmeasure_trans_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 20)	
			scatterplot3d::scatterplot3d(fms_trans[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. f-measure - trans",xlab=xname,ylab=yname,zlab="f-measure")	
			dev.off()
			save_flag <- 1
		}
		cat("\ncreating a plot for the sum of intersections data for 'trans'\nplease wait...\n\n")
		while(dev.cur() != 1)
		{
			dev.off()
		}
		scatterplot3d::scatterplot3d(sum_ints_trans[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. sumOfIntersections - trans",xlab=xname,ylab=yname,zlab="sum of intersections (bp)")
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_sumOfIntersections_trans_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_sumOfIntersections_trans_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 20)	
			scatterplot3d::scatterplot3d(sum_ints_trans[,1:3],pch=16,highlight.3d=TRUE,type="h",main="RE window vs. sumOfIntersections - trans",xlab=xname,ylab=yname,zlab="sum of intersections (bp)")
			dev.off()
			save_flag <- 1
		}
	}
	
	ans3 <- readline(prompt=cat("\nwould you like to plot each chromosome separately?\ny/n\n\n"))

	#chroms:
	prec_chroms_all <- c()
	rcll_chroms_all <- c()
	fms_chroms_all <- c()
	sum_ints_chroms_all <- c()
	for(q in 1:num_of_chroms)
	{
		if(ans3 == "y")
		{
			#plotting and saving each chromosome seperately:
			
			cat("\ncreating a plot for the precision data for chromosome",q,"\nplease wait...\n\n")
			while(dev.cur() != 1)
			{
				dev.off()
			}		
			scatterplot3d::scatterplot3d(prec_chroms[[q]][,1:3],pch=16,highlight.3d=TRUE,type="h",main=paste("RE window vs. precision - chroms",q),xlab=xname,ylab=yname,zlab="precision")
			ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
			if(ans_save == "y")
			{
				#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
				png_name <- paste("REwindowCompare_precision_chrom",q,"_",DandT2,".png",sep="")
				#png_name <- paste("REwindowCompare_precision_chrom",q,"_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
				#if(nchar(png_name) > 255)
				#{
				#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
				#	png_name <- substring(png_name,nm_length_rmv)
				#}
				png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 25)	
				scatterplot3d::scatterplot3d(prec_chroms[[q]][,1:3],pch=16,highlight.3d=TRUE,type="h",main=paste("RE window vs. precision - chroms",q),xlab=xname,ylab=yname,zlab="precision")
				dev.off()
				save_flag <- 1
			}
			cat("\ncreating a plot for the recall data for chromosome",q,"\nplease wait...\n\n")
			while(dev.cur() != 1)
			{
				dev.off()
			}		
			scatterplot3d::scatterplot3d(rcll_chroms[[q]][,1:3],pch=16,highlight.3d=TRUE,type="h",main=paste("RE window vs. recall - chroms",q),xlab=xname,ylab=yname,zlab="recall")
			ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
			if(ans_save == "y")
			{
				#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
				png_name <- paste("REwindowCompare_recall_chrom",q,"_",DandT2,".png",sep="")
				#png_name <- paste("REwindowCompare_recall_chrom",q,"_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
				#if(nchar(png_name) > 255)
				#{
				#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
				#	png_name <- substring(png_name,nm_length_rmv)
				#}
				png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 25)			
				scatterplot3d::scatterplot3d(rcll_chroms[[q]][,1:3],pch=16,highlight.3d=TRUE,type="h",main=paste("RE window vs. recall - chroms",q),xlab=xname,ylab=yname,zlab="recall")
				dev.off()
				save_flag <- 1
			}
			cat("\ncreating a plot for the f-measure data for chromosome",q,"\nplease wait...\n\n")
			while(dev.cur() != 1)
			{
				dev.off()
			}		
			scatterplot3d::scatterplot3d(fms_chroms[[q]][,1:3],pch=16,highlight.3d=TRUE,type="h",main=paste("RE window vs. f-measure - chroms",q),xlab=xname,ylab=yname,zlab="f-measure")
			ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
			if(ans_save == "y")
			{
				#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
				png_name <- paste("REwindowCompare_Fmeasure_chrom",q,"_",DandT2,".png",sep="")
				#png_name <- paste("REwindowCompare_Fmeasure_chrom",q,"_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
				#if(nchar(png_name) > 255)
				#{
				#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
				#	png_name <- substring(png_name,nm_length_rmv)
				#}
				png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 20)			
				scatterplot3d::scatterplot3d(fms_chroms[[q]][,1:3],pch=16,highlight.3d=TRUE,type="h",main=paste("RE window vs. f-measure - chroms",q),xlab=xname,ylab=yname,zlab="f-measure")
				dev.off()
				save_flag <- 1
			}
			cat("\ncreating a plot for the sum of intersections data for chromosome",q,"\nplease wait...\n\n")
			while(dev.cur() != 1)
			{
				dev.off()
			}		
			scatterplot3d::scatterplot3d(sum_ints_chroms[[q]][,1:3],pch=16,highlight.3d=TRUE,type="h",main=paste("RE window vs. sumOfIntersections - chroms",q),xlab=xname,ylab=yname,zlab="sum of intersections (bp)")
			ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
			if(ans_save == "y")
			{
				#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
				png_name <- paste("REwindowCompare_sumOfIntersections_chrom",q,"_",DandT2,".png",sep="")
				#png_name <- paste("REwindowCompare_sumOfIntersections_chrom",q,"_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
				#if(nchar(png_name) > 255)
				#{
				#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
				#	png_name <- substring(png_name,nm_length_rmv)
				#}
				png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 20)		
				scatterplot3d::scatterplot3d(sum_ints_chroms[[q]][,1:3],pch=16,highlight.3d=TRUE,type="h",main=paste("RE window vs. sumOfIntersections - chroms",q),xlab=xname,ylab=yname,zlab="sum of intersections (bp)")
				dev.off()
				save_flag <- 1
			}
		}
		
		#binding all chromosomes to one data frame
		prec_chroms_all <- rbind(prec_chroms_all,prec_chroms[[q]])
		rcll_chroms_all <- rbind(rcll_chroms_all,rcll_chroms[[q]])
		fms_chroms_all <- rbind(fms_chroms_all,fms_chroms[[q]])
		sum_ints_chroms_all <- rbind(sum_ints_chroms_all,sum_ints_chroms[[q]])	
	}
	
	#plotting and saving all chromosomes together:
	ans4 <- readline(prompt=cat("\nwould you like to plot all chromosomes together?\ny/n\n\n"))
	if(ans4 == "y")
	{	
		colors <- c("red","green","yellow","blue","black")

		prec_colors <- colors[as.numeric(factor(prec_chroms_all$chromosome))]	
		cat("\ncreating a plot for the precision data for all the chromosomes\nplease wait...\n\n")
		while(dev.cur() != 1)
		{
			dev.off()
		}	
		scatterplot3d::scatterplot3d(prec_chroms_all[,1:3],pch=16,color=prec_colors,type="h",main="RE window vs. precision - all chroms",xlab=xname,ylab=yname,zlab="precision")
		legend("bottom", legend = levels(factor(prec_chroms_all$chromosome)),col = colors, pch = 16, inset=-0.18, xpd = TRUE, horiz = TRUE)
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_precision_chromsAll_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_precision_chromsAll_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 25)		
			scatterplot3d::scatterplot3d(prec_chroms_all[,1:3],pch=16,color=prec_colors,type="h",main="RE window vs. precision - all chroms",xlab=xname,ylab=yname,zlab="precision")
			legend("bottom", legend = levels(factor(prec_chroms_all$chromosome)),col = colors, pch = 16, inset=-0.18, xpd = TRUE, horiz = TRUE)
			dev.off()
			save_flag <- 1
		}
		rcll_colors <- colors[as.numeric(factor(rcll_chroms_all$chromosome))]		
		cat("\ncreating a plot for the recall data for all the chromosomes\nplease wait...\n\n")
		while(dev.cur() != 1)
		{
			dev.off()
		}		
		scatterplot3d::scatterplot3d(rcll_chroms_all[,1:3],pch=16,color=rcll_colors,type="h",main="RE window vs. recall - all chroms",xlab=xname,ylab=yname,zlab="recall")	
		legend("bottom", legend = levels(factor(rcll_chroms_all$chromosome)),col = colors, pch = 16, inset=-0.18, xpd = TRUE, horiz = TRUE)
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_recall_chromsAll_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_recall_chromsAll_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 25)		
			scatterplot3d::scatterplot3d(rcll_chroms_all[,1:3],pch=16,color=rcll_colors,type="h",main="RE window vs. recall - all chroms",xlab=xname,ylab=yname,zlab="recall")	
			legend("bottom", legend = levels(factor(rcll_chroms_all$chromosome)),col = colors, pch = 16, inset=-0.18, xpd = TRUE, horiz = TRUE)
			dev.off()
			save_flag <- 1
		}
		fms_colors <- colors[as.numeric(factor(fms_chroms_all$chromosome))]	
		cat("\ncreating a plot for the f-measure data for all the chromosomes\nplease wait...\n\n")	
		while(dev.cur() != 1)
		{
			dev.off()
		}	
		scatterplot3d::scatterplot3d(fms_chroms_all[,1:3],pch=16,color=fms_colors,type="h",main="RE window vs. f-measure - all chroms",xlab=xname,ylab=yname,zlab="f-measure")		
		legend("bottom", legend = levels(factor(fms_chroms_all$chromosome)),col = colors, pch = 16, inset=-0.18, xpd = TRUE, horiz = TRUE)
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_Fmeasure_chromsAll_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_Fmeasure_chromsAll_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 20)		
			scatterplot3d::scatterplot3d(fms_chroms_all[,1:3],pch=16,color=fms_colors,type="h",main="RE window vs. f-measure - all chroms",xlab=xname,ylab=yname,zlab="f-measure")		
			legend("bottom", legend = levels(factor(fms_chroms_all$chromosome)),col = colors, pch = 16, inset=-0.18, xpd = TRUE, horiz = TRUE)
			dev.off()
			save_flag <- 1
		}
		sum_ints_colors <- colors[as.numeric(factor(sum_ints_chroms_all$chromosome))]		
		cat("\ncreating a plot for the sum of intersections data for all the chromosomes\nplease wait...\n\n")
		while(dev.cur() != 1)
		{
			dev.off()
		}		
		scatterplot3d::scatterplot3d(sum_ints_chroms_all[,1:3],pch=16,color=sum_ints_colors,type="h",main="RE window vs. sumOfIntersections - all chroms",xlab=xname,ylab=yname,zlab="sum of intersections (bp)")
		legend("bottom", legend = levels(factor(sum_ints_chroms_all$chromosome)),col = colors, pch = 16, inset=-0.18, xpd = TRUE, horiz = TRUE)
		ans_save <- readline("\nwould you like to save this plot?\ny/n\n\n")	
		if(ans_save == "y")
		{
			#testing to see that the name isn't longer than 256 characters (the file name can't be longer than PATH_MAX which is 256 bytes) if it is longer than we shorten it from the head of the name till it's short enough
			png_name <- paste("REwindowCompare_sumOfIntersections_chromsAll_",DandT2,".png",sep="")
			#png_name <- paste("REwindowCompare_sumOfIntersections_chromsAll_",sp1_1[1],"_vs_",sp1_2[1],"_",DandT2,".png",sep="")
			#if(nchar(png_name) > 255)
			#{
			#	nm_length_rmv <- as.integer(nchar(png_name)-255+1)
			#	png_name <- substring(png_name,nm_length_rmv)
			#}
			png(paste("~/Analyze4C/plots/",png_name,sep=""),width=1000,height=1000,pointsize = 20)		
			scatterplot3d::scatterplot3d(sum_ints_chroms_all[,1:3],pch=16,color=sum_ints_colors,type="h",main="RE window vs. sumOfIntersections - all chroms",xlab=xname,ylab=yname,zlab="sum of intersections (bp)")
			legend("bottom", legend = levels(factor(sum_ints_chroms_all$chromosome)),col = colors, pch = 16, inset=-0.18, xpd = TRUE, horiz = TRUE)
			dev.off()
			save_flag <- 1
		}
	}
	
	while(dev.cur() != 1)
	{
		dev.off()
	}

	#getting the max values for each category and printing to file (if the max value appears multiple times, then all examples will be printed to the file)
	prec_max <- c()
	rcll_max <- c()
	fms_max <- c()
	sum_ints_max <- c()
	prec_max <- rbind(prec_max,prec_all[which(as.numeric(prec_all$precision)==max(as.numeric(prec_all$precision))),])
	rcll_max <- rbind(rcll_max,rcll_all[which(as.numeric(rcll_all$recall)==max(as.numeric(rcll_all$recall))),])
	fms_max <- rbind(fms_max,fms_all[which(as.numeric(fms_all$Fmeasure)==max(as.numeric(fms_all$Fmeasure))),])
	sum_ints_max <- rbind(sum_ints_max,sum_ints_all[which(as.numeric(sum_ints_all$intersections)==max(as.numeric(sum_ints_all$intersections))),])
	prec_max <- rbind(prec_max,prec_trans[which(as.numeric(prec_trans$precision)==max(as.numeric(prec_trans$precision))),])
	rcll_max <- rbind(rcll_max,rcll_trans[which(as.numeric(rcll_trans$recall)==max(as.numeric(rcll_trans$recall))),])
	fms_max <- rbind(fms_max,fms_trans[which(as.numeric(fms_trans$Fmeasure)==max(as.numeric(fms_trans$Fmeasure))),])
	sum_ints_max <- rbind(sum_ints_max,sum_ints_trans[which(as.numeric(sum_ints_trans$intersections)==max(as.numeric(sum_ints_trans$intersections))),])
	for(q in 1:num_of_chroms)
	{
		prec_max <- rbind(prec_max,prec_chroms[[q]][which(as.numeric(prec_chroms[[q]]$precision)==max(as.numeric(prec_chroms[[q]]$precision))),])
		rcll_max <- rbind(rcll_max,rcll_chroms[[q]][which(as.numeric(rcll_chroms[[q]]$recall)==max(as.numeric(rcll_chroms[[q]]$recall))),])
		fms_max <- rbind(fms_max,fms_chroms[[q]][which(as.numeric(fms_chroms[[q]]$Fmeasure)==max(as.numeric(fms_chroms[[q]]$Fmeasure))),])
		sum_ints_max <- rbind(sum_ints_max,sum_ints_chroms[[q]][which(as.numeric(sum_ints_chroms[[q]]$intersections)==max(as.numeric(sum_ints_chroms[[q]]$intersections))),])
	}

	list_max <- list(prec_max,rcll_max,fms_max,sum_ints_max)
	
	sink(paste("~/Analyze4C/plots/REwindowCompare_maxValues_",DandT2,".txt",sep=""))
	print(list_max)
	sink()
	
	#saving the data into REwindowIntersections_plots.txt
	if(save_flag ==  1)
	{
		#saving the parameters and details to REwindowIntersections_plots.txt
		REwindowIntersections_plots[nrow(REwindowIntersections_plots)+1,] <- c(file.name_1,min_RE1,max_RE1,step_RE1,prp_1,RE_gap_1,bp_gap_1,co_type1,co1,file.name_2,min_RE2,max_RE2,step_RE2,prp_2,RE_gap_2,bp_gap_2,co_type2,co2,beta,DandT1)
		#sorting the list of experiments by bait alphabetically (and sorting the row indices)
		REwindowIntersections_plots <- REwindowIntersections_plots[order(REwindowIntersections_plots$filename_1),]
		rownames(REwindowIntersections_plots) <- seq(length=nrow(REwindowIntersections_plots))
		#adding the new data to the file (by erasing the old one and creating a new one)
		system("rm REwindowIntersections_plots.txt")
		write.table(REwindowIntersections_plots,"REwindowIntersections_plots.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)	
	}
	
	#erase all the files
	system("rm ~/Analyze4C/temp/*.bed")				
}
