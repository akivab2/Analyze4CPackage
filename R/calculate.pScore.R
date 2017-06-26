#' @export

#this function calculates the p score, either with a window size of only RE sites or also with bp
#the user needs to decide with which method they would like to calculate the p scores
#the user needs to enter the name of the raw data sgr file that they want to use
#the user then needs to enter the cis chromosome and the baits position
#finally the user needs to enter the number of RE sites per window and the number of bp per window if that option is chosen

calculate.pScore <- function(Experiments_4C)
{
	cat("\ncalculating p-scores:\n\n")
	inp <- 0 #gets a default of 0
	#getting the name of the raw data file the user wants to use
	repeat
	{
		choice <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the raw data file from to use for the p-score calculations:\n\n1) original\n2) rearranged\n3) coverage removed\n\n")))		
		if(choice == 2)
		{
			ls_files <- system("ls ~/Analyze4C/rawData/rearranged",intern=TRUE)
		}
		else if(choice == 3)
		{
			ls_files <- system("ls ~/Analyze4C/rawData/coverage_removed",intern=TRUE)
		}
		else
		{
			ls_files <- system("ls ~/Analyze4C/rawData/original",intern=TRUE)
		}
		
		if(length(ls_files) != 0)
		{	
			file.name <- readline(prompt=cat("\nThese are the raw data files available\nwhich would you like to use?\n",ls_files,"\n",sep="\n"))
			ind_files <- pmatch(file.name,ls_files,nomatch=0)
			if(ind_files == 0)
			{
				cat("no such file exists.\nplease try again.\n\n")
			}
			else
			{
				break
			}
		}
		else
		{
			cat("\nthis folder is empty. please choose a different folder from which you would like to use a file\n\n")
		}
	}
	#importing the data from the file
	if(choice == 2)
	{
		rawData <- read.table(paste("~/Analyze4C/rawData/rearranged/",file.name,sep=""))
	}
	else if(choice == 3)
	{
		rawData <- read.table(paste("~/Analyze4C/rawData/coverage_removed/",file.name,sep=""))
	}
	else
	{
		rawData <- read.table(paste("~/Analyze4C/rawData/original/",file.name,sep=""))
	}
	
	#finding the specific raw data files information in Experiments_4C
	sp1 <- unlist(strsplit(file.name,"[.]"))
	sp2 <- strsplit(sp1[[1]][1],"_")
	out <- findIn.Experiments_4C(sp2[[1]][1],sp2[[1]][2],sp2[[1]][3],Experiments_4C)
	
	#getting the cis chromosome number
	vp.chrom <- as.numeric(out[2])
	
	#getting the bait index
	vp.pos <- as.numeric(out[3])
	
	#getting the amount of bp around the bait that will be erased
	erase <- 1e6 #for now this will be the default. later i could have the user asked how many bp they want to be removed
	
	#ask the user if they would like to calculate the average RE per bp using RE.per.bp.R
	inp2 <- readline(prompt=cat("\nwould you like to calculate the average RE per bp?\ny/n\n"))
	if(inp2=="y")
	{
		RE.per.bp(file.name,choice)
	}
	
	#checking if the raw data taken has its indices in sequence (they might not be when rearranged, for translocations for example)
	#this is in order to see if the user can calculate p-scores using windows with bp size
	#note: this will only be able to identify a translocation that was added that creates a sequence of indices out of order
	#but if the indices added are smaller than the ones they are added to, we will no detect this
	#in the future i could change the detection function to something else that will be able to find this (maybe by looking at the name of the file or by looking at the 'rearranged_rawData.txt' file)
	num_of_chr <- length(unique(rawData[,1])); #getting the number of chromosomes
	for(k in 1:num_of_chr)
	{
		if(is.unsorted(rawData[rawData[,1]==k,2]))
		{
			cat("\nthe raw data seems to be rearranged.\nyou may not calculate p-scores with window sizes by bp, only by RE sites.\n\n")
			inp <- 1
			break
		}
	}

	#ask if the user would like to calculate the p score with bp window size or just RE site window size
	#1 is for only RE site windows, 2 is for both RE site and bp window sizes
	if(inp == 0) #if it got a value of 1 we know the raw datas indices are not in sequence
	{
		inp <- as.integer(readline(prompt="\nif you would like to calculate the p score using only RE site window sizes press 1.\nif you would like to calculate the p score using RE site window sizes and bp window sizes press 2.\n\n"))
	}
		
	#get the window size in RE site
	REs <- as.integer(readline(prompt=cat("\nplease enter the size of window by RE sites (should be an even number):\n\n")))

	if(inp==1)
	{
		#using the paramaters above we calculate the p score with RE site window sizes only
		cat("\ncalculating the p-scores and creating the p-score file\nplease wait...\n\n")		
		ps <- pScore_nokb(rawData,vp.chrom,vp.pos,erase,REs/2)
		
		#creating the p scores file name
		ps.name <- paste(sp1[[1]][1],"_",REs,"RE_","pScore.sgr",sep="")

		#write p score to file
		write.table(ps,paste("~/Analyze4C/pScores/no_bp/",ps.name,sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	}
	else if(inp==2)
	{
		#get the window size in bp
		wind <- readline(prompt=cat("\nplease enter the size of window by bp (note that the size should probably be proportionate to the RE site window size):\n\n"))
		wind <- as.integer(wind)
		
		#using the paramaters above we calculate the p score with RE site and bp window sizes
		cat("\ncalculating the p-scores and creating the p-score file\nplease wait...\n\n")		
		ps <- pScore(rawData,vp.chrom,vp.pos,erase,wind,REs/2)
		
		#creating the p scores file name
		ps.name <- paste(sp1[[1]][1],"_",REs,"RE_",wind,"bp_pScore.sgr",sep="")
		
		#write p score to file
		write.table(ps,paste("~/Analyze4C/pScores/with_bp/",ps.name,sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	}
}
