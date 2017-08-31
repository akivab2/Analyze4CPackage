#' @export

coverage.changer.main_menu <- function(Experiments_4C,coverageChanged_Experiments)
{
	if((length(system("ls ~/Analyze4C/rawData/original",intern=TRUE)) == 0) & (length(system("ls ~/Analyze4C/rawData/rearranged",intern=TRUE)) == 0) & (length(system("ls ~/Analyze4C/rawData/coverage_removed",intern=TRUE)) == 0))
	{
		cat("\nThere are no raw data files available\n\n")
	}
	else
	{
		#getting the name of the raw data file the user wants to use
		repeat
		{
			choice <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the raw data file from:\n\n1) original\n2) rearranged\n3) coverage removed\n\n")))		
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
		
		#get the cis information from Experiments_4C using the name of the file chosen
		#finding the specific raw data files information in Experiments_4C
		sp1 <- unlist(strsplit(file.name,"[.]"))
		sp2 <- strsplit(sp1[[1]][1],"_")
		out <- findIn.Experiments_4C(sp2[[1]][1],sp2[[1]][2],sp2[[1]][3],Experiments_4C)
		
		#getting the cis chromosome number
		cis <- as.numeric(out[2])

		#ask the user if they would like to see the data stats (this is done before the coverage removal in order to know how much to remove if at all)
		ans1 <- readline(prompt=cat("\nwould you like to view the data stats?\ny/n\n\n"))
		if(ans1 == "y")
		{
			rawData.stats(Experiments_4C,rawData,file.name)		
		}
		
		#choosing the type of coverage removal
		ans4 <- as.integer(readline(prompt=cat("\nchoose the type of coverage removal you would like to do:\n1) by coverage percent\n2) by number of reads\n\n")))
		if(ans4==1)
		{
			type <- "by percentage"
			flag1 <- 0
			while(flag1==0)
			{				
				#get the coverage and reads stats for the original data
				ans5 <- readline(prompt=cat("\nwould you like to get coverage stats for the original data?\ny/n\n\n"))
				if(ans5=="y")
				{
					#getting the coverage stats
					rawData.stats(Experiments_4C,rawData,file.name)	
				}
				
				#getting the percentage of coverage to be removed
				rmv <- as.numeric(readline(prompt=cat("\nenter the percentage of coverage to be removed (between 0 and 1):\n\n")))
				
				#activate coverageChanger_chrChooser get the output into a variable (a list) and separate the list into the 2 outputs: the new data, and from what chromosomes there was a removal
				out2 <- coverageChanger_chrChooser(rawData,rmv,cis)
				newData <- out2[[1]]
				location_text <- out2[[2]]
				
				#ask the user if they would like to see the data stats for the original data and the removed coverage data, in order to compare them
				ans6 <- readline(prompt=cat("\nwould you like to view the data stats?\ny/n\n\n"))
				if(ans6 == "y")
				{
					repeat
					{
						ans7 <- as.integer(readline(prompt=cat("\nchoose the data you would like to see the stats for:\n\n1) original data\n2) removed coverage data\n3) exit\n\n")))
						if(ans7 == 1)
						{
							rawData.stats(Experiments_4C,rawData,file.name)
						}
						else if(ans7 == 2)
						{
							rawData.stats(Experiments_4C,newData,file.name)
						}
						else
						{
							break
						}	
					}			
				}
				
				ans8 <- readline(prompt=cat("\nare you satisfied with this coverage removal?\ny/n\n\n"))
				if(ans8=="y")
				{
					flag1 <- 1
				}
			}	
		}
		else
		{
			type <- "by reads"		
			out2 <- remove_coverage_byReads(rawData,Experiments_4C,file.name,cis)
			newData <- out2[[1]]
			location_text <- out2[[2]]			
			min_reads <- out2[[3]]
		}
		
		#ask user if they would like to create a file with the data
		ans2 <- readline(prompt=cat("\nwould you like to create a file with this new data?\ny/n\n\n"))
		if(ans2 == "y")
		{	
			#getting the date and time of coverage removal
			DandT1 <- toString(Sys.time())
			
			#adding this data to the file coverageChanged_Experiments, which includes the name of the experiment, the date and time of coverage removal, and the percentage removed and from which chromosomes
			if(ans4==1) #if removed coverage by percentage
			{
				rmv2 <- paste((rmv*100),"% removed",sep="")
				min_reads <- "NA"
			}
			else #if removed coverage by number of reads
			{
				rmv2 <- "NA"
			}	
			coverageChanged_Experiments[nrow(coverageChanged_Experiments)+1,] <- c(sp1[[1]][1],DandT1,type,location_text,rmv2,min_reads)
			#sorting the list of experiments by bait alphabetically (and sorting the row indices)
			coverageChanged_Experiments <- coverageChanged_Experiments[order(coverageChanged_Experiments$Experiment),]
			rownames(coverageChanged_Experiments) <- seq(length=nrow(coverageChanged_Experiments))
			#adding the new data to the file (by erasing the old one and creating a new one)
			system("rm coverageChanged_Experiments.txt")
			write.table(coverageChanged_Experiments,"coverageChanged_Experiments.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
			
			#creating the new file name
			DandT2 <- gsub(" ","_",DandT1)
			DandT2 <- gsub(":","",DandT2)
			new.name <- paste(sp1[[1]][1],"_covRemoved_",DandT2,".sgr",sep="")
			
			#creating the file
			write.table(newData,paste("~/Analyze4C/rawData/coverage_removed/",new.name,sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
		}
	}	
}
