#' @export

#the function gets two raw data files and combines them by adding the reads
#this function works only with files from 'normal or the 'removed coverage' folders (not 'rearranged')

combineLanes <- function(Experiments_4C)
{
	#choosing the 1st raw data file to work with:
	#getting the name of the raw data file the user wants to use
	repeat
	{
		choice1 <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the 1st raw data file from:\n\n1) original\n2) coverage removed\n\n")))
		if(choice1 == 2)
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
			break
		}
	}
	#importing the data from the file
	if(choice1 == 2)
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
	
	#choosing the 2nd raw data file to work with:
	#getting the name of the raw data file the user wants to use
	repeat
	{
		choice2 <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the 1st raw data file from:\n\n1) original\n2) coverage removed\n\n")))
		if(choice2 == 2)
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
			break
		}
	}
	#importing the data from the file
	if(choice2 == 2)
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
	
	#make sure the files have the same indices and come from the same bait
	if(sp2_1[[1]][1]==sp2_2[[1]][1] & sp2_1[[1]][2]==sp2_2[[1]][2])
	{
		if(all(rawData_1[,1]==rawData_2[,1]) & all(rawData_1[,2]==rawData_2[,2]))
		{
			#getting the time and date
			DandT1 <- toString(Sys.time())
			DandT2 <- gsub(" ","_",DandT1)
			DandT2 <- gsub(":","",DandT2)
			
			#creating the new raw data by adding the first and second files reads
			new_data <- cbind(rawData_1[,1],rawData_1[,2],(rawData_1[,3]+rawData_2[,3]))
			
			combinedLane <- paste(sp2_1[[1]][3],sp2_2[[1]][3],sep="")
			#creating the correct name for the new file and placing it in the correct folder
			if(choice1 == 2)
			{
				nm <- paste("~/Analyze4C/rawData/original/",sp2_1[[1]][1],"_",sp2_1[[1]][2],"_",combinedLane,sep="")
				for(j in 4:length(sp2_1[[1]]))
				{
					nm <- paste(nm,"_",sp2_1[[1]][j],sep="")
				} 
			}
			else if(choice2 == 2)
			{
				nm <- paste("~/Analyze4C/rawData/coverage_removed/",sp2_2[[1]][1],"_",sp2_2[[1]][2],"_",combinedLane,sep="")
				for(j in 4:length(sp2_2[[1]]))
				{
					nm <- paste(nm,"_",sp2_2[[1]][j],sep="")
				}
			}
			else
			{
				nm <- paste("~/Analyze4C/rawData/original/",sp2_1[[1]][1],"_",sp2_1[[1]][2],"_",combinedLane,sep="")
				for(j in 4:length(sp2_1[[1]]))
				{
					nm <- paste(nm,"_",sp2_1[[1]][j],sep="")
				} 
			}
			nm <- paste(nm,"_",DandT2,".sgr",sep="")
			
			write.table(new_data,nm,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

			#adding the new data to the 'Experiments_4C'
			Experiments_4C[nrow(Experiments_4C)+1,] <- c(sp2_1[[1]][1],sp2_1[[1]][2],combinedLane,out_1[[2]],out_1[[3]],out_1[[4]],out_1[[5]])
			#sorting the list of experiments by bait alphabetically (and sorting the row indices)
			Experiments_4C <- Experiments_4C[order(Experiments_4C$Bait),]
			rownames(Experiments_4C) <- seq(length=nrow(Experiments_4C))
			#adding the new data to the file (by erasing the old one and creating a new one)
			system("rm Experiments_4C.txt")
			write.table(Experiments_4C,"Experiments_4C.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
			
			combined_files <- read.table("combined_files.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			#add to 'combined_files' which 2 files were combined 
			combined_files[nrow(combined_files)+1,] <- c(file.name_1,file.name_2,DandT1)
			#sorting the list of experiments by raw1 alphabetically (and sorting the row indices)
			combined_files <- combined_files[order(combined_files$Raw1),]
			rownames(combined_files) <- seq(length=nrow(combined_files))
			#adding the new data to the file (by erasing the old one and creating a new one)
			system("rm combined_files.txt")
			write.table(combined_files,"combined_files.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)			
		}
	}
}
