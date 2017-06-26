#' @export

#this function rearranges and/or removes sections of the raw data and/or p-score files.
#the rearrangement is necesarry when for example a translocation occured and the way the genome is organized in reality is different than the known genome assembly.
#1) the function first requests from the user to choose a raw data to rearrange.
#2) then the user is asked to enter the number of sections the genome will be divided into in order to be rearranged, for this the user should view the raw data beforehand
#and know what line numbers of the file to enter that will divide the genome into sections, the amount of pairs of line numbers (first and last) will determine how many sections there will be.
#3) the user will then enter the line numbers and the sections will be arranged in the order of the line numbers the user input.
#4) the user can then choose to change the chromosome numbers of sections, this is mainly meant for the purpose of calculating p-scores (since the algorithm requires it).
#for this the user should view the the 'temp' file that is created from the rearranged data and decide if and where the chromosome numbers should be changed.
#the 'temp' file is then removed.
#5) asking user if they would like to calculate the p-score with the new rearranged data. The method of calculation uses only RE site windows (since the order of indices is rearranged).
#once the p-scores are calculated the user could decide if to remove the rearranged sections from the p-score file, this is in order to be able to view the data in a browser.
#the line numbers of the of where to remove could be found if the p-score file is viewed first in 'pScores/no_bp/'. A new p-score file with the removed sections is then created and the other
#removed. The names of the p-score files have the date and time in them and say if they are of a translocation rearrangement. The data about the file is entered in the file 'pScores_of_rearranged'
#which contains all the rearranged p-score files and information about each (the user inputs some of these), and should also write if the file was erased.
#6) creating the rearranged raw data file, removing the rearranged raw data sections as requested by user (the data could and should be first viewed in '/rawData/rearranged/').
#once the sections are removed the file of the unremoved data could be removed and the removed data file is created. All this information about the files is stored in the file 'rearranged_rawData',
#including the date and time of creation and if it was a translocation (this data is written in the file name as well)

#in short the user could rearrange raw data files, remove sections from it, and calculate p-scores using these files

rawData_rearranger <- function(Experiments_4C,pScores_of_rearranged,rearranged_rawData)
{
	#introduction
	cat("\n      ****raw data rearranger****      \n")
	cat("\nthis function rearranges the raw data file per request")
	cat("\nexamples of times when this is needed include when the raw data doesn't comply to the properties of the genome browser or when a translocation occured")
	cat("\n\nin order to do this there are a few steps necessary:\n\n")
	cat("1) choose a raw data file you would like to rearrange\n")
	cat("2) choose the amount of sections the raw data should be divided into\n")
	cat("3) inputting the line numbers of each section in the order that they will be combined in (the function will then combine the sections)\n")
	cat("4) deciding if to leave the original chromosome numbers or change those as well\n")
	cat("5) deciding if to calculate p scores with the rearranged data, if yes then also deciding if to remove sections of p-score file\n")
	cat("6) deciding if to remove certain sections of the raw data\n")
	
	#removing sections without rearranging

	#1)
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
			file.name <- readline(prompt=cat("These are the raw data files available\nwhich would you like to use?\n",ls_files,"\n",sep="\n"))
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
	}
	
	cat("\nfor the next steps it is advised to manually view the raw data original file in a genome browser and take note of the actual line numbers where you would like the data be splitted\n")
	cat("then for the amount of sections just add 1 to the amount of line numbers you took note of\n")
	cat("note that the amount of sections refers to the whole genome (meaning, take into account that the whole genome file is considered one block at first and will be divided into sections as requested by user)\n")
	
	#2)
	secs <- readline(prompt="\ninto how many sections would you like the raw data be divided?\n")
	secs <- as.integer(secs)
	
	#3)
	flag <- 0
	cat("\nnote: the sections will be arranged in the order of the entered line numbers\n\n")
	len <- nrow(rawData)
	while(flag == 0)
	{
		new <- c()
		for(i in 1:secs)
		{
			cat("\nsection number ",i," :\n",sep="")
			first1 <- as.integer(readline(prompt="\nenter the first line number of the section:\n"))
			second1 <- as.integer(readline(prompt="\nenter the second line number of the section:\n"))
				
			new <- rbind(new,rawData[first1:second1,])	
		}
		rownames(new) <- seq(length=nrow(new))

		#showing the sizes of the original and rearranged data, if no data was meant to be removed they should be the same
		cat("\n\noriginal data size: ",len,"\n\nnew rearranged data size: ",nrow(new),"\n",sep="")
		ans0 <- readline(prompt="\nwould you like to retry entering the sections?\ny/n\n")
		if(ans0 != "y")
		{
			flag <- 1
		}
	}
	
	#4)
	#creating a temp file to view and see if and which chromosome numbers to change and which line numbers to input for this
	write.table(new,"~/Analyze4C/rawData/rearranged/temp.sgr",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

	#creating a vector ('added_secs') of indices that were of areas that were added to the chromosomes
	cat("\nview the 'temp file in 'rearranged' folder, and enter the range of line numbers that contain the areas that were added to chromosomes\n")
	cat("meaning the sections that won't be compatible if the data is viewed in a browser (the sections that don't appear as in the original genome sequence and were added in 'impossible' places)\n\n")
	added_secs <- c()
	ans5 <- readline(prompt=cat("\nare there any sections of these type?\ny/n\n\n"))
	while(ans5 == "y")
	{
		first5 <- as.integer(readline(prompt="\nenter the first line number of the section:\n"))
		second5 <- as.integer(readline(prompt="\nenter the second line number of the section:\n"))
		added_secs <- c(added_secs,first5,second5)
		ans5 <- readline(prompt="\nare there anymore sections?\ny/n\n\n")
	}

	#changing chromosome numbers by either entering them manually or by using the indices entered earlier as added sections to chromosomes
	cat("\n\nview the 'temp' file in the 'rearranged' folder and decide if you would like to change any chromosome numbers. If yes enter the line numbers as they appear in the 'temp' file.\n")
	ans1 <- readline(prompt="\npress 'y' if you would like to change the chromosome numbers of any sections of the genome\ny/n\n\n")
	if(ans1 == "y")
	{
		ans15 <- readline(prompt=cat("\nshould we use the sections entered as 'added to the genome in 'impossible' places'?\ny/n\n"))
		if(ans15=="y")
		{
			for(r in seq(1,length(added_secs),by=2))
			{
				first2 <- added_secs[r]
				second2 <- added_secs[r+1]
				new_chr <- as.integer(readline(prompt=cat("\nthe lines are ",first2,":",second2,"\nenter the new chromosome number:\n",sep="")))
				new[first2:second2,1] <- new_chr				
			}
		}
		else
		{
			temp1 <- "y"
			while(temp1 == "y")
			{
				first2 <- as.integer(readline(prompt="\nenter the first line number of the section:\n"))
				second2 <- as.integer(readline(prompt="\nenter the second line number of the section:\n"))
				new_chr <- as.integer(readline(prompt="\nenter the new chromosome number:\n"))
				new[first2:second2,1] <- new_chr
				temp1 <- readline(prompt="\nwould you like to change any other chromosome numbers?\ny/n\n\n")
			}	
		}
	}
	#removing the 'temp' file because it is not needed anymore
	system("rm ~/Analyze4C/rawData/rearranged/temp.sgr")
	
	#getting the date and time to add to the names of the files and the 'pScores_of_rearranged' and 'rearranged_rawData' files
	DandT1 <- toString(Sys.time())
	DandT2 <- gsub(" ","_",DandT1)
	DandT2 <- gsub(":","",DandT2)
	
	#calculate p scores (no kb)
	#5)
	ans2 <- readline(prompt=cat("\nwould you like to calculate p-scores with the rearranged data?\n(note: since the index numbers are all rearranged the only method available for p-score calculation is only with RE site windows)\ny/n\n\n"))
	if(ans2 == "y")
	{
		#getting the amount of bp around the bait that will be erased
		erase <- 1e6 #for now this will be the default. later i could have the user asked how many bp they want to be removed

		#ask the user if they would like to calculate the average RE per bp using RE.per.bp.R
		ans3 <- readline(prompt=cat("\nwould you like to calculate the average RE per bp? (note: in order to check this you must use the original raw data file)\ny/n\n"))
		if(ans3=="y")
		{
			RE.per.bp()
		}
		
		#get the window size in RE site
		REs <- as.integer(readline(prompt=cat("\nplease enter the size of window by RE sites (should be an even number):\n\n")))
	
		#telling the user that a p-score file is being made
		cat("\na p-score file is being created.\nplease wait...\n\n")

		#using the paramaters above we calculate the p score with RE site window sizes only
		ps <- pScore_nokb(new,vp.chrom,vp.pos,erase,REs/2)
		
		#adding this data to the file pScores_of_rearranged, which includes the name of the experiment, the date and time when file was made, and a description of the rearrangement
		des <- readline(prompt="\nenter a short description of the p-score file:\n\n")
		ind_temp1 <- nrow(pScores_of_rearranged)+1
		note <- NA
		pScores_of_rearranged[nrow(pScores_of_rearranged)+1,] <- c(sp1[[1]][1],DandT1,des,note)
		#sorting the list by experiment name alphabetically (and sorting the row indices)
		pScores_of_rearranged <- pScores_of_rearranged[order(pScores_of_rearranged$Experiment),]
		rownames(pScores_of_rearranged) <- seq(length=nrow(pScores_of_rearranged))
		#adding the new data to the file (by erasing the old one and creating a new one)
		system("rm pScores_of_rearranged.txt")
		write.table(pScores_of_rearranged,"pScores_of_rearranged.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
		
		#creating the p scores file name
		ans6 <- readline(prompt="\nwas this rearrangement a switching of a translocation?\ny/n\n")
		if(ans6 == "y")
		{
			ps.name <- paste(sp1[[1]][1],"_translocatedSwitched_rearranged_",DandT2,"_",REs,"RE_pScore.sgr",sep="")
		}
		else
		{
			ps.name <- paste(sp1[[1]][1],"_rearranged_",DandT2,"_",REs,"RE_pScore.sgr",sep="")
		}
		
		#creating the p-score file
		write.table(ps,paste("~/Analyze4C/pScores/no_bp/",ps.name,sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
		
		ans7 <- readline(prompt=cat("\nwould you like to remove the rearranged sections?\n(this is mostly advised when the rearrangement 'messed up' the datas indices so it can't be viewed in a genome browser, for example when a translocation has been switched)\ny/n\n"))
		if(ans7 == "y")
		{
			len2 <- nrow(ps)
			ps_temp <- cbind(ps,rep(1,len2))
			ans17 <- readline(prompt=cat("\nshould we use the sections entered as 'added to the genome in 'impossible' places' as the sections to be removed?\nnote: this would be advised only if there were no removed sections from the cis in the p-score calculation)\ny/n\n"))
			if(ans17=="y")
			{
				for(r in seq(1,length(added_secs),by=2))
				{
					first3 <- added_secs[r]
					second3 <- added_secs[r+1]
					ps_temp[first3:second3,ncol(ps_temp)] <- 0			
				}
			}
			else
			{
				temp <- 1
				while(temp == 1)
				{
					cat("\nenter the line number of the section that you would like be removed (best to view the file first in the 'pScores/no_bp/' folder)\n\n")
					first3 <- as.integer(readline(prompt="\nenter the first line number of the section:\n"))
					second3 <- as.integer(readline(prompt="\nenter the second line number of the section:\n"))
					ps_temp[first3:second3,ncol(ps_temp)] <- 0
					ans8 <- readline(prompt="\nwould you like to remove another section?\ny/n\n")
					if(ans8 != "y")
					{
						temp <- 0
					}
				}
			}	
			ps2 <- ps_temp[ps_temp[,ncol(ps_temp)]==1,1:3]
			rownames(ps2) <- seq(length=nrow(ps2))

			#asking user if they would like to create a file with the p-score data after sections were removed
			ans8 <- readline(prompt=cat("\nwould you like to create a file with this new p-score data after the sections were removed?\ny/n\n\n"))
			if(ans8 == "y")
			{					
				#adding this data to the file pScores_of_rearranged, which includes the name of the experiment, the date and time when file was made, and a description of the rearrangement
				des <- readline(prompt="\nenter a short description of the p-score file:\n\n")
				note <- NA
				pScores_of_rearranged[nrow(pScores_of_rearranged)+1,] <- c(sp1[[1]][1],DandT1,des,note)
				#sorting the list by experiment name alphabetically (and sorting the row indices)
				pScores_of_rearranged <- pScores_of_rearranged[order(pScores_of_rearranged$Experiment),]
				rownames(pScores_of_rearranged) <- seq(length=nrow(pScores_of_rearranged))
				#adding the new data to the file (by erasing the old one and creating a new one)
				system("rm pScores_of_rearranged.txt")
				write.table(pScores_of_rearranged,"pScores_of_rearranged.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
				
				#creating the  p scores file name
				if(ans6 == "y")
				{
					ps.name2 <- paste(sp1[[1]][1],"_translocatedSwitched_removed_fromPS_rearranged_",DandT2,"_",REs,"RE_pScore.sgr",sep="")
				}
				else
				{
					ps.name2 <- paste(sp1[[1]][1],"_sections_removed_fromPS_rearranged_",DandT2,"_",REs,"RE_pScore.sgr",sep="")
				}
				
				#creating the p-score file
				write.table(ps2,paste("~/Analyze4C/pScores/no_bp/",ps.name2,sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
				
				ans10 <- readline(prompt="\nwould you like to erase the p-score file before sections were removed?\ny/n\n")
				if(ans10 == "y")
				{
					system(paste("rm ~/Analyze4C/pScores/no_bp/",ps.name,sep=""))
					pScores_of_rearranged[ind_temp1,]$Note <- "erased"
					system("rm pScores_of_rearranged.txt")
					write.table(pScores_of_rearranged,"pScores_of_rearranged.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
				}
			}			
		}
	}
	
	#creating the rearranged file, removing sections that would mess up the browser, and creating a file of after removal
	#6)
	
	#telling user that a raw data file of the rearranged data is going to be created
	cat("\ncreating raw data of rearranged raw data:\n\n")
	
	#adding this data to the file rearranged_rawData, which includes the name of the experiment, the date and time when file was made, and a description of the rearrangement
	des <- readline(prompt="\nenter a short description of the rearranged raw data file:\n\n")
	ind_temp2 <- nrow(rearranged_rawData)+1
	added_secs2 <- toString(added_secs)
	note <- NA
	rearranged_rawData[nrow(rearranged_rawData)+1,] <- c(sp1[[1]][1],DandT1,des,added_secs2,note)
	#sorting the list by experiment name alphabetically (and sorting the row indices)
	rearranged_rawData <- rearranged_rawData[order(rearranged_rawData$Experiment),]
	rownames(rearranged_rawData) <- seq(length=nrow(rearranged_rawData))
	#adding the new data to the file (by erasing the old one and creating a new one)
	system("rm rearranged_rawData.txt")
	write.table(rearranged_rawData,"rearranged_rawData.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
	
	#creating the new file name
	ans11 <- readline(prompt="\nwas this rearrangement a switching of a translocation?\ny/n\n")
	if(ans11 == "y")
	{
		ra.name <- paste(sp1[[1]][1],"_translocatedSwitched_rearranged_",DandT2,".sgr",sep="")
	}
	else
	{
		ra.name <- paste(sp1[[1]][1],"_rearranged_",DandT2,".sgr",sep="")
	}
	
	#telling the user a file with the rearranged raw data is being created
	cat("\na file with the rearranged raw data is being created.\nplease wait...\n\n")
	
	#creating the file
	write.table(new,paste("~/Analyze4C/rawData/rearranged/",ra.name,sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

	#ask if the user would like to remove sections of the new rearranged data
	ans4 <- readline(prompt=cat("\nwould you like to remove the rearranged sections?\n(this is mostly advised when the rearrangement 'messed up' the datas indices so it can't be viewed in a genome browser, for example when a translocation has been switched)\ny/n\n"))
	if(ans4 == "y")
	{
		len3 <- nrow(new)
		new_temp <- cbind(new,rep(1,len3))
		ans16 <- readline(prompt=cat("\nshould we use the sections entered as 'added to the genome in 'impossible' places'?\ny/n\n"))
		if(ans16=="y")
		{
			for(r in seq(1,length(added_secs),by=2))
			{				
				first4 <- added_secs[r]
				second4 <- added_secs[r+1]
				new_temp[first4:second4,ncol(new_temp)] <- 0				
			}		
		}
		else	
		{
			temp2 <- 1
			while(temp2 == 1)
			{
				cat("\nenter the line numbers of the section that you would like be removed (best to the view the file that was already created first in the folder '/rawData/rearranged/')\n\n")
				first4 <- as.integer(readline(prompt="\nenter the first line number of the section:\n"))
				second4 <- as.integer(readline(prompt="\nenter the second line number of the section:\n"))
				new_temp[first4:second4,ncol(new_temp)] <- 0
				ans12 <- readline(prompt="\nwould you like to remove another section?\ny/n\n")
				if(ans12 != "y")
				{
					temp2 <- 0
				}
			}
		}	
		new2 <- new_temp[new_temp[,ncol(new_temp)]==1,1:3]
		rownames(new2) <- seq(length=nrow(new2))
		
		#asking user if they would like to create a file with the rearranged data after sections were removed
		ans13 <- readline(prompt=cat("\nwould you like to create a file with the rearranged data after the sections were removed?\ny/n\n\n"))
		if(ans13 == "y")
		{				
			#adding this data to the file rearranged_rawData, which includes the name of the experiment, the date and time when file was made, and a description of the rearrangement
			des <- readline(prompt="\nenter a short description of the rearranged raw data file:\n\n")
			added_secs2 <- NA
			note <- NA
			rearranged_rawData[nrow(rearranged_rawData)+1,] <- c(sp1[[1]][1],DandT1,des,added_secs2,note)
			#sorting the list by experiment name alphabetically (and sorting the row indices)
			rearranged_rawData <- rearranged_rawData[order(rearranged_rawData$Experiment),]
			rownames(rearranged_rawData) <- seq(length=nrow(rearranged_rawData))
			#adding the new data to the file (by erasing the old one and creating a new one)
			system("rm rearranged_rawData.txt")
			write.table(rearranged_rawData,"rearranged_rawData.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
			
			#creating the new file name
			if(ans11 == "y")
			{
				ra.name2 <- paste(sp1[[1]][1],"_translocatedSwitched_removed_fromRD_rearranged_",DandT2,".sgr",sep="")
			}
			else
			{
				ra.name2 <- paste(sp1[[1]][1],"_sections_removed_fromRD_rearranged_",DandT2,".sgr",sep="")
			}
			
			#creating the file
			write.table(new2,paste("~/Analyze4C/rawData/rearranged/",ra.name2,sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
			
			#there was a possibility to have an option of erasing the raw data file with the translocation rearranged
			#this probably shouldn't be an option since the raw data file with translocation still in it is useful for other things (like calculating p-scores when comparing to rna-seq)
			#this is what i could use if i want to re-add it:
			
			#ans14 <- readline(prompt="\nwould you like to erase the raw data file before sections were removed?\ny/n\n")
			#if(ans14 == "y")
			#{
			#	system(paste("rm ~/Analyze4C/rawData/rearranged/",ra.name,sep=""))
			#	rearranged_rawData[ind_temp2,]$Note <- "erased"
			#	system("rm rearranged_rawData.txt")
			#	write.table(rearranged_rawData,"rearranged_rawData.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
			#}
		}		
	}
}
