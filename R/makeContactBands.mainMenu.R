#' @export

# i need to test this more, especially make sure that contactBands and contactBands_byChrom work

#note: if the folder of p-scores chosen is without bp, then there is a chance that the raw data was rearranged, and so we check that
#but if the folder is with bp then there is no chance that there was a rearrangement, since all the p-score files that their respective raw data files were rearranged
#had to have been calculated without bp

makeContactBands.mainMenu <- function(Experiments_4C,rearranged_rawData)
{
	ls_genomes <- system("ls ~/Analyze4C/genomes/",intern=TRUE)
	genome_file.name <- readline(prompt=cat("\nchoose the file with genome sizes that you would like to use:\n",ls_genomes,"\n",sep="\n"))
	genome_sizes <- read.table(paste("~/Analyze4C/genomes/",genome_file.name,sep=""))
	
	cat("\nin order to create contact bands you will first need to choose a p-score file to use,")
	cat("\ninput the max number of negative RE sites that won't constitut a gap,")
	cat("\ninput the max number of bp that are allowed to be in between RE sites and not be considered a gap,")
	cat("\nand finally input a cutoff that will differentiate between the positive and negative RE sites\n\n")
	
	{
		#getting the name of the p-score file the user wants to use
		choice <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the p-score file from:\n\n1) calculated using bps\n2) calculated without using bps\n\n")))
		repeat
		{
			if(choice == 1)
			{
				ls_files <- system("ls ~/Analyze4C/pScores/with_bp",intern=TRUE)
			}
			else
			{
				ls_files <- system("ls ~/Analyze4C/pScores/no_bp",intern=TRUE)
			}	
			file.name <- readline(prompt=cat("These are the p-score files available\nwhich would you like to use?\n",ls_files,"\n",sep="\n"))
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
		if(choice == 1)
		{
			data <- read.table(paste("~/Analyze4C/pScores/with_bp/",file.name,sep=""))
		}
		else
		{
			data <- read.table(paste("~/Analyze4C/pScores/no_bp/",file.name,sep=""))
		}
		
		#finding the specific p-score files information in Experiments_4C
		sp <- unlist(strsplit(file.name,"_"))
		out <- findIn.Experiments_4C(sp[1],sp[2],sp[3],Experiments_4C)
		
		#getting the cis chromosome number
		cis <- as.numeric(out[2])
	}
		
	#asking user if to create contact bands using a p-score cutoff for all chromosomes or a percent cutoff
	choice2 <- as.integer(readline(prompt=cat("\nenter the number of what type of contact band making you would like to do:\n\n1) all chromosomes with the same cutoff (in p-score)\n2) all chromosomes with the same cutoff percentage\n\n")))
	
	#getting the different paramaters and inputs from user
	RE_gap <- as.numeric(readline(prompt=cat("\nplease choose up to how many negative RE sites in between positive RE sites it is not considered a break in the contact band:\n\n")))
	bp_gap <- as.numeric(readline(prompt=cat("\nplease choose the max number of bp that are allowed to be in between RE sites and not be considered a gap:\n\n")))

	name_temp <- unlist(strsplit(file.name,"[.]"))[1]
	
	#creating the contact bands by asking for cutoffs first
	if(choice2 == 1)
	{
		co <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))
		cat("\ncalculating contact bands\nplease wait...\n\n")
		
		#checking if we are dealing with a file that has a translocation that isn't removed
		name_spl <- unlist(strsplit(name_temp,"_")) #splits the p-score file name. #note: i might not need this since we split the file name above
		#removing the translocated sections:
		if(("rearranged" %in% name_spl) & !("removed" %in% name_spl)) #if the p-score filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
		{
			#checking if the data had coverage removed, if yes then we remove the part of the name of the file that has to do with that (for technical reasons) and take the file from the correct folder
			if("covRemoved" %in% name_spl)
			{
				adrs <- "~/Analyze4C/rawData/coverage_removed/"
			}
			else
			{
				adrs <- "~/Analyze4C/rawData/rearranged/"
			}
			#finding which is the row in the 'rearranged_rawData' that contains information on the particular p-score data:
			
			rear_ind <- grep("rearranged",name_spl)
			col_num1 <- which(names(rearranged_rawData)=="Date_and_Time") #getting the number of column of the line numbers to remove
			dts <- rearranged_rawData[,col_num1] #gets all the dates and times in the file
			
			#getting the date information of the p-score file
			dte <- name_spl[rear_ind+1]
			
			#getting the time information of the p-score file, and converting it to look like the time in this format: hh:mm:ss
			tim <- name_spl[rear_ind+2]
			timspl <- strsplit(tim,"")
			timfin <- c(timspl[[1]][1])
			for(p in 2:length(timspl[[1]]))
			{
				if(p%%2!=0)
				{
					timfin <- paste(timfin,":",timspl[[1]][p],sep="")
				}
				else
				{
					timfin <- paste(timfin,timspl[[1]][p],sep="")
				}
			}

			#pasting the date and time together
			dandt <- paste(dte,timfin,sep=" ")
			
			#getting the row number of the p-score file in 'rearranged_rawData'
			row_num <- which(dts %in% dandt)
			
			#getting the number of column of the line numbers to remove
			col_num2 <- which(names(rearranged_rawData)=="Added_Sections_Lines")

			#using the line number we got from above we can extract the line numbers to remove
			row_num <- row_num[which(!is.na(rearranged_rawData[row_num,col_num2]))]
			lin_nums1 <- rearranged_rawData[row_num,col_num2]
			lin_nums2 <- as.integer(strsplit(lin_nums1,", ")[[1]])
			
			#removing the sections from the data
			#data <- cbind(data,rep(1,nrow(data)))
			#for(r in seq(1,length(lin_nums2),by=2))
			#{
			#	first <- lin_nums2[r]
			#	second <- lin_nums2[r+1]
			#	data[first:second,ncol(data)] <- 0			
			#}
			#data <- data[data[,ncol(data)]==1,1:3]
			#rownames(data) <- seq(length=nrow(data))
			
			#getting the raw data file that the p-score file was created from:
			ls_raws <- system(paste("ls",adrs),intern=TRUE)
			dAndtspl <- paste(name_spl[rear_ind+1],"_",name_spl[rear_ind+2],sep="")
			raws_inds <- grep(dAndtspl,ls_raws)
			for(o in raws_inds)
			{
				if(identical(grep("removed",ls_raws[o]),integer(0)))
				{
					raw_dat <- read.table(paste(adrs,ls_raws[o],sep=""))
					break
				}
			}

			#removing the sections from the raw data
			raw_dat <- cbind(raw_dat[,1],raw_dat[,2],rep(1,nrow(raw_dat)))
			for(r in seq(1,length(lin_nums2),by=2))
			{
				first <- lin_nums2[r]
				second <- lin_nums2[r+1]
				raw_dat[first:second,ncol(raw_dat)] <- 0			
			}
			raw_dat <- raw_dat[raw_dat[,ncol(raw_dat)]==1,1:2]
			rownames(raw_dat) <- seq(length=nrow(raw_dat))
			
			#adding 1 to each index of the raw data and p-score data, to create a bed file
			#there might be an issue if the index is the last one on the chromosome but i will leave it for now
			raw_dat <- cbind(raw_dat,raw_dat[,2]+1)
			data <- cbind(data[,1],data[,2],data[,2]+1,data[,3])
			
			#creating bed files
			write.table(raw_dat,"~/Analyze4C/temp/rawBed.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
			write.table(data,"~/Analyze4C/temp/psBed.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
			
			#intersecting the bed files so that we get the indices that aren't removed (since the way we removed was with line numbers of the files, we can't just remove them from the p-score data since we removed some lines from the cis when calculating the p-scores)
			system("bedtools intersect -a ~/Analyze4C/temp/psBed.bed -b ~/Analyze4C/temp/rawBed.bed > ~/Analyze4C/temp/out.bed")
			system("bedtools intersect -a ~/Analyze4C/temp/psBed.bed -b ~/Analyze4C/temp/rawBed.bed -c > ~/Analyze4C/temp/out2.bed") #this output will have all the same lines as psBed.bed but with an indication if the indices intersected or not (0 not intersected, 1 intersected)
			#there was a problem with intersected repeats, where if the same indices appeared in the translocated area that needs to be removed and an area that doesn'table
			#this will cause repeats in the intersections and some indices to be out of order
			#i couldn't figure out an automatic way to fix this
			#for now the user will fix it manually
			cat("\nview the 'out.bed' file for any repeats of intesections.\nthis happens if there is an index that appears in the translocated area and in areas that aren't.\n")
			cat("bedtools doesn't know how to distinguish between 'good' intersections like these and 'bad' ones.\nthese repeats should be at the beginning of the file (if the translocation occured on the left telomere of the chromosome.\n")
			cat("you could confirm that these are the lines that need to be erased by looking at the 'out2.bed' which shows that there are some indices that appear early on and in translocated areas that have been intersected.\n")
			cat("most chances are that you will see a '1' at the line of some indices with other indices that have a '0' surrounding them\n")
			cat("erase these lines\n")
			ent <- readline(prompt=cat("press ENTER once this is done in order to continue\n\n"))

			data <- read.table("~/Analyze4C/temp/out.bed")
			data <- cbind(data[,1],data[,2],data[,4])
			
			#removing the bed files created
			system("rm ~/Analyze4C/temp/psBed.bed")
			system("rm ~/Analyze4C/temp/rawBed.bed")
			system("rm ~/Analyze4C/temp/out.bed")
			system("rm ~/Analyze4C/temp/out2.bed")
			
			#getting the contact bands
			conts <- contactBands(data,RE_gap,bp_gap,cis,co,genome_sizes)
			
			#creating the contact bands file
			name <- paste(name_temp,"_",RE_gap,"maxREgap_",bp_gap,"maxbpgap_",co,"cutoff_contactSectionsRemoved_ContactBands.bed",sep="") #the added part of 'contactSectionsRemoved' to the name indicates	that the sections were removed and then the contact bands were found		
		}
		else #if the p-score file is of either a non translocation type or that the translocation was removed
		{
			conts <- contactBands(data,RE_gap,bp_gap,cis,co,genome_sizes)
			name <- paste(name_temp,"_",RE_gap,"maxREgap_",bp_gap,"maxbpgap_",co,"cutoff_ContactBands.bed",sep="")
		}
	}
	else
	{
		co_per <- as.numeric(readline(prompt=cat("\nplease enter a p-score cutoff percent (between 0 and 1) that will be applied to all chromosomes. Anything above the cutoff will be considered a positive RE site:\n\n")))
		cat("\ncalculating contact bands\nplease wait...\n\n")
		
		#checking if we are dealing with a file that has a translocation that isn't removed
		name_spl <- unlist(strsplit(file.name,"_")) #splits the p-score file name
		if(("rearranged" %in% name_spl) & !("removed" %in% name_spl)) #if the p-score filename contains the word 'rearranged' and not 'removed', meaning that it was rearranged and the sections that were rearranged still intact
		{
			conts <- contactBands_byChrom(data,RE_gap,bp_gap,co_per,genome_sizes,rearranged_rawData,flag=1,file.name=file.name) #the function 'contactBands_byChrom' with the input 1 will remove the added sections after applying the cutoff and then get contact bands
			conts <- conts[[1]]
			name <- paste(name_temp,"_",RE_gap,"maxREgap_",bp_gap,"maxbpgap_",co_per*100,"percentCutoff_contactSectionsRemoved_ContactBands.bed",sep="") #the added part of 'contactSectionsRemoved' to the name indicates that the after the cutoff was applied to the p-scores, then the sections were removed and contact bands were found
		}
		else
		{
			conts <- contactBands_byChrom(data,RE_gap,bp_gap,co_per,genome_sizes,rearranged_rawData) #the function 'contactBands_byChrom' will apply the cutoff without removing any sections
			conts <- conts[[1]]
			name <- paste(name_temp,"_",RE_gap,"maxREgap_",bp_gap,"maxbpgap_",co_per*100,"percentCutoff_ContactBands.bed",sep="")
		}
	}
	
	#asking user if to create a file
	ans <- readline(prompt=cat("\nwould you like to create a file of the contact bands?\ny/n\n\n"))
	if(ans == "y")
	{
		write.table(conts,paste("~/Analyze4C/contact_bands/",name,sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	}
	name #returns the name of the contact bands created
}
