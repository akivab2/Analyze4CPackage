#' @export

#this function gets p-score data and a cutoff percentage and gets the contact bands of the p-scores that are above the value of the p-score at the cutoff percentile
#this is done per chromosome - meaning that the percentile is calculated for each chromosome separately
#if the p-score data is taken from a rearranged genome then we first calculate the cutoff percentile, then remove those sections that were added, and then find the contact bands (for each chromosome separately of course)
#the reason that this is the order is because the p-scores are calculated with these sections, and so when we look at the percentiles we need to consider all the p-scores that exist when calculated
#once we calculate the p-score cutoff we can then dispose of these added sections and find the contact bands (we can't find them with these sections since the idnices will be all messed up and out of order, os it is best to remove them) 

contactBands_byChrom <- function(dat,RE_gap,bp_gap,co_per,genome_sizes,rearranged_rawData,flag_rmv=0,in_lns=0,flag=0,file.name="0")
{
	#ln <- length(unique(dat[,1]));
	rmv_lst <- list()
	ind_counter <- 0
	
	#contains all the contacts
	conts <- c()
	
	#removing the added sections:
	if(flag==1 & file.name!="0")
	{	
		#getting the raw data file that the p-score file was created from:
		
		#splits the p-score file name
		name_spl <- unlist(strsplit(file.name,"[.]"))
		name_spl <- unlist(strsplit(name_spl[1],"_")) 
		
		#checking if the data had coverage removed, if yes then we remove the part of the name of the file that has to do with that (for technical reasons) and take the file from the correct folder		
		if(("covRemoved" %in% name_spl))
		{
			adrs <- "~/Analyze4C/rawData/coverage_removed/"
		}
		else
		{
			adrs <- "~/Analyze4C/rawData/rearranged/"
		}
		rear_ind <- grep("rearranged",name_spl)
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

		#finding which is the row in the 'rearranged_rawData' that contains information on the particular p-score data:
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
		
		#creating bed files
		write.table(raw_dat,"~/Analyze4C/temp/rawBed.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)		
	}
			
	#for(j in 1:ln)
	for(j in unique(dat[,1]))
	{
		#getting the 
		data.chr <- dat[dat[,1]==j,]
				
		#removing the translocated sections
		if(flag==1 & file.name!="0")
		{				
			#removing the sections from the data
			#for(r in seq(1,length(lin_nums2),by=2))
			#{
			#	first0 <- lin_nums2[r]
			#	second0 <- lin_nums2[r+1]
			#	sq <- c(first0:second0,rep(0,nrow(data.chr)-length(first0:second0))) #creates a sequence with the section line numbers and completes it to be the size of data.chr for tecnincal reasons
			#	data.chr <- data.chr[!rownames(data.chr)==sq,]			
			#}
			#rownames(data.chr) <- seq(length=nrow(data.chr))

			#removing the sections from the data:
			
			data.chr <- cbind(data.chr[,1],data.chr[,2],data.chr[,2]+1,data.chr[,3])
			write.table(data.chr,"~/Analyze4C/temp/psBed.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	
			#intersecting the bed files so that we get the indices that aren't removed (since the way we removed was with line numbers of the files, we can't just remove them from the p-score data since we removed some lines from the cis when calculating the p-scores)
			system("bedtools intersect -a ~/Analyze4C/temp/psBed.bed -b ~/Analyze4C/temp/rawBed.bed > ~/Analyze4C/temp/out.bed")
			system("bedtools intersect -a ~/Analyze4C/temp/psBed.bed -b ~/Analyze4C/temp/rawBed.bed -c > ~/Analyze4C/temp/out2.bed") #this output will have all the same lines as psBed.bed but with an indication if the indices intersected or not (0 not intersected, 1 intersected)
			#there was a problem with intersected repeats, where if the same indices appeared in the translocated area that needs to be removed and an area that doesn'table
			#this will cause repeats in the intersections and some indices to be out of order
			#i couldn't figure out an automatic way to fix this
			#for now the user will fix it manually
			if(flag_rmv != 2)
			{
				cat("\nchromosome ",j,":\n",sep="")
				cat("\nview the 'out.bed' file for any repeats of intesections.\nthis happens if there is an index that appears in the translocated area and in areas that aren't.\n")
				cat("bedtools doesn't know how to distinguish between 'good' intersections like these and 'bad' ones.\nthese repeats should be at the beginning of the file (if the translocation occured on the left telomere of the chromosome.\n")
				cat("you could confirm that these are the lines that need to be erased by looking at the 'out2.bed' which shows that there are some indices that appear early on and in translocated areas that have been intersected.\n")
				cat("most chances are that you will see a '1' at the line of some indices with other indices that have a '0' surrounding them\n")
				#cat("erase these lines\n")
				#ent <- readline(prompt=cat("press ENTER once this is done in order to continue\n\n"))
				lns <- c()
				fl <- 1
				cat("\n\nentering the lines that need to be removed:\n")
				while(fl == 1)
				{
					ans <- readline(prompt=cat("\nare there any (more) lines that need to be erased?\ny/n\n\n"))
					if(ans == "y")
					{
						new_ln <- as.integer(readline(prompt=cat("\nenter the line number that needs to be removed, as it appears in 'out.bed':\n")))
						lns <- c(lns,new_ln)
					}
					else
					{
						fl <- 0
					}
				}
			}
			else
			{
				lns <- in_lns[[j]][[2]]
			}
			ind_counter <- ind_counter + 1
			data.chr <- read.table("~/Analyze4C/temp/out.bed")
			#if(is.null(lns))
			#{
			#	lns <- NA
			#}
			#else
			if(!is.null(lns))
			{
				data.chr <- data.chr[-1*lns,] #removing the specific lines
			}	
			data.chr <- cbind(data.chr[,1],data.chr[,2],data.chr[,4])

			lst_temp <- list(j,lns)
			rmv_lst[[ind_counter]] <- lst_temp

			#removing the bed files created
			system("rm ~/Analyze4C/temp/psBed.bed")
			system("rm ~/Analyze4C/temp/out.bed")
			system("rm ~/Analyze4C/temp/out2.bed")			
		}

		#getting the quantile p-score
		topPer <- unname(quantile(data.chr[,3],co_per))
	
		#gets the size of the data
		len <- nrow(data.chr)
		
		#first gets the beginning of each band
		first <- 0
		
		#has the previous site with positive RE site
		prev <- 0
		
		#counts number of RE sites with 0 in a row
		zero_REs <- 0

		for(i in 1:len)
		{
			#if this is the first positive RE site dealt with in the genome
			if(data.chr[i,3] > topPer && first == 0)
			{
				first <- data.chr[i,]
				prev <- first
				last <- first
				zero_REs <- 0
			}
			else if(all(first != 0)) #the 'all' makes it easier to test if 'first' is true
			{
					#if we reached a positive RE (above cutoff) and we are less than RE_gap from the last one and less than bp_gap bp from the last one
					if((data.chr[i,3] > topPer) && (zero_REs <= RE_gap) && ((data.chr[i,2]-prev[2]) <= bp_gap))
					{
						last <- data.chr[i,]
						prev <- last
					}
					#if we reached a positive RE (above cutoff) and we are more than RE_gap from the last one or more than bp_gap bp from the last one
					else if((data.chr[i,3] > topPer) && ((zero_REs > RE_gap) || ((data.chr[i,2]-prev[2]) > bp_gap)))
					{
						if((first[2] == last[2]) && ((first[2] + 1) <= genome_sizes[j,2]))
						{
							last[2] <- first[2] + 1
						}
						band <- c(first[1],first[2],last[2])
						conts <- rbind(conts,band)
						first <- data.chr[i,]
						prev <- first
						last <- first
						zero_REs <- 0
					}
					#if we reached a band that isn't above cutoff
					else if(data.chr[i,3] <= topPer)
					{
						zero_REs <-  zero_REs + 1
					}
			}
		}

		band <- c(first[1],first[2],last[2])
		conts <- rbind(conts,band)
	}
	
	if(flag==1 & file.name!="0")
	{	
		system("rm ~/Analyze4C/temp/rawBed.bed")
	}	
	return(list(conts,rmv_lst))
}
