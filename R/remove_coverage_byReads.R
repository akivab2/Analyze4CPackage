#' @export

#this function removes coverage by changing the minimum number of reads that will be considered as positive reads, according the the users request

remove_coverage_byReads <- function(raw_dat,Experiments_4C,file.name,cis)
{
	flag1 <- 0
	while(flag1 == 0)
	{
		#get the coverage and reads stats for the original data
		ans1 <- readline(prompt=cat("\nwould you like to get coverage and reads stats for the original data?\ny/n\n\n"))
		if(ans1=="y")
		{
			#coverage vs. minimum reads line plot
			cat("\nplot of minimum reads against coverages:\n\n")
			coverageVSmin_linePlot(raw_dat,file.name,cis)

			#getting the distribution of number of reads
			cat("\ndistribution plots: (note that the data at first is with the minimum reads input from before. if the user will choose to change file then the number of minimum reads will be erased as long as you are in the main menu of the distribution plots)\n\n")
			distribution_plotter(Experiments_4C,raw_dat,file.name,1)
			
			flag2 <- 0
			while(flag2==0)
			{
				#getting the minimum number of reads that will consider the RE site as positive
				ans2 <- as.integer(readline(prompt=cat("\nenter a minimum number of reads that will consider an RE site as positive (1 read is considered the smallest option):\n\n")))
				#getting the coverage stats
				raw_dat_temp <- raw_dat
				raw_dat_temp[raw_dat_temp[,3]<ans2,3] <- 0
				rawData.stats(Experiments_4C,raw_dat_temp,file.name)
												
				ans3 <- readline(prompt=cat("\nwould you like to choose another minimum number of reads to get the stats for?\ny/n\n\n"))
				if(ans3=="n")
				{
					flag2 <- 1
				}
			}	
		}
		
		new_dat <- raw_dat
		#choose what sections of the data we want to perform this on
		choice <- as.integer(readline(prompt=cat("\nchoose what sections you would like to perform the coverage removal on:\n1) all\n2) trans\n3) by chromosome\n\n")))
		if(choice==2) #trans
		{			
			#get the number of minimum reads, after looking at stats
			min_reads <- as.integer(readline(prompt=cat("\nenter the minimum number of reads that will be considered positive:\n")))
		
			new_dat[new_dat[,1]!=cis & new_dat[,3]<min_reads,3] <- 0

			sec <- "trans chromosomes coverage changed"
		}
		else if(choice==3) #by chromosome
		{	
			#choosing which chromosome to remove coverage from. we assume that the chromosome number is legal
			choice2 <- as.integer(readline(prompt=print("\nenter the chromosome number you would like to perform the coverage removal on (cis is chromosome",cis,"):\n\n")))

			#get the number of minimum reads, after looking at stats
			min_reads <- as.integer(readline(prompt=cat("\nenter the minimum number of reads that will be considered positive:\n")))
		
			new_dat[new_dat[,1]!=choice2 & new_dat[,3]<min_reads,3] <- 0
			
			if(choice2==cis)
			{
				sec <- "cis chromosome coverage changed"
			}
			else
			{
				sec <- paste("chromosome",ch,"coverage removed")
			}
		}
		else
		{
			#get the number of minimum reads, after looking at stats
			min_reads <- as.integer(readline(prompt=cat("\nenter the minimum number of reads that will be considered positive:\n")))
		
			new_dat[new_dat[,3]<min_reads,3] <- 0
				
			sec <- "all chromosomes coverage changed"
		}
			
		#get the coverage and reads stats for the new data
		ans4 <- readline(prompt=cat("\nwould you like to get coverage and reads stats for the data after removal of reads?\ny/n\n\n"))
		if(ans4=="y")
		{
				#getting the coverage stats
				rawData.stats(Experiments_4C,new_dat,file.name)	
		}
		
		ans7 <- readline(prompt=cat("\nare you satisfied with this coverage removal?\ny/n\n\n"))
		if(ans7=="y")
		{
			flag1 <- 1
		}
	}
	out <- list(new_dat,sec,min_reads)
	return(out)
}
