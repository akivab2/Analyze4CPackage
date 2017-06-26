#' @export

#this function calculates the different stats of the raw data files.
#including the coverages of cis, trans, all, any specific chromosome, and also could determine if the self has too many reads compared to the rest of the data

rawData.stats <- function(Experiments_4C,rawData=-1,file.name=NA) 
{
	cat("\n\n********* Raw data stats *********\n\n\n");

	if(all(rawData == -1))
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
	}
	
	#finding the specific raw data files information in Experiments_4C
	sp1 <- unlist(strsplit(file.name,"[.]"))
	sp2 <- strsplit(sp1[[1]][1],"_")
	out <- findIn.Experiments_4C(sp2[[1]][1],sp2[[1]][2],sp2[[1]][3],Experiments_4C)
	
	#getting the cis chromosome number
	cis.chrom <- as.numeric(out[2])
	
	#getting the bait index
	bait.ind <- as.numeric(out[3])
	
	#cis.chrom <- readline(prompt="Please enter the cis chromosome number:\n");
	#cis.chrom <- as.numeric(cis.chrom);
	
	#bait.ind <- readline(prompt="Please enter the bait position:\n");
	#bait.ind <- as.numeric(bait.ind);

	#repeating loop until the user decides to end the program
	flag <- 1;
	
	while(flag)
	{
		#getting the choice number from user
		choice <- readline(prompt="\nRaw stats menu\nPlease enter number of choice:\n\n1 - full coverage\n2 - trans coverage\n3 - cis coverage\n4 - coverage by chromosome\n5 - self\n6 - end\n\n");
		choice <- as.numeric(choice);

		switch(choice,
		#prints out the whole coverage of the rawData
		{
			all.cov <- sum(rawData[,3]>0)/sum(rawData[,3]>=0);
			cat("Whole coverage:\n",all.cov,"\n");
		},
		#prints out the trans coverage
		{
			trans.cov <- sum(rawData[rawData[,1]!=cis.chrom,3]>0)/sum(rawData[rawData[,1]!=cis.chrom,3]>=0)
			cat("trans coverage:\n",trans.cov,"\n");
		},
		#prints out the cis coverage
		{
			cis.cov <- sum(rawData[rawData[,1]==cis.chrom,3]>0)/sum(rawData[rawData[,1]==cis.chrom,3]>=0)
			cat("cis coverage:\n",cis.cov,"\n");			
		},
		#prints out the coverage of chromosome by choice
		{
			ch.num <- readline(prompt="please enter chromosome number that you would like to calculate coverage for:\n")
			ch.num <- as.numeric(ch.num);
			ch.cov <- sum(rawData[rawData[,1]==ch.num,3]>0)/sum(rawData[rawData[,1]==ch.num,3]>=0);
			cat("chromosome",ch.num,"coverage:\n",ch.cov,"\n");
		},
		#getting the number of reads for all chromosomes, for self, the percentage of self out of all, and determines if the self is too much (based on the self.cutoff) 
		{
			#checking if the self RE site exists in the data or was removed
			self.reads <- rawData[rawData[,1]==cis.chrom & rawData[,2]==bait.ind,3]
			if(identical(self.reads,integer(0))) #if the self RE site was removed from the data
			{
				cat("\nerror:\nthe self cannot be examined on this data since it does not exist here\ntry checking it with the original raw data file, before the removal of the RE sites\n\n")
			}
			else #if the self exists in the data
			{
				#getting the self cutoff percentage number from user
				self.cutoff <- as.numeric(readline(prompt="Please enter the self cutoff percentage (between 0 and 100):\n"));
				all.reads <- sum(rawData[,3]);
				cat("The number of reads in all chromosomes:\n",all.reads,"\n");
				cat("The number of reads on the self are:\n",self.reads,"\n");
				self.per <- (self.reads/all.reads)*100;
				cat("The percentage of self reads out of all reads are:\n",self.per,"%\n",sep="");
				if(self.per>=self.cutoff)
				{
					cat("The self has more reads than it should. This is a sign that there were too many self ligations and the restriction enzyme was not efficient enough\n");
				}
				else
				{cat("The number of self reads seems reasonable\n");}
			}
		},	
		#end the program
		{flag <- 0;}
		)
	}
	cat("Good Luck :)\n\n");	
}
