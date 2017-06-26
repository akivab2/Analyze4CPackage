#' @export

#this function calculates the value at a specific quantile the user chooses,
#or calculates the quantile of a specific value
#the user choose for which type of data file this should be calculated, and for which chromosomes

calculate_quantiles <- function(Experiments_4C)
{
	cat("\ncalculating quantiles:\n\n")
	#ask if to do it for raw data, p-score data, rnaseq data, or chipseq data
	ans1 <- as.integer(readline(prompt=cat("\nchoose the type of file you would like to calculate the quantiles for:\n1) raw data\n2) p-scores\n3) RNA-seq data\n4) ChIP-seq data\n\n")))
	#get the data
	if(ans1 == 1) #raw data 
	{
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
			dat <- read.table(paste("~/Analyze4C/rawData/rearranged/",file.name,sep=""))
		}
		else if(choice == 3)
		{
			dat <- read.table(paste("~/Analyze4C/rawData/coverage_removed/",file.name,sep=""))
		}
		else
		{
			dat <- read.table(paste("~/Analyze4C/rawData/original/",file.name,sep=""))
		}
		
		#finding the specific raw data files information in Experiments_4C
		sp1 <- unlist(strsplit(file.name,"[.]"))
		sp2 <- strsplit(sp1[[1]][1],"_")
		out <- findIn.Experiments_4C(sp2[[1]][1],sp2[[1]][2],sp2[[1]][3],Experiments_4C)
		
		#getting the cis chromosome number
		cis <- as.numeric(out[2])
		
		#cis_available gets 1 since there is a cis
		cis_available <- 1
		
		#getting the index of the column with the values of the data
		val_col <- 3	
	}
	else if(ans1 == 2) #p-scores
	{
		repeat
		{
			#getting the name of the p-score file the user wants to use
			choice <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the p-score file from:\n\n1) calculated using bps\n2) calculated without using bps\n\n")))
			if(choice == 1)
			{
				ls_files <- system("ls ~/Analyze4C/pScores/with_bp",intern=TRUE)
			}
			else
			{
				ls_files <- system("ls ~/Analyze4C/pScores/no_bp",intern=TRUE)
			}
						
			if(length(ls_files) != 0)
			{
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
			else
			{
				cat("\nthis folder is empty. please choose a different folder from which you would like to use a file\n\n")
			}	
		}
		
		#importing the data from the file
		if(choice == 1)
		{
			dat <- read.table(paste("~/Analyze4C/pScores/with_bp/",file.name,sep=""))
		}
		else
		{
			dat <- read.table(paste("~/Analyze4C/pScores/no_bp/",file.name,sep=""))
		}
		
		#finding the specific p-score files information in Experiments_4C
		sp <- unlist(strsplit(file.name,"_"))
		out <- findIn.Experiments_4C(sp[1],sp[2],sp[3],Experiments_4C)
		
		#getting the cis chromosome number
		cis <- as.numeric(out[2])
		
		#cis_available gets 1 since there is a cis
		cis_available <- 1	
		
		#getting the index of the column with the values of the data
		val_col <- 3		
	}
	else if(ans1 == 3) #rna-seq
	{
		repeat
		{
			ls_files <- system("ls ~/Analyze4C/RNAseq/FPKM | grep 'FPKM.bed$'",intern=TRUE)
			if(length(ls_files) != 0)
			{	
				file.name <- readline(prompt=cat("\nThese are the FPKM files available\nwhich would you like to use?\n",ls_files,"\n",sep="\n"))
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
				cat("\nthis folder is empty. there are no RNA-seq files available\n\n")
				return()
			}
		}	
		#importing the data from the file
		dat <- read.table(paste("~/Analyze4C/RNAseq/FPKM/",file.name,sep=""))
				
		#cis_available gets 0 since there is no cis with this type of data
		cis_available <- 0	

		#getting the index of the column with the values of the data
		val_col <- 4		
	}
	else if(ans1 == 4) #chip-seq
	{
		repeat
		{
			choice <- as.integer(readline(prompt=cat("\nWhich type of ChIP-seq file you would to choose (what values should it contain)?:\n\n1) tags\n2) p-scores\n\n")))
			if(choice == 1)
			{
				ls_files <- system("ls ~/Analyze4C/ChIPseq/peaks | grep '_tags.bed$'",intern=TRUE)
			}
			else if(choice == 2)
			{
				ls_files <- system("ls ~/Analyze4C/ChIPseq/peaks | grep '_pscores.bed$'",intern=TRUE)
			}
			
			if(length(ls_files) != 0)
			{	
				file.name <- readline(prompt=cat("\nThese are the ChIP-seq (peaks) files available\nwhich would you like to use?\n",ls_files,"\n",sep="\n"))
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
				cat("\nthere are no ChIP-seq files of this sort available\n\n")
				return()
			}
		}
		
		#importing the data from the file
		dat <- read.table(paste("~/Analyze4C/ChIPseq/peaks/",file.name,sep=""))			
				
		#cis_available gets 0 since there is no cis with this type of data
		cis_available <- 0

		#getting the index of the column with the values of the data
		val_col <- 4		
	}	

	#ask if to calculate for each chromosome individually, all, trans, or cis
	if(cis_available == 1)
	{
		ans4 <- as.integer(readline(prompt=cat("\nchoose which chromosomes you would like to calculate for?\n1) all chromosomes together\n2) all trans together\n3) only cis\n4) each chromosome separately\n\n")))
		if(ans4 == 1)
		{
			dat_temp <- dat
			chosen_chroms <- "all chromosomes together"	
		}
		else if(ans4 == 2)
		{
			dat_temp <- dat[dat[,1]!=cis,]
			chosen_chroms <- "trans chromosomes together"
		}
		else if(ans4 == 3)
		{
			dat_temp <- dat[dat[,1]==cis,]
			chosen_chroms <- "just the cis chromosome"
		}	
		else if(ans4 == 4)
		{
			dat_temp <- dat
			chosen_chroms <- "each chromosome separately"
		}
	}
	else
	{
		ans4 <- as.integer(readline(prompt=cat("\nchoose which chromosomes you would like to calculate for?\n1) all chromosomes together\n2) each chromosome separately\n\n")))
		if(ans4 == 1)
		{
			dat_temp <- dat
			chosen_chroms <- "all chromosomes together"	
		}	
		else if(ans4 == 2)
		{
			dat_temp <- dat
			chosen_chroms <- "each chromosome separately"
			ans4 <- 4 #this is for later use, so it will be similar to the options above
		}	
	}
	
	ans <- 0
	while(ans != 4)
	{
		#main menu
		cat("\nthe chosen file is",file.name,"\n")
		cat("\nthe chosen chromosomes to look at are",chosen_chroms,"\n\n")
		cat("\ncalculating quantiles menu:\n\n")
		cat("1) choose method of testing\n")
		cat("2) choose for which chromosomes to check\n")
		cat("3) choose a different file\n")		
		cat("4) exit\n\n")
		ans <- as.integer(readline())
		
		if(ans == 1)
		{
			#ask if to calculate the quantiles by giving a quantile number or using a cutoff get the quantiles
			ans3 <- as.integer(readline(prompt=cat("\nplease choose an option:\n1) calculate the value at a specific quantile\n2) get the quantile of a specific value\n\n")))
			if(ans3 == 1)
			{
				quant <- as.numeric(readline(prompt=cat("\nenter the quantile you would like to check (between 0 and 1):\n\n")))
				#calculate and print quantile
				if(ans4 == 4)
				{
					for(u in unique(dat[,1]))
					{
						cat("\nchromosome ",u,":\n",sep="")
						cat(quantile(dat[dat[,1]==u,val_col],quant))
						cat("\n")
					}
				}
				else
				{
					cat("\n")
					cat(quantile(dat[,val_col],quant))
					cat("\n")
				}								
			}
			else
			{
				val <- as.numeric(readline(prompt=cat("\nenter the value at which you would like to check the quantile for:\n\n")))
				#calculate and print quantile
				if(ans4 == 4)
				{
					for(u in unique(dat[,1]))
					{
						cat("\nchromosome ",u,":\n",sep="")
						cat(ecdf(dat[dat[,1]==u,val_col])(val))
						cat("\n")
					}
				}
				else
				{
					cat("\n")
					cat(ecdf(dat[,val_col])(val))
					cat("\n")
				}	
			}			
		}
		else if(ans == 2)
		{
			#ask if to calculate for each chromosome individually, all, trans, or cis
			if(cis_available == 1)
			{
				ans4 <- as.integer(readline(prompt=cat("\nchoose which chromosomes you would like to calculate for?\n1) all chromosomes together\n2) all trans together\n3) only cis\n4) each chromosome separately\n\n")))
				if(ans4 == 1)
				{
					dat_temp <- dat
					chosen_chroms <- "all chromosomes together"	
				}
				else if(ans4 == 2)
				{
					dat_temp <- dat[dat[,1]!=cis,]
					chosen_chroms <- "trans chromosomes together"
				}
				else if(ans4 == 3)
				{
					dat_temp <- dat[dat[,1]==cis,]
					chosen_chroms <- "just the cis chromosome"
				}	
				else if(ans4 == 4)
				{
					dat_temp <- dat
					chosen_chroms <- "each chromosome separately"
				}
			}
			else
			{
				ans4 <- as.integer(readline(prompt=cat("\nchoose which chromosomes you would like to calculate for?\n1) all chromosomes together\n2) each chromosome separately\n\n")))
				if(ans4 == 1)
				{
					dat_temp <- dat
					chosen_chroms <- "all chromosomes together"	
				}	
				else if(ans4 == 2)
				{
					dat_temp <- dat
					chosen_chroms <- "each chromosome separately"
					ans4 <- 4 #this is for later use, so it will be similar to the options above
				}	
			}	
		}
		else if(ans == 3)
		{
			#ask if to do it for raw data, p-score data, rnaseq data, or chipseq data
			ans1 <- as.integer(readline(prompt=cat("\nchoose the type of file you would like to calculate the quantiles for:\n1) raw data\n2) p-scores\n3) RNA-seq data\n4) ChIP-seq data\n\n")))
			#get the data
			if(ans1 == 1) #raw data 
			{
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
					dat <- read.table(paste("~/Analyze4C/rawData/rearranged/",file.name,sep=""))
				}
				else if(choice == 3)
				{
					dat <- read.table(paste("~/Analyze4C/rawData/coverage_removed/",file.name,sep=""))
				}
				else
				{
					dat <- read.table(paste("~/Analyze4C/rawData/original/",file.name,sep=""))
				}
				
				#finding the specific raw data files information in Experiments_4C
				sp1 <- unlist(strsplit(file.name,"[.]"))
				sp2 <- strsplit(sp1[[1]][1],"_")
				out <- findIn.Experiments_4C(sp2[[1]][1],sp2[[1]][2],sp2[[1]][3],Experiments_4C)
				
				#getting the cis chromosome number
				cis <- as.numeric(out[2])
				
				#cis_available gets 1 since there is a cis
				cis_available <- 1	

				#getting the index of the column with the values of the data
				val_col <- 3				
			}
			else if(ans1 == 2) #p-scores
			{
				repeat
				{
					#getting the name of the p-score file the user wants to use
					choice <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the p-score file from:\n\n1) calculated using bps\n2) calculated without using bps\n\n")))
					if(choice == 1)
					{
						ls_files <- system("ls ~/Analyze4C/pScores/with_bp",intern=TRUE)
					}
					else
					{
						ls_files <- system("ls ~/Analyze4C/pScores/no_bp",intern=TRUE)
					}
								
					if(length(ls_files) != 0)
					{
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
					else
					{
						cat("\nthis folder is empty. please choose a different folder from which you would like to use a file\n\n")
					}	
				}
				
				#importing the data from the file
				if(choice == 1)
				{
					dat <- read.table(paste("~/Analyze4C/pScores/with_bp/",file.name,sep=""))
				}
				else
				{
					dat <- read.table(paste("~/Analyze4C/pScores/no_bp/",file.name,sep=""))
				}
				
				#finding the specific p-score files information in Experiments_4C
				sp <- unlist(strsplit(file.name,"_"))
				out <- findIn.Experiments_4C(sp[1],sp[2],sp[3],Experiments_4C)
				
				#getting the cis chromosome number
				cis <- as.numeric(out[2])
				
				#cis_available gets 1 since there is a cis
				cis_available <- 1	

				#getting the index of the column with the values of the data
				val_col <- 3				
			}
			else if(ans1 == 3) #rna-seq
			{
				repeat
				{
					ls_files <- system("ls ~/Analyze4C/RNAseq/FPKM | grep 'FPKM.bed$'",intern=TRUE)
					if(length(ls_files) != 0)
					{	
						file.name <- readline(prompt=cat("\nThese are the FPKM files available\nwhich would you like to use?\n",ls_files,"\n",sep="\n"))
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
						cat("\nthis folder is empty. there are no RNA-seq files available\n\n")
						break
					}
				}
				
				if(length(ls_files) != 0)
				{
					#importing the data from the file
					dat <- read.table(paste("~/Analyze4C/RNAseq/FPKM/",file.name,sep=""))
							
					#cis_available gets 0 since there is no cis with this type of data
					cis_available <- 0
					
					#getting the index of the column with the values of the data
					val_col <- 4					
				}	
			}
			else if(ans1 == 4) #chip-seq
			{
				repeat
				{
					choice <- as.integer(readline(prompt=cat("\nWhich type of ChIP-seq file you would to choose (what values should it contain)?:\n\n1) tags\n2) p-scores\n\n")))
					if(choice == 1)
					{
						ls_files <- system("ls ~/Analyze4C/ChIPseq/peaks | grep '_tags.bed$'",intern=TRUE)
					}
					else if(choice == 2)
					{
						ls_files <- system("ls ~/Analyze4C/ChIPseq/peaks | grep '_pscores.bed$'",intern=TRUE)
					}
					
					if(length(ls_files) != 0)
					{	
						file.name <- readline(prompt=cat("\nThese are the ChIP-seq (peaks) files available\nwhich would you like to use?\n",ls_files,"\n",sep="\n"))
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
						cat("\nthere are no ChIP-seq files of this sort available\n\n")
					}
				}
				
				if(length(ls_files) != 0)
				{				
					#importing the data from the file
					dat <- read.table(paste("~/Analyze4C/ChIPseq/peaks/",file.name,sep=""))			
							
					#cis_available gets 0 since there is no cis with this type of data
					cis_available <- 0
					
					#getting the index of the column with the values of the data
					val_col <- 4					
				}		
			}			
		}	
	}	
}
