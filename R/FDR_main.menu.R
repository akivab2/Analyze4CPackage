#' @export

#need to add an option of testing the data only on cis or trans, or on specific chromosomes
FDR_main.menu <- function(Experiments_4C)
{
	#create menu: before menu you need to choose a file to use, then the menu asks what action you want to do, one option at the end asks if you want to change files to use
	
	cat("\n*******************\n   FDR Main Menu   \n*******************\n\n")
	
	#maybe create a function for all this
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
	}
	
	#checking if the raw data taken has its indices in sequence (they might not be when rearranged, for translocations for example)
	#this is in order to see if the user can calculate p-scores using windows with bp size
	inp <- 0 #gets a default of 0
	num_of_chr <- length(unique(rawData[,1])); #getting the number of chromosomes
	for(k in 1:num_of_chr)
	{
		if(is.unsorted(rawData[rawData[,1]==k,2]))
		{
			cat("\nthe raw data seems to be rearranged.\nyou may only do actions that when calculating p-scores it does so only with windows by RE sites (and not bp).\n\n")
			inp <- 1
			break
		}
	}

	
	flag0 <- 1
	while(flag0)
	{	
		choice1 <- as.integer(readline(prompt=cat("\nFDR Menu:\n1) get the average FDR for each chromosome using cutoffs using windows of bp and RE sites\n2) get the average FDR for each chromosome using cutoffs using only windows of RE sites\n3) get the lowest cutoff percentiles (in average and median) that is closest to a constant FDR value given by user, for windows of bp and RE sites\n4) get the lowest cutoff percentiles (in average and median) that is closest to a constant FDR value given by user, only for windows RE sites\n5) calculate the window size that will give the lowest FDR for a specific cutoff percent, for windows of bp and RE sites\n6) calculate the window size that will give the lowest FDR for a specific cutoff percent, only for windows of RE sites\n7) calculate FDRs with a range of cutoffs and range of window sizes, for windows of bp and RE sites\n8) calculate FDRs with a range of cutoffs and range of window sizes, only for windows of RE sites\n9) choose a different raw data file\n10) exit menu\n\n")))
		switch(choice1,
			{#1) FDR
				if(inp == 1) #if the raw data files indices are not in order
				{
					cat("\nyou cannot do this action\nthe raw data chosen is probably rearranged with indices out of order\nany action on this raw data that includes creating sliding windows, must do so only by RE sites and bp\nin order to get the average FDRs press 2\n\n")
				}
				else
				{
					#get the p score file that the user wants to use
					#with kb
					ps_files <- system("ls ~/Analyze4C/pScores/with_bp",intern=TRUE)
					if((identical(ps_files,character(0)) == "FALSE"))
					{	
						{
							#getting the name of the p score file the user wants to use
							repeat
							{
								ps.file <- readline(prompt=cat("These are the p-score files available\nwhich would you like to use? (make sure to choose one that is of the same experiment as the raw data)\n",ps_files,"\n",sep="\n"))
								ind_files2 <- pmatch(ps.file,ps_files,nomatch=0)
								if(ind_files2 == 0)
								{
									cat("no such file exists.\nplease try again.\n\n")
								}
								else
								{
									break
								}
							}

							#importing the p-score data from the file 	
							ps <- read.table(paste("~/Analyze4C/pScores/with_bp/",ps.file,sep=""))

							#finding the specific raw data files information in Experiments_4C
							sp3 <- unlist(strsplit(ps.file,"[.]"))
							sp4 <- strsplit(sp3[[1]][1],"_")
							
							#getting the size of the window per RE sites
							sp5 <- sp4[[1]][which(grepl("^[[:digit:]]+RE$",sp4[[1]]))]
							RE_wind <- strsplit(sp5[[1]][1],"RE")
							RE_wind <- as.numeric(RE_wind)
							
							#getting the size of the window per bp sites
							sp6 <- sp4[[1]][which(grepl("^[[:digit:]]+bp$",sp4[[1]]))]
							bp_wind <- strsplit(sp6[[1]][1],"bp")
							bp_wind <- as.numeric(bp_wind)
						}
										
						#get the number of cycles from user
						cyc <- readline(prompt=cat("\nplease enter the number of cycles you would like to have:\n"))
						cyc <- as.numeric(cyc)
						#get the percentile cutoff from user
						per <- readline(prompt=cat("\nplease enter cutoff percentile between 0 and 1:\n"))
						per <- as.numeric(per)
						#getting the FDR
						FDR(ps,rawData,vp.chrom,vp.pos,erase,bp_wind,RE_wind/2,cyc,per)
					}
					else
					{
						cat("\nthere are no available files in the 'with_bp' folder of the p-scores\ntry a different option\n\n")
					}	
				}	
			},
			{#2) FDR_nokb
				#get the p score file that the user wants to use
				#no kb
				ps_files <- system("ls ~/Analyze4C/pScores/no_bp",intern=TRUE)
				if((identical(ps_files,character(0)) == "FALSE"))
				{
					{
						#getting the name of the p score file the user wants to use
						repeat
						{
							ps.file <- readline(prompt=cat("These are the p-score files available\nwhich would you like to use? (make sure to choose one that is of the same experiment as the raw data)\n",ps_files,"\n",sep="\n"))
							ind_files2 <- pmatch(ps.file,ps_files,nomatch=0)
							if(ind_files2 == 0)
							{
								cat("no such file exists.\nplease try again.\n\n")
							}
							else
							{
								break
							}
						}

						#importing the p-score data from the file 	
						ps <- read.table(paste("~/Analyze4C/pScores/no_bp/",ps.file,sep=""))

						#finding the specific raw data files information in Experiments_4C
						sp3 <- unlist(strsplit(ps.file,"[.]"))
						sp4 <- strsplit(sp3[[1]][1],"_")
						
						#getting the size of the window per RE sites
						sp5 <- sp4[[1]][which(grepl("^[[:digit:]]+RE$",sp4[[1]]))]
						RE_wind <- strsplit(sp5[[1]][1],"RE")
						RE_wind <- as.numeric(RE_wind)
					}

					#get the number of cycles from user
					cyc <- readline(prompt=cat("\nplease enter the number of cycles you would like to have:\n"))
					cyc <- as.numeric(cyc)
					#get the percentile cutoff from user
					per <- readline(prompt=cat("\nplease enter cutoff percentile between 0 and 1:\n"))
					per <- as.numeric(per)
					#calculating the FDR
					FDR_nokb(ps,rawData,vp.chrom,vp.pos,erase,RE_wind/2,cyc,per)
				}
				else
				{
					cat("\nthere are no available files in the 'no_bp' folder of the p-scores\ntry a different option\n\n")
				}
			},
			{#3) FDR_COper_finder
				if(inp == 1) #if the raw data files indices are not in order
				{
					cat("\nyou cannot do this action\nthe raw data chosen is probably rearranged with indices out of order\nany action on this raw data that includes creating sliding windows, must do so only by RE sites and bp\nin order to get the average FDRs press 2\n\n")
				}
				else
				{
					#get the p score file that the user wants to use
					#with kb
					ps_files <- system("ls ~/Analyze4C/pScores/with_bp",intern=TRUE)
					if((identical(ps_files,character(0)) == "FALSE"))
					{						
						{
							#getting the name of the p score file the user wants to use
							repeat
							{
								ps.file <- readline(prompt=cat("These are the p-score files available\nwhich would you like to use? (make sure to choose one that is of the same experiment as the raw data)\n",ps_files,"\n",sep="\n"))
								ind_files2 <- pmatch(ps.file,ps_files,nomatch=0)
								if(ind_files2 == 0)
								{
									cat("no such file exists.\nplease try again.\n\n")
								}
								else
								{
									break
								}
							}

							#importing the p-score data from the file 	
							ps <- read.table(paste("~/Analyze4C/pScores/with_bp/",ps.file,sep=""))

							#finding the specific raw data files information in Experiments_4C
							sp3 <- unlist(strsplit(ps.file,"[.]"))
							sp4 <- strsplit(sp3[[1]][1],"_")
							
							#getting the size of the window per RE sites
							sp5 <- sp4[[1]][which(grepl("^[[:digit:]]+RE$",sp4[[1]]))]
							RE_wind <- strsplit(sp5[[1]][1],"RE")
							RE_wind <- as.numeric(RE_wind)
							
							#getting the size of the window per bp sites
							sp6 <- sp4[[1]][which(grepl("^[[:digit:]]+bp$",sp4[[1]]))]
							bp_wind <- strsplit(sp6[[1]][1],"bp")
							bp_wind <- as.numeric(bp_wind)
						}
										
						#get the number of cycles from user
						cyc <- readline(prompt=cat("\nplease enter the number of cycles you would like to have:\n"))
						cyc <- as.numeric(cyc)
						#get the percentile cutoff from user
						per <- readline(prompt=cat("\nplease enter cutoff percentile between 0 and 1:\n"))
						per <- as.numeric(per)
						#get the FDR from user
						fdr <- readline(prompt=cat("\nplease enter the FDR between 0 and 1:\n"))
						fdr <- as.numeric(fdr)
						#getting the percentage of removal from the initial cutoff that will be used everytime the cutoff is lowered to get closer to the FDR
						rem <- as.numeric(readline(prompt=cat("\nplease enter a percentage removal (between 0 and 1) that will be used to reduce the cutoff until the FDR is surpassed:\n")))
						#getting the best cutoffs for each chromosome
						cat("\ncalculating the cutoffs and FDRs\nplease wait...\n\n")	
						FDR_COper_finder(ps,rawData,vp.chrom,vp.pos,erase,bp_wind,RE_wind/2,cyc,per,fdr,rem)
					}
					else
					{
						cat("\nthere are no available files in the 'with_bp' folder of the p-scores\ntry a different option\n\n")
					}	
				}	
			},
			{#4) FDR_COper_finder_nokb
				#get the p score file that the user wants to use
				#no kb
				ps_files <- system("ls ~/Analyze4C/pScores/no_bp",intern=TRUE)
				if((identical(ps_files,character(0)) == "FALSE"))
				{				
					{
						#getting the name of the p score file the user wants to use
						repeat
						{
							ps_files <- system("ls ~/Analyze4C/pScores/no_bp",intern=TRUE)
							ps.file <- readline(prompt=cat("These are the p-score files available\nwhich would you like to use? (make sure to choose one that is of the same experiment as the raw data)\n",ps_files,"\n",sep="\n"))
							ind_files2 <- pmatch(ps.file,ps_files,nomatch=0)
							if(ind_files2 == 0)
							{
								cat("no such file exists.\nplease try again.\n\n")
							}
							else
							{
								break
							}
						}

						#importing the p-score data from the file 	
						ps <- read.table(paste("~/Analyze4C/pScores/no_bp/",ps.file,sep=""))

						#finding the specific raw data files information in Experiments_4C
						sp3 <- unlist(strsplit(ps.file,"[.]"))
						sp4 <- strsplit(sp3[[1]][1],"_")
						
						#getting the size of the window per RE sites
						sp5 <- sp4[[1]][which(grepl("^[[:digit:]]+RE$",sp4[[1]]))]
						RE_wind <- strsplit(sp5[[1]][1],"RE")
						RE_wind <- as.numeric(RE_wind)		
					}

					#get the number of cycles from user
					cyc <- readline(prompt=cat("\nplease enter the number of cycles you would like to have:\n"))
					cyc <- as.numeric(cyc)
					#get the percentile cutoff from user
					per <- readline(prompt=cat("\nplease enter cutoff percentile between 0 and 1:\n"))
					per <- as.numeric(per)
					#get the FDR from user
					fdr <- readline(prompt=cat("\nplease enter the FDR between 0 and 1:\n"))				
					fdr <- as.numeric(fdr)
					#getting the percentage of removal from the initial cutoff that will be used everytime the cutoff is lowered to get closer to the FDR
					rem <- as.numeric(readline(prompt=cat("\nplease enter a percentage removal (between 0 and 1) that will be used to reduce the cutoff until the FDR is surpassed:\n")))				
					#getting the best cutoffs for each chromosome
					cat("\ncalculating the cutoffs and FDRs\nplease wait...\n\n")					
					FDR_COper_finder_nokb(ps,rawData,vp.chrom,vp.pos,erase,RE_wind/2,cyc,per,fdr,rem)
				}
				else
				{
					cat("\nthere are no available files in the 'no_bp' folder of the p-scores\ntry a different option\n\n")
				}	
			},
			{#5) FDR_window_cmp
				
				if(inp == 1) #if the raw data files indices are not in order
				{
					cat("\nyou cannot do this action\nthe raw data chosen is probably rearranged with indices out of order\nany action on this raw data that includes creating sliding windows, must do so only by RE sites and bp\nin order to get the average FDRs press 2\n\n")
				}
				else
				{
					#getting the smallest window size in bp
					small.bp <- readline(prompt=cat("\nenter the smallest window size in bp:\n"))
					small.bp <- as.numeric(small.bp)
					#getting the biggest window size in bp
					big.bp <- readline(prompt=cat("\nenter the largest window size in bp:\n"))	
					big.bp <- as.numeric(big.bp)				
					#getting the step for the window in bp
					step.bp	<- readline(prompt=cat("\nenter the step of the growing window sizes in bp:\n"))
					step.bp <- as.numeric(step.bp)				
					#get the ratio between RE site and bp				
					ratio <- readline(prompt=cat("\nplease enter the ratio between the RE sites and bp (every how many bp is there an RE site?):\n"))
					ratio <- as.numeric(ratio)				
					#get the number of cycles from user
					cyc <- readline(prompt=cat("\nplease enter the number of cycles you would like to have:\n"))
					cyc <- as.numeric(cyc)				
					#get the percentile cutoff from user
					per <- readline(prompt=cat("\nplease enter cutoff percentile between 0 and 1:\n"))
					per <- as.numeric(per)
					#gets the best window size that will give the lowest FDR for a specific cutoff
					cat("\ncalculating the FDRs and window sizes\nplease wait...\n\n")
					FDR_window_cmp(rawData,small.bp,big.bp,step.bp,vp.chrom,vp.pos,erase,ratio,cyc,per)
				}	
			},
			{#6) FDR_window_cmp_nokb
			
				#getting the smallest window size in RE sites
				small.RE <- readline(prompt=cat("\nenter the smallest window size in RE sites:\n"))
				small.RE <- as.numeric(small.RE)
				#getting the biggest window size in RE sites
				big.RE <- readline(prompt=cat("\nenter the largest window size in RE sites:\n"))	
				big.RE <- as.numeric(big.RE)				
				#getting the step for the window in RE sites
				step.RE	<- readline(prompt=cat("\nenter the step of the growing window sizes in RE sites:\n"))
				step.RE <- as.numeric(step.RE)							
				#get the number of cycles from user
				cyc <- readline(prompt=cat("\nplease enter the number of cycles you would like to have:\n"))
				cyc <- as.numeric(cyc)				
				#get the percentile cutoff from user
				per <- readline(prompt=cat("\nplease enter cutoff percentile between 0 and 1:\n"))
				per <- as.numeric(per)
				#gets the best window size that will give the lowest FDR for a specific cutoff
				cat("\ncalculating the FDRs and window sizes\nplease wait...\n\n")
				FDR_window_cmp_nokb(rawData,small.RE,big.RE,step.RE,vp.chrom,vp.pos,erase,cyc,per)
			},			
			{#7) FDR_calculator
				
				if(inp == 1) #if the raw data files indices are not in order
				{
					cat("\nyou cannot do this action\nthe raw data chosen is probably rearranged with indices out of order\nany action on this raw data that includes creating sliding windows, must do so only by RE sites and bp\nin order to get the average FDRs press 2\n\n")
				}
				else
				{
					#getting the smallest window size in bp
					small.bp <- readline(prompt=cat("\nenter the smallest window size in bp:\n"))
					small.bp <- as.numeric(small.bp)
					#getting the biggest window size in bp
					big.bp <- readline(prompt=cat("\nenter the largest window size in bp:\n"))	
					big.bp <- as.numeric(big.bp)
					#getting the step for the window in bp
					step.bp	<- readline(prompt=cat("\nenter the step of the growing window sizes in bp:\n"))
					step.bp <- as.numeric(step.bp)
					#getting the smallest cutoff size
					small.CO <- readline(prompt=cat("\nenter the smallest cutoff:\n"))
					small.CO <- as.numeric(small.CO)
					#getting the biggest cutoff size
					big.CO <- readline(prompt=cat("\nenter the largest cutoff:\n"))	
					big.CO <- as.numeric(big.CO)
					#getting the step for the cutoff
					step.CO	<- readline(prompt=cat("\nenter the step of the growing cutoff:\n"))
					step.CO <- as.numeric(step.CO)
					#get the ratio between RE site and bp				
					ratio <- readline(prompt=cat("\nplease enter the ratio between the RE sites and bp (every how many bp is there an RE site?):\n"))				
					ratio <- as.numeric(ratio)
					#get the number of cycles from user
					cyc <- readline(prompt=cat("\nplease enter the number of cycles you would like to have:\n"))
					cyc <- as.numeric(cyc)
					#calculating all the data for the FDR per window sizes and cutoffs
					cat("\ncalculating the cutoffs and FDRs\nplease wait...\n\n")
					FDR_cal <- FDR_calculator(rawData,small.bp,big.bp,step.bp,small.CO,big.CO,step.CO,vp.chrom,vp.pos,erase,ratio,cyc)
					DandT1 <- toString(Sys.time())
					DandT2 <- gsub(" ","_",DandT1)
					DandT2 <- gsub(":","",DandT2)
					#creating a table file of the results
					write.csv(FDR_cal,paste("~/Analyze4C/FDR/",sp1[[1]][1],"_FDRs_",cyc,"cycles_",DandT2,".csv",sep=""),row.names=FALSE)
					#creating a plot for the trans FDRs					
					cat("\ncreating average FDRs for trans plot\nplease wait...\n\n")
					#print(ggplot2::ggplot(data=FDR_cal,ggplot2::aes(x=window_RE,y=fdr.avg.trans,group=CO)) + ggplot2::geom_line(ggplot2::aes(color=CO)) + ggplot2::geom_point() + ggplot2::ggtitle(paste("window size (by RE sites) VS. trans FDRs, per cutoff percentage -",sp1[[1]][1])) + ggplot2::labs(x="window size per RE sites",y="average FDR of trans chromosomes") + ggplot2::theme(plot.title = ggplot2::element_text(size=11)))
					print(ggplot2::ggplot(data=FDR_cal,ggplot2::aes(x=window_RE,y=fdr.avg.trans,group=CO)) + ggplot2::geom_line(ggplot2::aes(color=factor(CO))) + ggplot2::geom_point(size=0.5) + ggplot2::ggtitle(paste("window size (by RE sites) VS. trans FDRs, per cutoff percentage -",sp1[[1]][1])) + ggplot2::labs(x="window size per RE sites",y="average FDR of trans chromosomes") + ggplot2::theme(plot.title = ggplot2::element_text(size=11)))
					ans <- readline(prompt=cat("\nwould you like to save the plot?\ny/n\n\n"))
					if(ans == "y")
					{
						wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
						ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
						ggplot2::ggsave(paste("~/Analyze4C/plots/",sp1[[1]][1],"_transFDRs_",cyc,"cycles_",DandT2,"_lineplot.png",sep=""),width=wd,height=ht)
					}		

					#i might want to add an option of creating plots of the data here
					#e.g.: print(ggplot2::ggplot(data=FDR_cal,ggplot2::aes(x=CO,y=chr1,group=window_RE)) + ggplot2::geom_line(ggplot2::aes(color=window_RE)) + ggplot2::geom_point())
					#or: print(ggplot2::ggplot(data=FDR_cal,ggplot2::aes(x=CO,y=chr1,group=window_RE)) + ggplot2::geom_line(ggplot2::aes(color=window_RE)) + ggplot2::geom_point())	
				}	
			},
			{#8) FDR_calculator_nokb
			
				#getting the smallest window size in RE sites
				small.RE <- readline(prompt=cat("\nenter the smallest window size in RE sites:\n"))
				small.RE <- as.numeric(small.RE)
				#getting the biggest window size in RE sites
				big.RE <- readline(prompt=cat("\nenter the largest window size in RE sites:\n"))	
				big.RE <- as.numeric(big.RE)				
				#getting the step for the window in RE sites
				step.RE	<- readline(prompt=cat("\nenter the step of the growing window sizes in RE sites:\n"))
				step.RE <- as.numeric(step.RE)
				#getting the smallest cutoff size
				small.CO <- readline(prompt=cat("\nenter the smallest cutoff:\n"))
				small.CO <- as.numeric(small.CO)
				#getting the biggest cutoff size
				big.CO <- readline(prompt=cat("\nenter the largest cutoff:\n"))	
				big.CO <- as.numeric(big.CO)
				#getting the step for the cutoff
				step.CO	<- readline(prompt=cat("\nenter the step of the growing cutoff:\n"))
				step.CO <- as.numeric(step.CO)
				#get the number of cycles from user
				cyc <- readline(prompt=cat("\nplease enter the number of cycles you would like to have:\n"))
				cyc <- as.numeric(cyc)
				#calculating all the data for the FDR per window sizes and cutoffs
				cat("\ncalculating the cutoffs and FDRs\nplease wait...\n\n")
				FDR_cal <- FDR_calculator_nokb(rawData,small.RE,big.RE,step.RE,small.CO,big.CO,step.CO,vp.chrom,vp.pos,erase,cyc)
				DandT1 <- toString(Sys.time())
				DandT2 <- gsub(" ","_",DandT1)
				DandT2 <- gsub(":","",DandT2)
				#creating a table file of the results
				write.csv(FDR_cal,paste("~/Analyze4C/FDR/",sp1[[1]][1],"_FDRs_",cyc,"cycles_",DandT2,".csv",sep=""),row.names=FALSE)
				#creating a plot for the trans FDRs					
				cat("\ncreating average FDRs for trans plot\nplease wait...\n\n")
				#print(ggplot2::ggplot(data=FDR_cal,ggplot2::aes(x=window_RE,y=fdr.avg.trans,group=CO)) + ggplot2::geom_line(ggplot2::aes(color=CO)) + ggplot2::geom_point() + ggplot2::ggtitle(paste("window size (by RE sites) VS. trans FDRs, per cutoff percentage -",sp1[[1]][1])) + ggplot2::labs(x="window size per RE sites",y="average FDR of trans chromosomes") + ggplot2::theme(plot.title = ggplot2::element_text(size=11)))
				print(ggplot2::ggplot(data=FDR_cal,ggplot2::aes(x=window_RE,y=fdr.avg.trans,group=CO)) + ggplot2::geom_line(ggplot2::aes(color=factor(CO))) + ggplot2::geom_point(size=0.5) + ggplot2::ggtitle(paste("window size (by RE sites) VS. trans FDRs, per cutoff percentage -",sp1[[1]][1])) + ggplot2::labs(x="window size per RE sites",y="average FDR of trans chromosomes") + ggplot2::theme(plot.title = ggplot2::element_text(size=11)))
				ans <- readline(prompt=cat("\nwould you like to save the plot?\ny/n\n\n"))
				if(ans == "y")
				{
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",sp1[[1]][1],"_transFDRs_",cyc,"cycles_",DandT2,"_lineplot.png",sep=""),width=wd,height=ht)
				}		
				
				#i might want to add an option of creating plots of the data here
				#e.g.: print(ggplot2::ggplot(data=FDR_cal,ggplot2::aes(x=CO,y=chr1,group=window_RE)) + ggplot2::geom_line(ggplot2::aes(color=window_RE)) + ggplot2::geom_point())
				#or: print(ggplot2::ggplot(data=FDR_cal,ggplot2::aes(x=CO,y=chr1,group=window_RE)) + ggplot2::geom_line(ggplot2::aes(color=window_RE)) + ggplot2::geom_point())				
			},
			{#9) use a different raw data file
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
				}
				
				#checking if the raw data taken has its indices in sequence (they might not be when rearranged, for translocations for example)
				#this is in order to see if the user can calculate p-scores using windows with bp size
				inp <- 0 #gets a default of 0
				num_of_chr <- length(unique(rawData[,1])); #getting the number of chromosomes
				for(k in 1:num_of_chr)
				{
					if(is.unsorted(rawData[rawData[,1]==k,2]))
					{
						cat("\nthe raw data seems to be rearranged.\nyou may only do actions that when calculating p-scores it does so only with windows by RE sites (and not bp).\n\n")
						inp <- 1
						break
					}
				}
			},
			{#10) exit this menu
				flag0 <- 0;
			}
		)
	}
}
