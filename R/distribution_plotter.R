#' @export

#this function takes either raw data or p-score data and creates different plots according to the users request
#firstly the user chooses if to take a raw data or a p-score file
#then the function gets different data about the file, including name, cis etc.
#the function creates transliterations of the chromosome numbers, this is for later use in the plots
#the user then chooses a plot to do:
#boxplots, violin plots, density plots, and bar charts: the idea is to compare the spread of reads or p-scores between chromosomes
#One of the uses for this is to prove that when I remove coverage I can do so from all of trans and by chromosome. 
#The code I wrote that removes coverage looks at the reads on each RE site and according to that gives the site a score, the more reads the higher the score. 
#And then we remove coverage, the lower the score - the higher chance it will be removed. 
#Since we are working with probabilities and scores then we need to make sure that there isn't any specific chromosome that will influence the scores more than any other. 
#If the distributions of read numbers for RE sites are similar in all chromosomes (in trans) then we can include all in the same calculation of scoring when  removing coverage. 
#normal distribution: the idea is to see if the data distributes normally. an example of when this is important is in order to know what type of correlation to use.
#if they distribute normally then you can use pearson, if not then only spearman.
#the way the plot works is by plotting the data against a known data which distributes normally, if our data is also normal then they should correlate with a line going diagonally
#in the middle.
#the minumum number of reads or p-score by default is 0 but could be changed.
#the file could always be changed.

#note: when the minimum reads is 0 and the option of transforming the scale by log10, then all the values of 0 will be turned into infinity and cannot be used
#these values will be removed and a warning will be created (which is fine). in order to prevent this, either don't transform the axis, or use a minimum of above 0
#also, currently when a plot is made it prints it to screen, this could be time consuming. i could maybe ask user if they want this or if just to save the plot without printing to screen
distribution_plotter <- function(Experiments_4C,dat=0,file.name="NA",tp1=0)
{
	cat("\ndistribution plotter:\n\nhere you choose a type of plot for distribution viewing purposes.\nthe boxplots,violin plots,density plots,and bar charts are meant to show the similarities or differences of distribution between chromosomes.\nthe normal distribution plot shows if data distributes normally\nthiscurrently works only for trans\n\n")
	flag <- 1
	while(flag == 1)
	{
		if(tp1==0) #if the original input has data and file name, the user should also input tp1 as 1 or 2 (if raw data or p-score data respectivally)
		{
			tp1 <- as.integer(readline(prompt=cat("\nwhat type of file would you like to use (enter the number):\n\n1) raw data\n2) p-score\n\n")))
			if(tp1 == 1) #choosing a raw data file
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
				
				tmp_nm <- sp1[1]
				
				#the default minimum number of reads is 0
				min_reads <- 0
				min_any <- min_reads
			}
			else #choosing a p-score file
			{
				tp2 <- as.integer(readline(prompt=cat("\nchoose the name of the folder you would like to use a p-score file from (choose the number):\n\n1) no_bp\n2) with_bp\n\n")))
				if(tp2 == 1)
				{
					fld <- "no_bp"
				}
				else
				{
					fld <- "with_bp"
				}
				#getting the name of the p score file the user wants to use
				repeat
				{
					ps_files <- system(paste("ls ~/Analyze4C/pScores/",fld,sep=""),intern=TRUE)
					ps.file <- readline(prompt=cat("These are the p-score files available\nwhich would you like to use?\n",ps_files,"\n",sep="\n"))
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
				dat <- read.table(paste("~/Analyze4C/pScores/",fld,"/",ps.file,sep=""))

				#finding the specific raw data files information in Experiments_4C
				sp3 <- unlist(strsplit(ps.file,"[.]"))
				sp4 <- strsplit(sp3[[1]][1],"_")
				
				#getting the size of the window per RE sites
				sp5 <- sp4[[1]][4]
				RE_wind <- strsplit(sp5[[1]][1],"RE")
				RE_wind <- as.numeric(RE_wind)
			
				if(tp2 == 2)
				{
					#getting the size of the window per bp sites
					sp6 <- sp4[[1]][5]
					bp_wind <- strsplit(sp6[[1]][1],"bp")
					bp_wind <- as.numeric(bp_wind)
				}
				
				#finding the specific raw data files information in Experiments_4C			
				out <- findIn.Experiments_4C(sp4[[1]][1],sp4[[1]][2],sp4[[1]][3],Experiments_4C)
				
				#getting the cis chromosome number
				cis <- as.numeric(out[2])	
		
				tmp_nm <- sp3[1]
				
				#the default minimum p-score is 0
				min_pScore <- 0 
				min_any <- min_pScore
			}
		}
		else if(tp1 == 1) #getting the raw data datas information
		{				
			#finding the specific raw data files information in Experiments_4C
			sp1 <- unlist(strsplit(file.name,"[.]"))
			sp2 <- strsplit(sp1[[1]][1],"_")
			out <- findIn.Experiments_4C(sp2[[1]][1],sp2[[1]][2],sp2[[1]][3],Experiments_4C)
			
			#getting the cis chromosome number
			cis <- as.numeric(out[2])	
			
			tmp_nm <- sp1[1] 
			
			#the default minimum number of reads is 0
			min_reads <- 0
			min_any <- min_reads	
		}
		else if(tp1 == 2) #getting the p-scores datas information
		{
			#finding the specific raw data files information in Experiments_4C
			sp3 <- unlist(strsplit(ps.file,"[.]"))
			sp4 <- strsplit(sp3[[1]][1],"_")
		
			#finding the specific raw data files information in Experiments_4C			
			out <- findIn.Experiments_4C(sp4[[1]][1],sp4[[1]][2],sp4[[1]][3],Experiments_4C)
			
			#getting the cis chromosome number
			cis <- as.numeric(out[2])	
	
			tmp_nm <- sp3[1]
			
			#the default minimum p-score is 0
			min_pScore <- 0 
			min_any <- min_pScore
		}

		#changing the column names
		names(dat) <- c("chr","ind","reads")
		
		#putting the transliteration of each chromosome number in the position of the number
		num_names <- c("one","two","three","four","five","six","seven","eight","nine","ten","eleven","twelve","thirteen","fourteen","fifteen","sixteen","seventeen","eighteen","nineteen","twenty","twentyOne","twentyTwo","twentyThree","twentyFour","twentyFive","twentySix","twentySeven","twentyEight","twentyNine","Thirty")
		
		num_of_chroms <- length(unique(dat[,1]))
		for(k in 1:num_of_chroms)
		{
			dat[dat[,1]==k,1] <- num_names[k]
		}
		
		#getting the translitaration of the cis number
		cis_name <- num_names[cis]
	
		#the default of the data is all chromosomes
		sections <- "all"	
		
		#the default of the data is all the RE sites
		dat_temp <- dat[dat[,3]>=min_any,]
		
		#choice2 represents the number from the menu chosen by user
		choice2 <- 0
		#the while loop stops when 8 (choose a different file) or 9 (exit) are chosen
		while((choice2 != 9) && (choice2 != 10))
		{	
			#printing to screen the current file, the chromosomes looked at, and the current minimum number of reads/p-score chosen
			cat("\n\nthe current file chosen is ",tmp_nm,".sgr\n",sep="")
			
			if(sections=="all")
			{
				cat("\nthe sections looked at are all the chromosomes\n")
			}
			else if(sections=="cis")
			{
				cat("\nthe section looked at is the cis chromosome - chromosome",cis,"\n")
			}
			else if(sections=="trans")
			{
				cat("\nthe sections looked at are all the trans chromosomes\n")
			}
			else
			{
				cat("\nthe section looked at is",sections,"\n")
			}
			
			if(tp1 == 1)
			{
				cat("\nthe RE sites that have minimum ",min_reads,"reads are being used\n")
				#min_any <- min_reads
			}
			else
			{
				cat("\nthe RE sites that have the minimum p-score of",min_pScore,"are being used\n")
				#min_any <- min_pScore
			}
						
			#choosing what to do by user
			choice2 <- as.integer(readline(prompt=cat("\nchoose the type of plot you would like to produce (enter the number of the option):\n\n1) coverage against minimum reads line plot\n2) boxplots\n3) violin plots\n4) density plots\n5) bar chart\n6) test for normal distribution\n7) choose a different minimum for each RE site\n8) choose different sections\n9) choose a different file\n10) exit\n\n")))

			#coverage vs. minimum reads line plot
			if(choice2 == 1)
			{
				if(min_any > 0)
				{
					ans0 <- readline(prompt=cat("\nwarning:\nthe minimum number of reads/p-score is more than 0\nthis function should be done when the minimum is 0\nwould you like to procceed anyhow?\ny/n\n\n"))
					if(ans0 == "y")
					{
						if(tp1 == 1)
						{
							name <- file.name
						}
						else if(tp1 == 2)
						{
							name <- ps.file
						}
						cat("\nplot of minimum reads against coverages:\n\n")
						coverageVSmin_linePlot(dat_temp,name,cis_name)
					}
				}
				else
				{
					if(tp1 == 1)
					{
						name <- file.name
					}
					else if(tp1 == 2)
					{
						name <- ps.file
					}
					cat("\nplot of minimum reads against coverages:\n\n")
					coverageVSmin_linePlot(dat_temp,name,cis_name)
				}
			}
			#boxplots
			else if(choice2 == 2)
			{
				ans1 <- readline(prompt=cat("\nwould you like the y axis to be a log10 axis?\ny/n\n\n"))
				if(ans1 == "y")
				{
					tt <- paste(tmp_nm,"_",sections,"_Boxplots_log10",sep="")
					print(ggplot2::ggplot(dat_temp, ggplot2::aes(x=chr,y=reads,fill=chr)) + ggplot2::geom_boxplot() + ggplot2::scale_y_log10() + ggplot2::ggtitle(tt))
					if(min_any > 0)
					{
						nm <- paste(tmp_nm,"_",sections,"_Boxplots_log10_minimum",min_any,".png",sep="")
					}
					else
					{
						nm <- paste(tmp_nm,"_",sections,"_Boxplots_log10.png",sep="")
					}
				}
				else
				{
					tt <- paste(tmp_nm,"_",sections,"_Boxplots",sep="")
					print(ggplot2::ggplot(dat_temp, ggplot2::aes(x=chr,y=reads,fill=chr)) + ggplot2::geom_boxplot() + ggplot2::ggtitle(tt))
					if(min_any > 0)
					{
						nm <- paste(tmp_nm,"_",sections,"_Boxplots_minimum",min_any,".png",sep="")
					}
					else
					{
						nm <- paste(tmp_nm,"_",sections,"_Boxplots.png",sep="")
					}
				}
				
				ans2 <- readline(prompt=cat("\nwould you like to save the plot?\ny/n\n\n"))
				if(ans2 == "y")
				{
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
				}
			}
			#violin plots
			else if(choice2 == 3)
			{
				ans3 <- readline(prompt=cat("\nwould you like the y axis to be a log10 axis?\ny/n\n\n"))
				if(ans3 == "y")
				{
					tt <- paste(tmp_nm,"_",sections,"_Violinplots_log10",sep="")
					print(ggplot2::ggplot(dat_temp, ggplot2::aes(x=chr,y=reads,fill=chr)) + ggplot2::geom_violin() + ggplot2::scale_y_log10() + ggplot2::ggtitle(tt))
					if(min_any > 0)
					{
						nm <- paste(tmp_nm,"_",sections,"_Violinplots_log10_minimum",min_any,".png",sep="")
					}
					else
					{
						nm <- paste(tmp_nm,"_",sections,"_Violinplots_log10.png",sep="")
					}
				}
				else
				{
					tt <- paste(tmp_nm,"_",sections,"_Violinplots",sep="")
					print(ggplot2::ggplot(dat_temp, ggplot2::aes(x=chr,y=reads,fill=chr)) + ggplot2::geom_violin() + ggplot2::ggtitle(tt))
					if(min_any > 0)
					{
						nm <- paste(tmp_nm,"_",sections,"_Violinplots_minimum",min_any,".png",sep="")
					}
					else
					{
						nm <- paste(tmp_nm,"_",sections,"_Violinplots.png",sep="")
					}
				}
				
				ans4 <- readline(prompt=cat("\nwould you like to save the plot?\ny/n\n\n"))
				if(ans4 == "y")
				{
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
				}
			}
			#density plots
			else if(choice2 == 4)
			{
				tt <- paste(tmp_nm,"_",sections,"_densityplots",sep="")
				print(ggplot2::ggplot(dat_temp, ggplot2::aes(x=reads,fill=chr)) + ggplot2::geom_density() + ggplot2::facet_grid(chr ~ .) + ggplot2::ggtitle(tt))
				if(min_any > 0)
				{
					nm <- paste(tmp_nm,"_",sections,"_densityplots_minimum",min_any,".png",sep="")
				}
				else
				{
					nm <- paste(tmp_nm,"_",sections,"_densityplots.png",sep="")
				}
		
				ans5 <- readline(prompt=cat("\nwould you like to save the plot?\ny/n\n\n"))
				if(ans5 == "y")
				{
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
				}
			}
			#bar charts
			else if(choice2 == 5)
			{
				tt <- paste(tmp_nm,"_",sections,"_barcharts",sep="")
				print(ggplot2::ggplot(dat_temp, ggplot2::aes(x=reads,fill=chr)) + ggplot2::geom_bar() + ggplot2::ggtitle(tt))
				if(min_any > 0)
				{
					nm <- paste(tmp_nm,"_",sections,"_barcharts_minimum",min_any,".png",sep="")
				}
				else
				{
					nm <- paste(tmp_nm,"_",sections,"_barcharts.png",sep="")
				}
			
				ans6 <- readline(prompt=cat("\nwould you like to save the plot?\ny/n\n\n"))
				if(ans6 == "y")
				{
					wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
					ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
					ggplot2::ggsave(paste("~/Analyze4C/plots/",nm,sep=""),width=wd,height=ht)
				}
			}
			#normal distribution
			else if(choice2 == 6)
			{
				#cat("\nnote: the plot here includes the minimum number of reads/p-score, as opposed to the other plots which test what is above this number\n\n")
				plot.new()
				if(sections=="all")
				{
					qqnorm(dat[dat[,3]>=min_any,3],main=paste(tmp_nm))
					qqline(dat[dat[,3]>=min_any,3])
				}
				else if(sections=="cis")
				{
					qqnorm(dat[dat[,1]==cis_name & dat[,3]>=min_any,3],main=paste(tmp_nm))
					qqline(dat[dat[,1]==cis_name & dat[,3]>=min_any,3])
				}
				else if(sections=="trans")
				{
					qqnorm(dat[dat[,1]!=cis_name & dat[,3]>=min_any,3],main=paste(tmp_nm))
					qqline(dat[dat[,1]!=cis_name & dat[,3]>=min_any,3])
				}
				else
				{
					qqnorm(dat[dat[,1]==sections_temp & dat[,3]>=min_any,3],main=paste(tmp_nm))
					qqline(dat[dat[,1]==sections_temp & dat[,3]>=min_any,3])
				}
								
				ans7 <- readline(prompt=cat("\nwould you like to save the plot?\ny/n\n\n"))
				if(ans7 == "y")
				{
					if(min_any > 0)
					{
						savePlot(filename=paste("~/Analyze4C/plots/",tmp_nm,"_",sections,"_normalDistribution_tester_minimum",min_any,".png",sep=""),type="png")
					}
					else
					{
						savePlot(filename=paste("~/Analyze4C/plots/",tmp_nm,"_",sections,"_normalDistribution_tester.png",sep=""),type="png")
					}
				}
			}
			#change minimum for each RE site
			else if(choice2 == 7)
			{
				if(tp1 == 1)
				{
					cat("\nthe RE sites that have",min_reads,"reads and above are being used\n")
					min_reads <- as.integer(readline(prompt=cat("\nenter the new value:\n\n")))
					min_any <- min_reads
				}
				else
				{
					cat("\nthe RE sites that have a",min_pScore,"p-score and above are being used\n")
					min_pScore <- as.numeric(readline(prompt=cat("\nenter the new value:\n\n")))
					min_any <- min_pScore
				}

				if(sections=="all")
				{
					dat_temp <- dat[dat[,3]>=min_any,]
				}
				else if(sections=="cis")
				{
					dat_temp <- dat[dat[,1]==cis_name & dat[,3]>=min_any,]				
				}
				else if(sections=="trans")
				{
					dat_temp <- dat[dat[,1]!=cis_name & dat[,3]>=min_any,]
				}
				else
				{
					dat_temp <- dat[dat[,1]==sections & dat[,3]>=min_any,]
				}				
			}
			#choose different sections to look at
			else if(choice2 == 8)
			{
				ans8 <- as.integer(readline(prompt=cat("\nchoose which chromosomes to look at:\n1) all\n2) trans\n3) cis\n4) specific chromosome\n\n")))
				if(ans8 == 1)
				{
					sections <- "all"
					dat_temp <- dat[dat[,3]>=min_any,]
				}
				else if(ans8 == 2)
				{
					sections <- "trans"
					dat_temp <- dat[dat[,1]!=cis_name & dat[,3]>=min_any,]
				}
				else if(ans8 == 3)
				{
					sections <- "cis"
					dat_temp <- dat[dat[,1]==cis_name & dat[,3]>=min_any,]
				}
				else if(ans8 == 4)
				{
					sections <- readline(prompt=cat("\nenter the chromosome number you would like to use (write the number in english):\n\n"))
					if(sections==cis_name)
					{
						sections <- "cis"
						dat_temp <- dat[dat[,1]==cis_name & dat[,3]>=min_any,]
					}
					else
					{
						dat_temp <- dat[dat[,1]==sections & dat[,3]>=min_any,]
					}
					sections_temp <- sections
					sections <- cat("chr",as.integer(grep(paste("^",sections,"$",sep=""),num_names)))
				}
			}
			else if(choice2 == 9)
			{
				tp1 <- 0
			}
			else if(choice2 == 10)
			{
				flag <- 0
			}
		}
	}	
}
