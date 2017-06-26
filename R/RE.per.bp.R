#' @export

#this function calculates the average RE sites per bp, or in other words every how many bp are there an RE site
#using a requested raw data sgr file we create a bed file, then a window size is decided by user, the genome is divided into that window size, a tab is added at the end of each line,
#the window size file and the bed file are intersected, giving us the number of RE sites per window
#the function then averages the number of RE sites per window in each chromosome and for the whole genome, and calculates in every how many bp is there an RE site (using the average from the whole genome and not each chromosome separately)

RE.per.bp <- function(file.name=0,choice=1)
{
#converting the raw data sgr file into a bed file. meaning that an extra index file is added (just added 1 to each first index file and removed the number of reads for each RE site)

	if(file.name==0)
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
	}	
	#a bed file is created
	if(choice == 2)
	{
		system(paste("perl /data/home/lab/akivab/Analyze4C/proxy/sgr_to_bed.pl /data/home/lab/akivab/Analyze4C/rawData/rearranged/",file.name," > ~/Analyze4C/temp/rawBed.bed",sep=""))
	}
	else if(choice == 3)
	{
		system(paste("perl /data/home/lab/akivab/Analyze4C/proxy/sgr_to_bed.pl /data/home/lab/akivab/Analyze4C/rawData/coverage_removed/",file.name," > ~/Analyze4C/temp/rawBed.bed",sep=""))
	}
	else
	{
		system(paste("perl /data/home/lab/akivab/Analyze4C/proxy/sgr_to_bed.pl /data/home/lab/akivab/Analyze4C/rawData/original/",file.name," > ~/Analyze4C/temp/rawBed.bed",sep=""))
	}	

############################################################################################################################################################################

#creating a file which divides the whole genome	into windows

	#getting the window size from user
	win.size <- readline(prompt="\nin order to calculate the average RE site per bp, the genome needs to be divided into windows.\nenter a window size in bp.\ne.g. 10000\n")
	win.size <- as.integer(win.size)
	
	#getting the name of the chromosome sizes in genome file the user wants to use
	repeat
	{
		ls_files2 <- system("ls ~/Analyze4C/genomes",intern=TRUE)
		file.name2 <- readline(prompt=cat("These are the chromosome sizes in genome files available\nwhich would you like to use?\n",ls_files2,"\n",sep="\n"))
		ind_files2 <- pmatch(file.name2,ls_files2,nomatch=0)
		if(ind_files2 == 0)
		{
			cat("no such file exists.\nplease try again.\n\n")
		}
		else
		{
			break
		}
	}
	
	#this removes 'chr' from the names of each chromosomes and only leaves the chromosome numbers
	#this is the format that the raw data files and bed files have
	system(paste("sed 's/^chr//g' ~/Analyze4C/genomes/",file.name2," > ~/Analyze4C/temp/genome_temp.txt",sep=""))
	
	#make windows
	system(paste("bedtools makewindows -g ~/Analyze4C/temp/genome_temp.txt -w ",win.size," > ~/Analyze4C/temp/windows.bed",sep=""))
	
	#adding a tab at the end of each line
	#system("sed -i 's/$/\t/' ~/Analyze4C/temp/windows.bed")
	
############################################################################################################################################################################

#intersecting the raw data bed file and the window file
	
	system("bedtools intersect -a ~/Analyze4C/temp/windows.bed -b ~/Analyze4C/temp/rawBed.bed -c > ~/Analyze4C/temp/RE_per_kb.bed")

############################################################################################################################################################################

#calculating the average RE sites in each window, for each chromosome separately and for the whole genome. this is the average number of RE sites per bp
	RE_per_kb <- read.delim("~/Analyze4C/temp/RE_per_kb.bed",header=FALSE,quote="")
	chroms <- read.delim(paste("~/Analyze4C/genomes/",file.name2,sep=""),header=FALSE,quote="")
	size <- nrow(chroms)
	
	cat("\nthe average number of RE sites per ",win.size," bp are:\n\n",sep="")
	for(i in 1:size)
	{
		mn <- round(mean(RE_per_kb[RE_per_kb[,1]==i,4]))
		cat(paste(chroms[i,1]),": ",mn,"\n",sep="")
	}
	
	mn2 <- round(mean(RE_per_kb[,4]))
	cat("all: ",mn2,"\n",sep="")

#calculating the average number of bp between every RE site
	
	final <- round(win.size/mn2)
	cat("\nif we divide the window size by the average number of RE sites per window size (for the whole genome, not each chromosome separately)\nwe get roughly 1 RE site every ",final," bps\n",sep="")
	
############################################################################################################################################################################

#removing the files that aren't needed anymore

	system("rm ~/Analyze4C/temp/windows.bed")
	system("rm ~/Analyze4C/temp/rawBed.bed")
	system("rm ~/Analyze4C/temp/RE_per_kb.bed")
	system("rm ~/Analyze4C/temp/genome_temp.txt")
}
