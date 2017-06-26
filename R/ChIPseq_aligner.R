#' @export

#note: if the function 'fastq-dump' isn't available but the package of sratoolkit is installed, you could create a symbolic link which lets you use the function
#e.g.: "ln -s /data/private/Software/sratoolkit.2.2.2b-centos_linux64/bin/fastq-dump /data/home/lab/akivab/bin"

#the parameters that are important:
#1) the tags - how many reads align to each region
#2) summits - the position that has the highest fragment pileup for each peak - i might want to use the peak to comapre since this has the position which most represents the peak
#3) -10*LOG10(pvalue) - the higher the number the stronger the peak result

ChIPseq_aligner <- function(ChIPseq_data,MACS_outputs)
{
	#ask if fastq-dump and MACS are available, if not then the function will not respond
	ans1 <- readline(prompt=cat("\nare the functions and package 'fastq-dump'(which is part of the sra toolkit package) and 'MACS' installed?\ny/n\n\n"))
	if(ans1 == "y")
	{
		#printing all the available sra and fastq data files and descriptions. if there aren't any or user wants new data, they can download it.
		flag <- 0 #this tells the function if the user will download new data
		sra_or_fastq_files <- system("ls ~/Analyze4C/ChIPseq/SRA_or_FASTQ",intern=TRUE)
		if(length(sra_or_fastq_files) != 0)
		{
			repeat
			{
				cat("\nthese are the sra and fastq data files of ChIPseq that should be available:\n\n")
				print(ChIPseq_data)
				ans2 <- as.integer(readline(prompt=cat("\nwould you like to:\n1) use data from the list\n2) download new ChIPseq data to server\n\n")))
				if(ans2 == 1) #using already existing data
				{
						sra_or_fastq_filename <- readline(prompt=cat("\nThese are the sra and fastq files of ChIPseq data available\nwhich would you like to use?\n",sra_or_fastq_files,"\n\n"))
						sra_or_fastq_Ind_files <- pmatch(sra_or_fastq_filename,sra_or_fastq_files,nomatch=0)
						if(sra_or_fastq_Ind_files == 0)
						{
							cat("no such file exists.\nplease try again.\n\n")
						}
						else
						{
							chip_ind <- which(ChIPseq_data$Name_of_File %in% sra_or_fastq_filename)
							DandT1 <- ChIPseq_data$Download_Date[chip_ind]
							DandT2 <- gsub(" ","_",DandT1)
							DandT2 <- gsub(":","",DandT2)
							break
						}
				}
				else if(ans2 == 2)
				{
					flag <- 1
					break
				}
			}	
		}	
		else
		{
			cat("\nthe folder is empty. It seems that there are no sra or fastq files to choose from.\nIn order to proceed you must download ChIPseq data in sra or fastq formats\n\n")
			ans3 <- readline(prompt=cat("\nwould you like to download new ChIPseq data?\ny/n\n\n"))
			if(ans3 ==  "y") #downloading new data
			{
				flag <- 1
			}	
			else
			{
				return()
			}
		}

		#downloading new data
		if(flag == 1)
		{
			cat("\nChIPseq details:\n\n")
			#name of organism
			org <- readline(prompt=cat("enter the name of the ChIPseqs organism:\n"))
			#name of tissue
			tissue <- readline(prompt=cat("enter the name of the ChIPseqs tissue:\n"))
			#age
			age <- readline(prompt=cat("enter the age of the cells when ChIPseq was extracted:\n"))
			#additional notes
			ans4 <- readline(prompt=cat("\nwould you like to add any additional information?\ny/n\n\n"))
			if(ans4 == "y")
			{
				note <- readline(prompt=cat("write additional information about the ChIPseq data:\n"))
			}
			else
			{
				note <- NA
			}
			
			#getting the link to the file
			link_sra_or_fastq <- readline(prompt=cat("\nenter the link to the sra or fastq data file:\n\n"))
			
			#checking that the 'temp' folder is empty, if not then empty it out
			ls_temp <- system("ls ~/Analyze4C/ChIPseq/temp/",intern=TRUE)
			#removing sra and fastq files from temp
			if(identical(grep(".sra",ls_temp),integer(0)) == "FALSE")
			{
				system("rm ~/Analyze4C/ChIPseq/temp/*.sra")
			}
			#removing fastq files from temp
			if(identical(grep(".fastq",ls_temp),integer(0)) == "FALSE")
			{
				system("rm ~/Analyze4C/ChIPseq/temp/*.fastq")
			}
			
			#getting the sra or fastq file into 'temp'
			system(paste("wget -P ~/Analyze4C/ChIPseq/temp",link_sra_or_fastq))
			
			#getting the date and time of downloading ChIPseq data
			DandT1 <- toString(Sys.time())
			DandT2 <- gsub(" ","_",DandT1)
			DandT2 <- gsub(":","",DandT2)

			#getting the name of the sra or fastq file
			sra_or_fastq_filename <- system("ls ~/Analyze4C/ChIPseq/temp | egrep '.sra$|.fastq$'",intern=TRUE)
			
			#moving the sra or fastq file from 'temp' to the 'SRA_or_FASTQ' folder
			system(paste("mv ~/Analyze4C/ChIPseq/temp/",sra_or_fastq_filename," ~/Analyze4C/ChIPseq/SRA_or_FASTQ",sep=""))
			
			#adding all the details from above to the file ChIPseq_data.txt
			ChIPseq_data[nrow(ChIPseq_data)+1,] <- c(sra_or_fastq_filename,org,tissue,age,link_sra_or_fastq,note,DandT1)
			#sorting the list of sra and fastq files by name (and sorting the row indices)
			ChIPseq_data <- ChIPseq_data[order(ChIPseq_data$Name_of_File),]
			rownames(ChIPseq_data) <- seq(length=nrow(ChIPseq_data))
			#adding the new data to the file (by erasing the old one and creating a new one)
			system("rm ChIPseq_data.txt")
			write.table(ChIPseq_data,"ChIPseq_data.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
			
			ans5 <- readline(prompt=cat("\nwould you like to continue with this current file that was just downloaded?\ny/n\n\n"))
			if(ans5 == "n")
			{
				stop
			}
		}
		
		sra_or_fastq_filename_temp <- strsplit(sra_or_fastq_filename,"[.]")
		prefix_name <- paste(sra_or_fastq_filename_temp[[1]][1],"_",DandT2,sep="")
		
		#converting sra to fastq if an sra file was chosen
		if(identical(grep(".fastq$",sra_or_fastq_filename),integer(0)))
		{
			#fastq dump the sra data, send the output to a temp folder
			cat("\n\nfastq dumping:\nretrieve fastq data from sra file\nplease wait...\n\n")
			system(paste("fastq-dump -O ~/Analyze4C/ChIPseq/temp/ ~/Analyze4C/ChIPseq/SRA_or_FASTQ/",sra_or_fastq_filename,sep=""))
					
			#getting the name of the fastq file in temp
			fastq_filename <- system("ls ~/Analyze4C/ChIPseq/temp | grep '.fastq$'",intern=TRUE)
		}
		else
		{
			fastq_filename <- sra_or_fastq_filename
		}
		
		#get the links to the genome files folder (the folder should be in TAIR10, it has a bunch of files with 'genome' and extensions of 'ebwt' and 'fa') and the genes annotation folder and file (genes.gtf)
		link_genome <- readline(prompt=cat("\nplease provide a link to the whole genomes bowties index file and add '/genome' at the end, this should be in the bowtie index folder of the genomes sequence folder\ne.g.: Arabidopsis/Arabidopsis_thaliana_TAIR10/Sequence/BowtieIndex/genome\n\n"))

		#bowtie
		cat("\nthe fastq data is being aligned to the genome via bowtie\nplease wait...\n\n")
		system(paste("bowtie -m 1 -q -S ",link_genome," ~/Analyze4C/ChIPseq/temp/",fastq_filename," ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam",sep=""))

		ans8 <- as.integer(readline(prompt=cat("\nchoose the genome type you are working with:\n1) Human\n2) c.elegans\n3) Mouse\n4) Drosophila\n5) other\n\n")))
		if(ans8 == 1)
		{
			cat("\nthe mappable size of the Human genome that appears in MACS is 2.7e9\n")
			genome_size_name <- "hs"
		}
		else if(ans8 == 2)
		{
			cat("\nthe mappable size of the c.elegans genome that appears in MACS is 9e7\n")
			genome_size_name <- "ce"
		}
		else if(ans8 == 3)
		{
			cat("\nthe mappable size of the Mouse genome that appears in MACS is 1.87e9\n")
			genome_size_name <- "mm"
		}
		else if(ans8 == 4)
		{
			cat("\nthe mappable size of the Drosophila genome that appears in MACS is 1.2e8\n")
			genome_size_name <- "dm"
		}

		if(ans8 == 5)
		{
			flag_size <- 1 #this tells us if MACS will be using an internal genome size or a custom one
			#choosing a genome size file
			ls_genomes <- system("ls ~/Analyze4C/genomes/",intern=TRUE)
			genome_file.name <- readline(prompt=cat("\nchoose the file with genome sizes that you would like to use:\n",ls_genomes,"",sep="\n"))
			genome_sizes <- read.table(paste("~/Analyze4C/genomes/",genome_file.name,sep=""))

			numOFchroms <- nrow(genome_sizes)
			genome_size <- 0
			for(t in 1:numOFchroms)
			{
				genome_size <- genome_size + genome_sizes[t,2]
			}
			cat("\nthe number of bp reached from summing all the chromosome sizes from",genome_file.name,"is",genome_size,"\n\n")

			ans10 <- as.integer(readline("\nwould you like to:\n1) use this size\n2) enter a different size\n\n"))
			if(ans10 == 2)
			{
				genome_size <- readline(prompt=cat("\nenter the mappable genome size (the altogether size that is non-repetitive and is mappable):\n\n"))
			}
		}
		else
		{
			flag_size <- 0 #this tells us if MACS will be using an internal genome size or a custom one
			ans9 <- as.integer(readline(prompt=cat("\nwould you like to:\n1) continue with this size\n2) sum the sizes from a genome file\n3) enter a size manually\n\n")))
			if(ans9 == 1)
			{
				flag_size <- 0 #this tells us if MACS will be using an internal genome size or a custom one
			}
			else if(ans9 == 2)
			{
				flag_size <- 1 #this tells us if MACS will be using an internal genome size or a custom one
				#choosing a genome size file
				ls_genomes <- system("ls ~/Analyze4C/genomes/",intern=TRUE)
				genome_file.name <- readline(prompt=cat("\nchoose the file with genome sizes that you would like to use:\n",ls_genomes,"",sep="\n"))
				genome_sizes <- read.table(paste("~/Analyze4C/genomes/",genome_file.name,sep=""))

				numOFchroms <- nrow(genome_sizes)
				genome_size <- 0
				for(t in 1:numOFchroms)
				{
					genome_size <- genome_size + genome_sizes[t,2]
				}
				cat("\nthe number of bp reached from summing all the chromosome sizes from",genome_file.name,"is",genome_sum,"\n\n")
				
				ans10 <- as.integer(readline("\nwould you like to:\n1) use this size\n2) enter a different size\n\n"))
				if(ans10 == 2)
				{
					genome_size <- readline(prompt=cat("\nenter the mappable genome size (the altogether size that is non-repetitive and is mappable):\n\n"))
				}	
			}
			else if(ans9 == 3)
			{
				flag_size <- 1 #this tells us if MACS will be using an internal genome size or a custom one
				genome_size <- readline(prompt=cat("\nenter the mappable genome size (the altogether size that is non-repetitive and is mappable):\n\n"))
			}	
		}

		setwd("~/Analyze4C/ChIPseq/temp/")
		if(flag_size == 0) #if we are using the internal genome size
		{
			#MACS - entering the genome files, and the name of the fastq file created from the fastq dump
			#if user wants to change the default parameters, then they can do so
			cat("\nThe default of the MACS running is as follows:\n\n")
			cat("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"-g",genome_size_name,"--nomodel --nolambda -B -S")
			cat("\n")
			ans6 <- readline(prompt=cat("\nwould you like to change the parameters?\ny/n\n\n"))
			if(ans6 == "y")
			{
				ans7 <- 0
				cat("\nnote: only the parameters are permitted to be changed or added to, not the files or output destinations\n")
				while(ans7 != "y")
				{
					new_parameters <- readline(prompt=cat("\nenter the new parameters freely. Make sure that they are legal:\n\n"))
					cat("\nthis is the MACS command as requested:\n\n")
					cat("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"-g",genome_size_name,new_parameters)
					cat("\n")
					MACS_command <- paste("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"_",DandT2,"-g",genome_size_name,new_parameters) #recording the full command used in order to put in MACS_outputs later
					ans7 <- readline(prompt=cat("\n\nwould you like to execute this command?\ny/n\n\n"))
				}
				system(paste("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"--nomodel --nolambda -g",genome_size_name,"-B -S"))
			}
			else
			{
				MACS_command <- paste("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"-g",genome_size_name,"--nomodel --nolambda -B -S") #recording the full command used in order to put in MACS_outputs later
				system(paste("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"-g",genome_size_name,"--nomodel --nolambda -B -S"))
			}
		}
		else #if we are using a non internal genome size
		{
			#MACS - entering the genome files, and the name of the fastq file created from the fastq dump
			#if user wants to change the default parameters, then they can do so
			cat("\nThe default of the MACS running is as follows:\n\n")
			cat("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"-g",genome_size,"--nomodel --nolambda -B -S")
			cat("\n")
			ans6 <- readline(prompt=cat("\nwould you like to change the parameters?\ny/n\n\n"))
			if(ans6 == "y")
			{
				ans7 <- 0
				cat("\nnote: only the parameters are permitted to be changed or added to, not the files or output destinations\n")
				while(ans7 != "y")
				{
					new_parameters <- readline(prompt=cat("\nenter the new parameters freely. Make sure that they are legal:\n\n"))
					cat("\nthis is the MACS command as requested:\n\n")
					cat("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"-g",genome_size,new_parameters)
					cat("\n")
					MACS_command <- paste("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"-g",genome_size,new_parameters) #recording the full command used in order to put in MACS_outputs later
					ans7 <- readline(prompt=cat("\n\nwould you like to execute this command?\ny/n\n\n"))
				}
				system(paste("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"--nomodel --nolambda -g",genome_size,"-B -S"))
			}
			else
			{
				MACS_command <- paste("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"-g",genome_size,"--nomodel --nolambda -B -S")  #recording the full command used in order to put in MACS_outputs later
				system(paste("macs14 -t ~/Analyze4C/ChIPseq/temp/ChIPseq_bowtie.sam --name",prefix_name,"-g",genome_size,"--nomodel --nolambda -B -S"))
			}			
		}
		setwd("~/Analyze4C/")
		
		#getting all the peaks data from the peaks excel file
		peaks <- read.table(paste("~/Analyze4C/ChIPseq/temp/",prefix_name,"_peaks.xls",sep=""),header=TRUE,stringsAsFactors=FALSE)

		#moving the peaks.xls file and MACS_bedGraph folder in a folder, adding the details of running MACS into 'MACS_outputs'
		cat("\nmoving the peaks.xls file and MACS_bedGraph folder into the 'peaks' folder...\n\n")
		system(paste("mv ~/Analyze4C/ChIPseq/temp/",prefix_name,"_peaks.xls ~/Analyze4C/ChIPseq/peaks",sep=""))
		system(paste("mv ~/Analyze4C/ChIPseq/temp/",prefix_name,"_MACS_bedGraph ~/Analyze4C/ChIPseq/peaks",sep=""))
		
		cat("\nadding the MACS commands used into MACS_outputs.txt...\n\n")
		MACS_outputs[nrow(MACS_outputs)+1,] <- c(sra_or_fastq_filename_temp[[1]][1],DandT1,MACS_command)
		MACS_outputs <- MACS_outputs[order(MACS_outputs$Name),]
		rownames(MACS_outputs) <- seq(length=nrow(MACS_outputs))
		system("rm MACS_outputs.txt")
		write.table(MACS_outputs,"MACS_outputs.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
		
		cat("\ncreating an p-score file, a tag file, and a combined file for the peaks, and placing them in the 'peaks' folder...\n\n")
		ans11 <- readline(prompt=cat("\nwould you like to remove 'Pt' and 'Mt' from the MACS output?\ny/n\n\n"))
		if(ans11 == "y")
		{
			peaks_temp <- peaks[peaks$chr!="Pt" & peaks$chr!="Mt",]
		}
		else
		{
			peaks_temp <- peaks
		}
		
		#creating the combined file
		peaks_tagsANDpscores <- data.frame(peaks_temp$chr,peaks_temp$start,peaks_temp$end,peaks_temp$tags,peaks_temp$X.10.log10.pvalue.)
		#names(peaks_tagsANDpscores) <- c("chr","start","end","tags","X.10.log10.pvalue.")
		write.table(peaks_tagsANDpscores,paste("~/Analyze4C/ChIPseq/peaks/",prefix_name,"_peaks_tagsANDpscores.bed",sep=""),sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
		
		#creating only the tags file
		peaks_tags <- data.frame(peaks_temp$chr,peaks_temp$start,peaks_temp$end,peaks_temp$tags)
		#names(peaks_tags) <- c("chr","start","end","tags")
		write.table(peaks_tags,paste("~/Analyze4C/ChIPseq/peaks/",prefix_name,"_peaks_tags.bedgraph",sep=""),sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)

		#creating only the p-score file
		peaks_pscores <- data.frame(peaks_temp$chr,peaks_temp$start,peaks_temp$end,peaks_temp$X.10.log10.pvalue.)
		#names(peaks_pscores) <- c("chr","start","end","X.10.log10.pvalue.")
		write.table(peaks_pscores,paste("~/Analyze4C/ChIPseq/peaks/",prefix_name,"_peaks_pscores.bedgraph",sep=""),sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
		
		#erasing all the outputs from temp (user could download the relevant files for themselves)
		cat("\nall the output files (except for the files placed in the 'peaks' folder) will be erased")
		cat("\nif you will like to save any for yourself, download the files to your personal computer\n")
		lin <- readline(prompt=cat("\nonce you are done downloading your files -\n\npress any key to continue and erase the files...\n\n"))
		system("rm -rf ~/Analyze4C/ChIPseq/temp/*")
	}
	else
	{
		cat("\ninstall the sra toolkit and MACS packages in order to use the function\n\n")
	}
}
