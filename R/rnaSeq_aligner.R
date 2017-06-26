#' @export

#note: if the function 'fastq-dump' isn't available but the package of sratoolkit is installed, you could create a symbolic link which lets you use the function
#e.g.: "ln -s /data/private/Software/sratoolkit.2.2.2b-centos_linux64/bin/fastq-dump /data/home/lab/akivab/bin"

rnaSeq_aligner <- function(RNAseq_data,tophat_cufflinks_outputs)
{
	#ask if fastq-dump, tophat, and cufflinks are available, if not then the function will not respond
	ans1 <- readline(prompt=cat("\nare the functions and packages 'fastq-dump'(which is part of the sra toolkit package),'tophat', and 'cufflinks' installed?\ny/n\n\n"))
	if(ans1 == "y")
	{
		#printing all the available sra data files and descriptions. if there aren't any or user wants new data, they can download it.
		flag <- 0 #this tells the function if the user will download new data
		sra_files <- system("ls ~/Analyze4C/RNAseq/sra",intern=TRUE)
		if(length(sra_files) != 0)
		{
			repeat
			{
				cat("\nthese are the sra data files of RNA-seq that should be available:\n\n")
				print(RNAseq_data)
				ans2 <- as.integer(readline(prompt=cat("\nwould you like to:\n1) use data from the list\n2) download new RNA-seq data to server\n\n")))
				if(ans2 == 1) #using already existing data
				{
						sra_filename <- readline(prompt=cat("\nThese are the sra files of RNA-seq data available\nwhich would you like to use?\n",sra_files,"\n\n"))
						sraInd_files <- pmatch(sra_filename,sra_files,nomatch=0)
						if(sraInd_files == 0)
						{
							cat("no such file exists.\nplease try again.\n\n")
						}
						else
						{
							rna_ind <- which(RNAseq_data$Name_of_File %in% sra_filename)
							DandT1 <- RNAseq_data$Download_Date[rna_ind]
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
			cat("\nthe folder is empty. It seems that there are no sra files to choose from.\nIn order to proceed you must download RNA-seq data in sra format\n\n")
			ans3 <- readline(prompt=cat("\nwould you like to download new RNA-seq data?\ny/n\n\n"))
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
			cat("\nRNA-seq details:\n\n")
			#name of organism
			org <- readline(prompt=cat("enter the name of the RNA-seqs organism:\n"))
			#name of tissue
			tissue <- readline(prompt=cat("enter the name of the RNA-seqs tissue:\n"))
			#age
			age <- readline(prompt=cat("enter the age of the cells when RNA-seq was extracted:\n"))
			#additional notes
			ans4 <- readline(prompt=cat("\nwould you like to add any additional information?\ny/n\n\n"))
			if(ans4 == "y")
			{
				note <- readline(prompt=cat("write additional information about the RNA-seq data:\n"))
			}
			else
			{
				note <- NA
			}
			
			#getting the link to the file
			link_sra <- readline(prompt=cat("\nenter the link to the sra data file:\n\n"))
			
			#checking that the 'temp' folder is empty, if not then empty it out
			ls_temp <- system("ls ~/Analyze4C/RNAseq/temp/",intern=TRUE)
			#removing sra files from temp
			if(identical(grep(".sra",ls_temp),integer(0)) == "FALSE")
			{
				system("rm ~/Analyze4C/RNAseq/temp/*.sra")
			}
			#removing fastq files from temp
			if(identical(grep(".fastq",ls_temp),integer(0)) == "FALSE")
			{
				system("rm ~/Analyze4C/RNAseq/temp/*.fastq")
			}
			
			#getting the sra file into 'temp'
			system(paste("wget -P ~/Analyze4C/RNAseq/temp",link_sra))
			
			#getting the date and time of downloading RNA-seq data
			DandT1 <- toString(Sys.time())
			DandT2 <- gsub(" ","_",DandT1)
			DandT2 <- gsub(":","",DandT2)

			#getting the name of the sra file
			sra_filename <- system("ls ~/Analyze4C/RNAseq/temp | grep '.sra$'",intern=TRUE)
			
			#moving the sra file from 'temp' to the 'sra' folder
			system(paste("mv ~/Analyze4C/RNAseq/temp/",sra_filename," ~/Analyze4C/RNAseq/sra",sep=""))
			
			#adding all the details from above to the file RNAseq_data.txt
			RNAseq_data[nrow(RNAseq_data)+1,] <- c(sra_filename,org,tissue,age,link_sra,note,DandT1)
			#sorting the list of sra files by name (and sorting the row indices)
			RNAseq_data <- RNAseq_data[order(RNAseq_data$Name_of_File),]
			rownames(RNAseq_data) <- seq(length=nrow(RNAseq_data))
			#adding the new data to the file (by erasing the old one and creating a new one)
			system("rm RNAseq_data.txt")
			write.table(RNAseq_data,"RNAseq_data.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
			
			ans5 <- readline(prompt=cat("\nwould you like to continue with this current file that was just downloaded?\ny/n\n\n"))
			if(ans5 == "n")
			{
				stop
			}
		}
		
		#fastq dump the sra data, send the output to a temp folder
		cat("\n\nfastq dumping:\nretrieve fastq data from sra file\nplease wait...\n\n")
		system(paste("fastq-dump -O  ~/Analyze4C/RNAseq/temp/ ~/Analyze4C/RNAseq/sra/",sra_filename,sep=""))
		
		#get the links to the genome files folder (the folder should be in TAIR10, it has a bunch of files with 'genome' and extensions of 'ebwt' and 'fa') and the genes annotation folder and file (genes.gtf)
		link_genome <- readline(prompt=cat("\nplease provide a link to the whole genomes bowties index file and add '/genome' at the end, this should be in the bowtie index folder of the genomes sequence folder\ne.g.: Arabidopsis/Arabidopsis_thaliana_TAIR10/Sequence/BowtieIndex/genome\n\n"))
		link_genes <- readline(prompt=cat("\nplease provide a link to the annotation of genes folder and add '/genes.gtf' at the end, this should be in the annotation archives folders of the genome\ne.g.: /Arabidopsis/Arabidopsis_thaliana_TAIR10/Annotation/Archives/archive-2013-03-06-09-54-25/Genes/genes.gtf\n\n"))
		
		#getting the name of the fastq file in temp
		fastq_filename <- system("ls ~/Analyze4C/RNAseq/temp | grep '.fastq$'",intern=TRUE)
		
		#tophat - entering the genome files, and the name of the fastq file created from the fastq dump
		#if user wants to change the default parameters, then they can do so
		cat("\nThe default of the tophat running is as follows:\n\n")
		cat("tophat -p 8 -G ",link_genes," -o ~/Analyze4C/RNAseq/temp/ ",link_genome," ~/Analyze4C/RNAseq/temp/",fastq_filename,sep="")
		cat("\n")
		ans6 <- readline(prompt=cat("\nwould you like to change the parameters?\ny/n\n\n"))
		if(ans6 == "y")
		{
			ans7 <- 0
			cat("\nnote: only the parameters are permitted to be changed or added to, not the files or output destinations\n")
			while(ans7 != "y")
			{
				new_parameters <- readline(prompt=cat("\nenter the new parameters freely. Make sure that they are legal:\n\n"))
				cat("\nthis is the tophat command as requested:\n\n")
				cat("tophat ",new_parameters," -G ",link_genes," -o ~/Analyze4C/RNAseq/temp/ ",link_genome," ~/Analyze4C/RNAseq/temp/",fastq_filename,sep="")
				cat("\n")
				tophat_command <- paste("tophat ",new_parameters," -G ",link_genes," -o ~/Analyze4C/RNAseq/temp/ ",link_genome," ~/Analyze4C/RNAseq/temp/",fastq_filename,sep="") #recording the full command used in order to put in tophat_cufflinks_outputs later
				ans7 <- readline(prompt=cat("\n\nwould you like to execute this command?\ny/n\n\n"))
			}
			system(paste("tophat ",new_parameters," -G ",link_genes," -o ~/Analyze4C/RNAseq/temp/ ",link_genome," ~/Analyze4C/RNAseq/temp/",fastq_filename,sep=""))
		}
		else
		{
			tophat_command <- paste("tophat -p 8 -G ",link_genes," -o ~/Analyze4C/RNAseq/temp/ ",link_genome," ~/Analyze4C/RNAseq/temp/",fastq_filename,sep="") #recording the full command used in order to put in tophat_cufflinks_outputs later
			system(paste("tophat -p 8 -G ",link_genes," -o ~/Analyze4C/RNAseq/temp/ ",link_genome," ~/Analyze4C/RNAseq/temp/",fastq_filename,sep=""))
		}
		
		#cufflinks - entering the genome files, and the accepted_hits.bam file from the previous output
		cufflinks_command <- paste("cufflinks -p 8 -G ",link_genes," -o ~/Analyze4C/RNAseq/temp/ ~/Analyze4C/RNAseq/temp/accepted_hits.bam",sep="") #recording the full command used in order to put in tophat_cufflinks_outputs later
		system(paste("cufflinks -p 8 -G ",link_genes," -o ~/Analyze4C/RNAseq/temp/ ~/Analyze4C/RNAseq/temp/accepted_hits.bam",sep=""))
		
		#getting genes.fpkm_tracking
		genes.fpkm <- read.delim("~/Analyze4C/RNAseq/temp/genes.fpkm_tracking",header=TRUE,stringsAsFactors=FALSE)
		
		#save the genes.fpkm_tracking file in a folder, add the date, time, and name of sra to it, add the details of running tophat and cufflinks into a file that records the name of the fastq file the date and these details
		cat("\nsaving and renaming the genes.fpkm_tracking file into the 'FPKM' folder...\n\n")
		system("mv ~/Analyze4C/RNAseq/temp/genes.fpkm_tracking ~/Analyze4C/RNAseq/FPKM")
		fpkm_filename_temp <- strsplit(sra_filename,"[.]")
		fpkm_filename <- paste(fpkm_filename_temp[[1]][1],"_",DandT2,"_genes.fpkm_tracking",sep="")
		system(paste("mv ~/Analyze4C/RNAseq/FPKM/genes.fpkm_tracking ~/Analyze4C/RNAseq/FPKM/",fpkm_filename,sep=""))
		
		cat("\nadding the tophat and cufflinks commands used into tophat_cufflinks_outputs.txt...\n\n")
		tophat_cufflinks_outputs[nrow(tophat_cufflinks_outputs)+1,] <- c(fpkm_filename_temp[[1]][1],DandT1,tophat_command,cufflinks_command)
		tophat_cufflinks_outputs <- tophat_cufflinks_outputs[order(tophat_cufflinks_outputs$Name),]
		rownames(tophat_cufflinks_outputs) <- seq(length=nrow(tophat_cufflinks_outputs))
		system("rm tophat_cufflinks_outputs.txt")
		write.table(tophat_cufflinks_outputs,"tophat_cufflinks_outputs.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
		
		#extracting the FPKMs and their corresponding indices from genes.fpkm_tracking
		cat("\ncreating an FPKM bed file and placing it in the 'FPKM' folder...\n\n")
		inds <- genes.fpkm$locus
		fpkms <- genes.fpkm$FPKM
		len1 <- length(inds)
		all.fpkms <- matrix(0,len1,4)
		counter1 <- 1
		counter2 <- 0
		for(j in 1:len1)
		{
			spl1 <- strsplit(toString(inds[j]),":")
			if(spl1[[1]][1] != "Pt" & spl1[[1]][1] != "Mt") #not taking the lines that include 'Mt' or 'Pt'
			{
				spl2 <- strsplit(spl1[[1]][2],"-")
				all.fpkms[counter1,] <- cbind(as.integer(spl1[[1]][1]),as.integer(spl2[[1]][1]),as.integer(spl2[[1]][2]),fpkms[j])
				counter1 <- counter1 + 1
			}
			else
			{
				counter2 <- counter2 + 1
			}
		}
		all.fpkms <- all.fpkms[1:(len1-counter2),]
		all.fpkms <- all.fpkms[order(all.fpkms[,1],all.fpkms[,2]),]

		#creating a bed file of the FPKMs and their corresponding indices
		fpkmBed_name <- paste(fpkm_filename_temp[[1]][1],"_",DandT2,"_FPKM.bed",sep="")
		write.table(all.fpkms,paste("~/Analyze4C/RNAseq/FPKM/",fpkmBed_name,sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
		
		#erasing all the outputs from temp (user could download the relevant files for themselves)
		cat("\nall the output files (except for the FPKM files placed in the 'FPKM' folder) will be erased")
		cat("\nif you will like to save any for yourself, download the files to your personal computer\n")
		lin <- readline(prompt=cat("\nonce you are done downloading your files -\n\npress any key to continue and erase the files...\n\n"))
		system("rm -rf ~/Analyze4C/RNAseq/temp/*")
	}
	else
	{
		cat("\ninstall the sra toolkit, tophat, and cufflink packages in order to use the function\n\n")
	}
}
