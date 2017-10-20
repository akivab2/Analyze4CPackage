#' @export

#note: maybe i should create a progress bar. i could look at this link for information:
# https://ryouready.wordpress.com/2009/03/16/r-monitor-function-progress-with-a-progress-bar/

main.menu <- function()
{	
	#not sure if to have this whole section here. this is important that ggplot2 is installed. if the user already has the package then all that needs to be done
	#is "library(ggplot2)", but maybe either are not necessary
	#install.packages("ggplot2")
	library(ggplot2)
	
	#####################################################################################################################
	
	#note: i need to add the proxy files folder with its contents
	#if "Analyze4C" is not created then create the folder from the root and create all the internal folder as well
	createRoot()
	
	#####################################################################################################################

	#loading all the functions and files in 'lib'
	
	pathnames <- list.files(pattern="[.]R$", path="~/Analyze4C/lib", full.names=TRUE)
	sapply(pathnames, FUN=source)

	#####################################################################################################################
	
	#note: this whole section could probably be removed
	
	#checking to see if the user is using the function in the correct directory - Analyze4C
	#the code will terminate if not
	cur_dir <- system("pwd | rev | cut -d '/' -f1 | rev",intern=TRUE)
	if(cur_dir != "Analyze4C")
		{stop("you are running the code in a wrong directory.\n please run in a directory called 'Analyze4C'.")}
	
	#if we are starting in the root folder then enter the folder "Analyze4C"
	#setwd("Analyze4C")

	
	#####################################################################################################################
	
	#creating the file "Experiments_4C" if it doesn't exist
	#the file contains the list of experiments which contains the bait, tissue, experiment name or lane
	#the list will be filled in by user while they use the program, every time they enter a new experiment it will add a new row
	Experiments_4C <- createORget_Experiments_4C()
	
	#####################################################################################################################
	
	
	flag0 <- 1
	while(flag0)
	{
		# i need to add options of printing the lists of the files in each of the categories (sgr files, fastq files etc.)
		#need to add an option of printing 'Experiments_4C' just to view the contents
		cat("\nAnalyze-4C menu\nPlease enter number of choice:\n\n")
		switch(menu(c("create a chromosome sizes file for genome","extract reads from sequencing","find experiment in Experiments_4C and add if non existant","align experiment to genome","create RE site sgr file",
		"create raw data sgr file","translocation fixing","barcode contamination fixing","calculate raw data stats","get bp to RE site ratio","create distribution plots and normal distribution testing","calculate P-scores",
		"False discover Rate (FDR)","change coverage","coverage VS f-measure, precision, and recall","RE windows VS precision, recall, f-measure, and intersections (static cutoff)","cutoff VS precision, recall, f-measure, and intersections (static window size)",
		"combine 2 raw data files into one","create contact bands","RNA-seq aligner and FPKM bed files creator","contacts VS. expression (FPKM from RNA-seq)","ChIP-seq aligner and peaks bed files creator","contacts VS. ChIPseq (peaks)","calculate values at quantiles and quantiles of values","randomize raw data","Exit"),graphics=TRUE) + 1,cat("Nothing done\n"),
		#creating chromosome sizes file for genome	
		{
			chrom_sizes_fileCreator()
		},
		#extract reads from sequencing
		{
			#inputs: 'Experiments_4C' which contains the list of experiments
			Experiments_4C <- extract_reads(Experiments_4C)
		},
		{
		#find or add in Experiments_4C
			cat("\nin order to find an experiment in Experiments_4C you must enter following details:\n\n")
			inp_bait <- readline(prompt="enter the name of the bait:\n")
			inp_tissue <- readline(prompt="enter the name of the tissue:\n")
			inp_lane <- readline(prompt="enter the name of the lane:\n")
			out1 <- findIn.Experiments_4C(inp_bait,inp_tissue,inp_lane,Experiments_4C)
			Experiments_4C <- out1[[1]]
		},
		{
		#aligning reads
			#inputs: 'Experiments_4C'
			aligner(Experiments_4C)
			#setTOroot()
		},
		{
		#RE file creator
			#getting the name of RE site and organism
			RE_name <- readline(prompt="\nplease provide the RE sites name:\n\n")
			RE_org <- readline(prompt="\nplease provide the organisms name:\n\n")
			
			RE.file.creator(RE_name,RE_org)
		},
		{
		#creating raw data files (sgr files)
			#inputs: 'Experiments_4C'
			rawData.creator(Experiments_4C)
			Experiments_4C <- createORget_Experiments_4C()
			#setTOroot()
		},
		{
		#translocation fixing
			#first we get the data from the file pScores_of_rearranged.txt
			pScores_of_rearranged <- read.table("pScores_of_rearranged.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			
			#first we get the data from the file rearranged_rawData.txt
			rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
						
			rawData_rearranger(Experiments_4C,pScores_of_rearranged,rearranged_rawData)
		},
		{
		#barcode contamination fixing
		
		},
		{
		#calculating the raw data stats
			#inputs: 'Experiments_4C'
			rawData.stats(Experiments_4C)
			#setTOroot()
		},
		{
		#get bp to RE site ratio
			RE.per.bp()
		},
		{
		#create distribution plots and normal distribution testing
			distribution_plotter(Experiments_4C)
		},
		{
		#calculating p scores
			calculate.pScore(Experiments_4C)
		},
		{
		#calculating FDR
			FDR_main.menu(Experiments_4C)
		},
		{
		#coverage changing
			#first we get the data from the file coverageChanged_Experiments.txt
			coverageChanged_Experiments <- read.table("coverageChanged_Experiments.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			
			coverage.changer.main_menu(Experiments_4C,coverageChanged_Experiments)
		},
		{
		#coverage VS f-measure, precision, and recall
			#first we get the data from the files coverageVSfmeasure_plots.txt and rearranged_rawData.txt
			coverageVSfmeasure_plots <- read.table("coverageVSfmeasure_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			
			ans1 <- as.integer(readline(prompt=cat("\nchoose the method you would like to use:\n\n1) remove coverage from the same file that is compared to\n2) remove coverage from a different file that is compared to\n\n")))
			if(ans1 == 1) #compare the removed coverage with the same original file
			{
				coverage_change_tester(Experiments_4C,coverageVSfmeasure_plots,rearranged_rawData)
			}
			else if(ans1 == 2) #compare the removed coverage with a different file
			{
				coverage_change_tester2(Experiments_4C,coverageVSfmeasure_plots,rearranged_rawData)
			}
		},
		{
		#RE windows VS precision, recall, f-measure, and intersections (static cutoff)
			#first we get the data from the files REwindowIntersections_plots.txt and rearranged_rawData.txt
			REwindowIntersections_plots <- read.table("REwindowIntersections_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			
			RE_window_compare(Experiments_4C,REwindowIntersections_plots,rearranged_rawData)
		},
		{
		#cutoff VS precision, recall, f-measure, and intersections (static window size)
			#first we get the data from the files coIntersections_plots.txt and rearranged_rawData.txt
			coIntersections_plots <- read.table("coIntersections_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			
			CO_compare(Experiments_4C,coIntersections_plots,rearranged_rawData)
		},		
		{
		#combining 2 raw data files to one
			combineLanes(Experiments_4C)
			Experiments_4C <- createORget_Experiments_4C()			
		},
		{
		#create contact bands
			#first we get the data from the file rearranged_rawData.txt
			rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			
			makeContactBands.mainMenu(Experiments_4C,rearranged_rawData)
		},
		{
		#RNA-seq aligner and FPKM bed files creator
			#first we get the data from the file RNAseq_data.txt and tophat_cufflinks_outputs.txt
			RNAseq_data <- read.table("RNAseq_data.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			tophat_cufflinks_outputs <- read.table("tophat_cufflinks_outputs.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			
			rnaSeq_aligner(RNAseq_data,tophat_cufflinks_outputs)
		},
		{
		#contacts VS. expression (FPKM from RNA-seq)
			Experiments_4C <- createORget_Experiments_4C()
		
			choice1 <- as.integer(readline(prompt=cat("\nchoose what type of comparison you would like to do:\n1) compare by quartiles \n2) correlate p-scores and FPKM\n3) correlate reads and FPKM\n4) compare one contacts band dataset with two expression datasets\n5) compare one expression dataset with two contact bands datasets\n\n")))
			if(choice1 == 1)
			{
				#first we get the data from the files rearranged_rawData.txt and expressionVScontacts_sumOFintersections_plots.txt
				rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				expressionVScontacts_sumOFintersections_plots <- read.table("expressionVScontacts_sumOFintersections_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
		
				rnaSeqVSContacts_quantiles(Experiments_4C,expressionVScontacts_sumOFintersections_plots,rearranged_rawData)
			}
			else if(choice1 ==2)
			{
				#first we get the data from the files rearranged_rawData.txt and expressionVScontacts_correlation_plots.txt
				rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				expressionVScontacts_correlation_plots <- read.table("expressionVScontacts_correlation_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				
				rnaSeqCorContacts(Experiments_4C,rearranged_rawData,expressionVScontacts_correlation_plots)
			}
			else if(choice1 ==3)
			{
				#first we get the data from the files rearranged_rawData.txt and expressionVScontacts_correlation_plots.txt
				rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				expressionVScontacts_correlation_plots <- read.table("expressionVScontacts_correlation_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				
				rnaSeqCorContacts_forReads(Experiments_4C,rearranged_rawData,expressionVScontacts_correlation_plots)
			}			
			else if(choice1 == 4)
			{
				#first we get the data from the files rearranged_rawData.txt and expressionVScontacts_plots.txt
				rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				expressionVScontacts_plots <- read.table("expressionVScontacts_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				
				tissue_Expression_comparison(Experiments_4C,expressionVScontacts_plots,rearranged_rawData)
			}
			else if(choice1 ==  5)
			{
				#first we get the data from the files rearranged_rawData.txt and expressionVScontacts_plots.txt
				rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				expressionVScontacts_plots <- read.table("expressionVScontacts_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")

				expression_Tissue_comparison(Experiments_4C,expressionVScontacts_plots,rearranged_rawData)
			}
		},
		{
		#ChIP-seq aligner and peaks bed files creator
			#first we get the data from the file ChIPseq_data.txt and MACS_outputs.txt
			ChIPseq_data <- read.table("ChIPseq_data.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			MACS_outputs <- read.table("MACS_outputs.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			
			ChIPseq_aligner(ChIPseq_data,MACS_outputs)
		},
		{
		#contacts VS. ChipSeq (peaks)
			Experiments_4C <- createORget_Experiments_4C()
		
			choice2 <- as.integer(readline(prompt=cat("\nchoose what type of comparison you would like to do:\n1) compare by quartiles \n2) correlate p-scores and peaks(p-scores or tags)\n3) correlate reads and peaks(p-scores or tags)\n4) compare one contacts band dataset with two ChIPseq datasets\n5) compare one ChIPseq dataset with two contact bands datasets\n\n")))
			if(choice2 == 1)
			{
				#first we get the data from the files rearranged_rawData.txt and ChIPseqVScontacts_sumOFintersections_plots.txt
				rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				ChIPseqVScontacts_sumOFintersections_plots <- read.table("ChIPseqVScontacts_sumOFintersections_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
		
				ChIPSeqVSContacts_quantiles(Experiments_4C,ChIPseqVScontacts_sumOFintersections_plots,rearranged_rawData)
			}
			else if(choice2 ==2)
			{
				#first we get the data from the files rearranged_rawData.txt and ChIPseqVScontacts_correlation_plots.txt
				rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				ChIPseqVScontacts_correlation_plots <- read.table("ChIPseqVScontacts_correlation_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				
				ChIPSeqCorContacts(Experiments_4C,rearranged_rawData,ChIPseqVScontacts_correlation_plots)
			}
			else if(choice2 ==3)
			{
				#first we get the data from the files rearranged_rawData.txt and ChIPseqVScontacts_correlation_plots.txt
				rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				ChIPseqVScontacts_correlation_plots <- read.table("ChIPseqVScontacts_correlation_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				
				ChIPSeqCorContacts_forReads(Experiments_4C,rearranged_rawData,ChIPseqVScontacts_correlation_plots)
			}			
			else if(choice2 == 4)
			{
				#first we get the data from the files rearranged_rawData.txt and ChIPseqVScontacts_plots.txt
				rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				ChIPseqVScontacts_plots <- read.table("ChIPseqVScontacts_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			
				tissue_ChIPseq_comparison(Experiments_4C,ChIPseqVScontacts_plots,rearranged_rawData)
			}
			else if(choice2 ==  5)
			{
				#first we get the data from the files rearranged_rawData.txt and ChIPseqVScontacts_plots.txt
				rearranged_rawData <- read.table("rearranged_rawData.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
				ChIPseqVScontacts_plots <- read.table("ChIPseqVScontacts_plots.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
			
				ChIPseq_Tissue_comparison(Experiments_4C,ChIPseqVScontacts_plots,rearranged_rawData)
			}
		},
		{
		#calculate values at quantiles and quantiles of values
			calculate_quantiles(Experiments_4C)
		},
		{
		#randomize raw data
			Randomize_rawData()	
		},			
		#exit the program
		{flag0 <- 0;}
		)
	}
}
