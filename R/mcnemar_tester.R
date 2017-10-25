#' @export

#this function performs the mcenmar test on three bed files, where we intersect two of them with the third and take the sum of intersections for each and compare them
#the idea of the test is to see how significantly different two binary paired groups are
#translating this to our case, since we are intersecting two files to the same query and since it is binary for each bp (intersected or not), it fits the category

mcnemar_tester <- function(query=0,fileA=0,fileB=0)
{
	if(query==0)
	{
		#choosing a query file (a contact band file). the file name and location is chosen
		repeat
		{
			ans1 <- readline(prompt=cat("\nwhat type of file would you like to choose for the query file:\n1) contact bands\n2) FPKM bed file\n\n"))
			if(ans1 == 1)
			{
				ls_conts <- system("ls ~/Analyze4C/contact_bands/",intern=TRUE)
				repeat
				{
					if(length(ls_conts) != 0)
					{	
						conts_filename <- readline(prompt=cat("\nchoose the contact bands file that you would like to use:\n",ls_conts,"",sep="\n"))
						ind_files1 <- pmatch(conts_filename,ls_conts,nomatch=0)
						if(ind_files1 == 0)
						{
							cat("no such file exists.\nplease try again.\n\n")
						}
						else
						{
							query <- paste("~/Analyze4C/contact_bands/",conts_filename,sep="")
							break
						}
					}
					else
					{
						cat("\nno contact bands files exist\n")
						cat("\nchoose a different type of file to use for the query\n\n")
					}
				}	
			}
			else if(ans1 == 2)
			{
				ls_FPKM <- system("ls ~/Analyze4C/RNAseq/FPKM | grep '.bed$'",intern=TRUE)
				repeat
				{
					if(length(ls_FPKM) != 0)
					{	
						FPKM_filename1 <- readline(prompt=cat("\nchoose the FPKM file that you would like to use:\n",ls_FPKM,"",sep="\n"))
						ind_files2 <- pmatch(FPKM_filename1,ls_FPKM,nomatch=0)
						if(ind_files2 == 0)
						{
							cat("no such file exists.\nplease try again.\n\n")
						}
						else
						{
							query <- paste("~/Analyze4C/RNAseq/FPKM/",FPKM_filename1,sep="")
							break
						}
					}
					else
					{
						cat("\nno FPKM files exist\n")
						cat("\nchoose a different type of file to use for the query\n\n")
					}					
				}				
			}
			
			if(query != 0)
			{
				break
			}
		}
	}
	
	if(fileA==0)
	{
		#choosing an intesection file (FPKM file). the file name and location is chosen
		repeat
		{
			ans1 <- readline(prompt=cat("\nwhat type of file would you like to choose for fileA:\n1) contact bands\n2) FPKM bed file\n\n"))
			if(ans1 == 1)
			{
				ls_conts <- system("ls ~/Analyze4C/contact_bands/",intern=TRUE)
				repeat
				{
					if(length(ls_conts) != 0)
					{	
						conts_filename <- readline(prompt=cat("\nchoose the contact bands file that you would like to use:\n",ls_conts,"",sep="\n"))
						ind_files1 <- pmatch(conts_filename,ls_conts,nomatch=0)
						if(ind_files1 == 0)
						{
							cat("no such file exists.\nplease try again.\n\n")
						}
						else
						{
							fileA <- paste("~/Analyze4C/contact_bands/",conts_filename,sep="")
							break
						}
					}
					else
					{
						cat("\nno contact bands files exist\n")
						cat("\nchoose a different type of file to use for fileA\n\n")
					}
				}	
			}
			else if(ans1 == 2)
			{
				ls_FPKM <- system("ls ~/Analyze4C/RNAseq/FPKM | grep '.bed$'",intern=TRUE)
				repeat
				{
					if(length(ls_FPKM) != 0)
					{	
						FPKM_filename1 <- readline(prompt=cat("\nchoose the FPKM file that you would like to use:\n",ls_FPKM,"",sep="\n"))
						ind_files2 <- pmatch(FPKM_filename1,ls_FPKM,nomatch=0)
						if(ind_files2 == 0)
						{
							cat("no such file exists.\nplease try again.\n\n")
						}
						else
						{
							fileA <- paste("~/Analyze4C/RNAseq/FPKM/",FPKM_filename1,sep="")
							break
						}
					}
					else
					{
						cat("\nno FPKM files exist\n")
						cat("\nchoose a different type of file to use for fileA\n\n")
					}					
				}				
			}
			
			if(fileA != 0)
			{
				break
			}
		}		
	}
	
	if(fileB==0)
	{
		#choosing an intesection file (FPKM file). the file name and location is chosen
		repeat
		{
			ans1 <- readline(prompt=cat("\nwhat type of file would you like to choose for fileB:\n1) contact bands\n2) FPKM bed file\n\n"))
			if(ans1 == 1)
			{
				ls_conts <- system("ls ~/Analyze4C/contact_bands/",intern=TRUE)
				repeat
				{
					if(length(ls_conts) != 0)
					{	
						conts_filename <- readline(prompt=cat("\nchoose the contact bands file that you would like to use:\n",ls_conts,"",sep="\n"))
						ind_files1 <- pmatch(conts_filename,ls_conts,nomatch=0)
						if(ind_files1 == 0)
						{
							cat("no such file exists.\nplease try again.\n\n")
						}
						else
						{
							fileB <- paste("~/Analyze4C/contact_bands/",conts_filename,sep="")
							break
						}
					}
					else
					{
						cat("\nno contact bands files exist\n")
						cat("\nchoose a different type of file to use for fileB\n\n")
					}
				}	
			}
			else if(ans1 == 2)
			{
				ls_FPKM <- system("ls ~/Analyze4C/RNAseq/FPKM | grep '.bed$'",intern=TRUE)
				repeat
				{
					if(length(ls_FPKM) != 0)
					{	
						FPKM_filename1 <- readline(prompt=cat("\nchoose the FPKM file that you would like to use:\n",ls_FPKM,"",sep="\n"))
						ind_files2 <- pmatch(FPKM_filename1,ls_FPKM,nomatch=0)
						if(ind_files2 == 0)
						{
							cat("no such file exists.\nplease try again.\n\n")
						}
						else
						{
							fileB <- paste("~/Analyze4C/RNAseq/FPKM/",FPKM_filename1,sep="")
							break
						}
					}
					else
					{
						cat("\nno FPKM files exist\n")
						cat("\nchoose a different type of file to use for fileB\n\n")
					}					
				}				
			}
			
			if(fileB != 0)
			{
				break
			}
		}
	}

	#merging the fragments the overlap, this prevents multiple intersections of the same sections
	system(paste("bedtools merge -i",query,"> ~/Analyze4C/temp/query_merged.bed"))
	system(paste("bedtools merge -i",fileA,"> ~/Analyze4C/temp/fileA_merged.bed"))
	system(paste("bedtools merge -i",fileB,"> ~/Analyze4C/temp/fileB_merged.bed"))
	
	#intersections:	
#	system("bedtools intersect -a ~/Analyze4C/temp/query_merged.bed -b ~/Analyze4C/temp/query_merged.bed -wo > ~/Analyze4C/temp/int_query.txt")
#	int_query <- read.table("~/Analyze4C/temp/int_query.txt",header=FALSE,stringsAsFactors=FALSE)
#	sum_int_query <- sum(int_query[,7])
#	system("bedtools intersect -a ~/Analyze4C/temp/query_merged.bed -b ~/Analyze4C/temp/fileA_merged.bed -wo > ~/Analyze4C/temp/intA.txt")
#	intA <- read.table("~/Analyze4C/temp/intA.txt",header=FALSE,stringsAsFactors=FALSE)
#	sum_intA <- sum(intA[,7])
#	system("bedtools intersect -a ~/Analyze4C/temp/query_merged.bed -b ~/Analyze4C/temp/fileB_merged.bed -wo > ~/Analyze4C/temp/intB.txt")
#	intB <- read.table("~/Analyze4C/temp/intB.txt",header=FALSE,stringsAsFactors=FALSE)
#	sum_intB <- sum(intB[,7])	
#	system("bedtools intersect -a ~/Analyze4C/temp/fileA_merged.bed -b ~/Analyze4C/temp/fileB_merged.bed > ~/Analyze4C/temp/intAB.bed")
#	intAB <- read.table("~/Analyze4C/temp/intAB.bed",header=FALSE,stringsAsFactors=FALSE)
#	system("bedtools intersect -a ~/Analyze4C/temp/query_merged.bed -b ~/Analyze4C/temp/intAB.bed -wo > ~/Analyze4C/temp/intABQ.txt")
#	intABQ <- read.table("~/Analyze4C/temp/intABQ.txt",header=TRUE,stringsAsFactors=FALSE)
#	sum_intABQ <- sum(intABQ[,7])
	
	#intersections:
	system("bedtools intersect -a ~/Analyze4C/temp/query_merged.bed -b ~/Analyze4C/temp/query_merged.bed -wo > ~/Analyze4C/temp/int_query.txt")	
	if(file.info("~/Analyze4C/temp/int_query.txt")$size != 0)
	{
		int_query <- read.table("~/Analyze4C/temp/int_query.txt",header=FALSE,stringsAsFactors=FALSE)
	}
	else
	{
		int_query <- data.frame(0,0,0,0,0,0,0,0)
	}
	sum_int_query <- sum(int_query[,7])

	system("bedtools intersect -a ~/Analyze4C/temp/query_merged.bed -b ~/Analyze4C/temp/fileA_merged.bed -wo > ~/Analyze4C/temp/intA.txt")
	if(file.info("~/Analyze4C/temp/intA.txt")$size != 0)
	{
		intA <- read.table("~/Analyze4C/temp/intA.txt",header=FALSE,stringsAsFactors=FALSE)
	}
	else
	{
		intA <- data.frame(0,0,0,0,0,0,0,0)
	}
	sum_intA <- sum(intA[,7])

	system("bedtools intersect -a ~/Analyze4C/temp/query_merged.bed -b ~/Analyze4C/temp/fileB_merged.bed -wo > ~/Analyze4C/temp/intB.txt")	
	if(file.info("~/Analyze4C/temp/intB.txt")$size != 0)
	{
		intB <- read.table("~/Analyze4C/temp/intB.txt",header=FALSE,stringsAsFactors=FALSE)
	}
	else
	{
		intB <- data.frame(0,0,0,0,0,0,0,0)
	}
	sum_intB <- sum(intB[,7])	

	system("bedtools intersect -a ~/Analyze4C/temp/fileA_merged.bed -b ~/Analyze4C/temp/fileB_merged.bed > ~/Analyze4C/temp/intAB.bed")
	if(file.info("~/Analyze4C/temp/intAB.bed")$size != 0)
	{
		intAB <- read.table("~/Analyze4C/temp/intAB.bed",header=FALSE,stringsAsFactors=FALSE)
	}
	else
	{
		intAB <- data.frame(0,0,0,0,0,0,0,0)
	}

	system("bedtools intersect -a ~/Analyze4C/temp/query_merged.bed -b ~/Analyze4C/temp/intAB.bed -wo > ~/Analyze4C/temp/intABQ.txt")
	if(file.info("~/Analyze4C/temp/intABQ.txt")$size != 0)
	{
		intABQ <- read.table("~/Analyze4C/temp/intABQ.txt",header=FALSE,stringsAsFactors=FALSE)
	}
	else
	{
		intABQ <- data.frame(0,0,0,0,0,0,0,0)
	}
	sum_intABQ <- sum(intABQ[,7])
		
	#calculating the sum of intersections of only intA, only intB, both, and neither:
	only_intA <- sum_intA - sum_intABQ
	only_intB <- sum_intB - sum_intABQ
	both <- sum_intABQ
	neither <- sum_int_query - (sum_intA + sum_intB - sum_intABQ)
	mat <- matrix(c(both,only_intA,only_intB,neither),2,2)

	#adding names to mat
	colnames(mat) <- c("fileA_intersected","fileA_NOTintersected")
	rownames(mat) <- c("fileB_intersected","fileB_NOTintersected")
	
	#printing mat
	cat("\nthe results of the chosen files intersected:\nquery - ",query,"\nfileA - ",fileA,"\nfileB - ",fileB,"\n\n")
	print(mat)
	
	#performing the mcnemar test
	mc <- mcnemar.test(mat)
	print(mc)
	
	#printing results
	if(only_intA > only_intB)
	{
		cat("intA has more intersections than intB\n")
	}
	else if(only_intA < only_intB)
	{
		cat("intB has more intersections than intA\n")
	}
	else if(only_intA == only_intB)
	{
		cat("the number of intersections of intA and intB are equal\n")
	}
	
	if(mc$p.value <= 0.0001)
	{	
		cat("\nthere is a significant difference, of p < 0.0001, between the two files intersections\n")
	}
	else if(mc$p.value <= 0.001)
	{
		cat("\there is a significant difference, of p < 0.001, between the two files intersections\n")
	}
	else if(mc$p.value <= 0.01)
	{
		cat("\there is a significant difference, of p < 0.01, between the two files intersections\n")
	}
	else if(mc$p.value <= 0.05)
	{
		cat("\there is a significant difference, of p < 0.05, between the two files intersections\n")
	}
	
	#removing all the files from the temp folder
	system("rm ~/Analyze4C/temp/query_merged.bed")
	system("rm ~/Analyze4C/temp/fileA_merged.bed")
	system("rm ~/Analyze4C/temp/fileB_merged.bed")
	system("rm ~/Analyze4C/temp/int_query.txt")
	system("rm ~/Analyze4C/temp/intA.txt")
	system("rm ~/Analyze4C/temp/intB.txt")
	system("rm ~/Analyze4C/temp/intAB.bed")	
	system("rm ~/Analyze4C/temp/intABQ.txt")
}	
