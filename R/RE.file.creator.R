#' @export

#this function first gets an input of the RE site name and organism from, it then checks to see if the sgr file txt file of the RE sites for that specific organism and RE site exists,
#if it doesn't it will create one by asking for the links to the genome annotations and fasta files and getting the sequence of the RE site.
#then the user will decide if to remove the 'Mt' and 'Pt' rows.
#the next part is to ask the user if they would like to remove the unaligned RE sites. The user must input the number of bps that the function will check on both sides of each
#first RE site index, this will then be attempted to be aligned, the idea is that if it doesn't align then these areas around the RE site cannot be aligned when the primer is roughly
#around the size input by the user before. usually the reason for unalignment is because of repeats in the sequence.
#eventually a txt file with the indices of only the aligned RE sites will be created. the format of the name is ["RE_name"."org_name".number of bp checked either sides of Re site bpUArm.txt]
#e.g. 'Dpn.arabidopsis.20bpUArm.txt'
#bp - base pairs. UA - unaligned. rm - removed
#the format of the file should be: +	1	310	GATC

#the last part asks if to remove all the 'blind' RE sites.
#these sites are those that are NOT situated adjacent to the primary RE site used in the 3C/4C process (in our case HindIII)
#we look at each chromosome separately and for every primary RE site, we mark the secondary RE sites (Dpn in our case) that are adjacent to them
#these secondary RE sites are the ones that are saved, all the others will be considered 'blind' and cannot be used
#the name of the file created will contain "_noBlindREs", this of course tells us that the 'blind' Dpn sites are removed

RE.file.creator <- function(RE_name,RE_org)
{
	#checking to see if the txt file for the specific RE site in an organism already exists, if no then it will create it
	name_txt <- paste(RE_name,".",RE_org,".txt",sep="")
	ls_txt <- system("ls ~/Analyze4C/sgr_orgs",intern=TRUE)
	mch <- pmatch(name_txt,ls_txt,nomatch=0)
	if(mch == 0)
	{
		cat("we first need to create an txt file that will later create sgr files for",RE_name,"in",RE_org,"\n\n")	
		
		#getting the links to the genome annotations and genome.fa files
		RE_seq <- readline(prompt="\nplease provide the RE sites sequence:\n\n")
		link1 <- readline(prompt="\nplease provide a link to the whole genomes bowties index file and add '/genome' at the end, this should be in the bowtie index folder of the genomes sequence folder\ne.g.: Arabidopsis/Arabidopsis_thaliana_TAIR10/Sequence/BowtieIndex/genome\n\n")
		link2 <- readline(prompt="\nplease provide a link to the whole genomes fasta file, which should be in the whole genome fasta folder\ne.g.: genome.fa\n\n")
		file.name <- paste(RE_name,".",RE_org,".txt",sep="")
		
		#creates a file with the positions of the RE site (the file will also contain the '+' sign and the sequence of the RE site)
		system(paste("bowtie -a -v 0 ",link1," --suppress 1,6,7 -c ",RE_seq," | awk '$1==\"+\"' | sort -k2,2 -k3,3n > ~/Analyze4C/temp/",file.name,sep=""))
		
		#remove 'Mt' and 'Pt' if user wants
		ans1 <- readline(prompt="\nwould you like to remove 'Mt' and 'Pt'?\ny/n\n")
		if(ans1=="y")
		{
			temp1 <- read.table(paste("~/Analyze4C/temp/",file.name,sep=""))
			system(paste("rm ~/Analyze4C/temp/",file.name,sep=""))
			write.table(temp1[temp1[,2]!='Mt' & temp1[,2]!='Pt',],paste("~/Analyze4C/temp/",file.name,sep=""),sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
		}
		
		#ask if to remove the unaligned RE sites, and how many bp to check from each side.
		#recommend the number of bps for each side should be close to the length of the primer.
		ans2 <- readline(prompt="\nwould you like to remove the RE sites that do not align to the genome?\ny/n\n")
		if(ans2=="y")
		{		
			name2 <- paste(RE_org,"_chromosome_sizes.txt",sep="")
			ls_files <- system("ls ~/Analyze4C/genomes",intern=TRUE)
			match <- pmatch(name2,ls_files,nomatch=0)
			if(match == 0)
			{
				cat("we first need to create a file with the",RE_org,"chromosome sizes\n\n")
				chrom_sizes_fileCreator(1,RE_org,name2)
			}
			
			temp2 <- read.table(paste("~/Analyze4C/temp/",file.name,sep=""))
			inds <- cbind(temp2[,2],temp2[,3])
			#getting the cutter number to add to the indices
			cutter <- readline(prompt="\nwhat number cutter is the RE site?\n")
			cutter <- as.integer(cutter)
			#creating the bed file for the RE site index file
			bed <- cbind(inds[,1],inds[,2],inds[,2] + cutter)
			bed.name <- paste(RE_name,".",RE_org,".bed",sep="")
			write.table(bed,paste("~/Analyze4C/temp/",bed.name,sep=""),sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
			
			#adding a number of bps before and after RE sites
			ans3 <- readline(prompt="\nhow many bps should be checked?\nthe number should be close to the size of the primer\nnote that the check is done both forward and backwards\n")
			ans3 <- as.integer(ans3)

			#checking to see if the txt file for the specific RE site in an organism already exists, if no then it will create it
			removed_name <- paste(RE_name,".",RE_org,".",ans3,"bpUArm.txt",sep="")
			ls_txt2 <- system("ls ~/Analyze4C/sgr_orgs",intern=TRUE)
			mch2 <- pmatch(removed_name,ls_txt2,nomatch=0)
			if(mch2 == 1)
			{
				paste("this file with this number of bps around this RE site already exists\n")
			}
			else
			{
				system(paste("perl ~/Analyze4C/proxy/add_inds.pl ~/Analyze4C/temp/",bed.name," ",ans3," > ~/Analyze4C/temp/temp_bed.bed",sep=""))
				
				#here we remove any indices that exceed the size of the genome
				rd1 <- read.table("~/Analyze4C/temp/temp_bed.bed")
				cond <- 0
				name_temp <- paste("~/Analyze4C/genomes/",name2,sep="")
				sizes <- read.table(name_temp)
				new_all <- c()
				for(j in unique(rd1[,1]))
				{
					tm <- rd1[rd1[,1]==j,]
					if(tm[2,2] < 1)
					{
						tm <- tm[-2,]
						cond <- 1
					}
					
					
					siz <- nrow(tm)
					if(tm[siz-1,3] > sizes[j,2])
					{
						tm <- tm[-(siz-1),]
						cond <- 1
					}
					new_all <- rbind(new_all,tm)
				}
				if(cond == 1)
				{
					write.table(new_all,"~/Analyze4C/temp/temp2_bed.bed",sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
				}
				
				#getting the fasta sequences for the indices in the temp2_bed.bed file
				system(paste("bedtools getfasta -fi",link2,"-bed ~/Analyze4C/temp/temp2_bed.bed -fo ~/Analyze4C/temp/temp_fa.fa"))
				
				#trying to align the sequences from "temp_fa.fa". The important file that will be used is the one with the unaligned sequences
				system(paste("bowtie -m 1 -S -f",link1,"~/Analyze4C/temp/temp_fa.fa --un ~/Analyze4C/temp/notAligned.fa ~/Analyze4C/temp/aligned.sam"))
				
				#editing the 'notAligned' fasta file and turning it into a bed file
				system("perl proxy/edit_notAligned.pl ~/Analyze4C/temp/notAligned.fa > ~/Analyze4C/temp/notAligned.bed")
				
				#comparing the notAligned sites with all the sites, those that don't align will get a 1 and those that do a 0
				system("bedtools intersect -a ~/Analyze4C/temp/temp_bed.bed -b ~/Analyze4C/temp/notAligned.bed -f 1.0 -c > ~/Analyze4C/temp/whichDidntAlign.bed")
				
				#getting the indices of the REs that aligned and those that didn't, 1 represents sites that didn't align and 0 ones that did
				#note: this part is instead of the function 'bad_sites.R'
				rd2 <- read.table("~/Analyze4C/temp/whichDidntAlign.bed")
				len <- nrow(temp2)
				fin <- rep(0,len)
				for (x in 1:len)
				{
					y <- x*2
					if(rd2[y-1,4]>=1 || rd2[y,4]>=1)
					{
						fin[x] <- 1
					}
				}

				#here i might want to have an option to keep the data of RE sites that don't align
				
				#creating an txt file of only the RE sites that align, this is the file that contains the final RE sites that will be used when creating sgr files
				#the name of the file will be ["RE_name"."org_name".number of bp checked either sides of Re site bpUArm.txt]
				#bp - base pairs. UA - unaligned. rm - removed
				#the format of the file should be: +	1	310	GATC
				aligned_sgr <- temp2[which(fin %in% 0),]
				write.table(aligned_sgr,paste("~/Analyze4C/sgr_orgs/",removed_name,sep=""),sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
				
				#erasing all the files that there is no use for anymore
				system("rm ~/Analyze4C/temp/temp_bed.bed")
				system("rm ~/Analyze4C/temp/temp_fa.fa")
				system("rm ~/Analyze4C/temp/notAligned.fa")
				system("rm ~/Analyze4C/temp/notAligned.bed")
				system("rm ~/Analyze4C/temp/aligned.sam")
				system("rm ~/Analyze4C/temp/temp2_bed.bed")
				system(paste("rm ~/Analyze4C/temp/",bed.name,sep=""))
				system("rm ~/Analyze4C/temp/whichDidntAlign.bed")
			}
		}
		
		ans4 <- readline(prompt="\nwould you like to remove the 'blind' RE sites?\ny/n\n")
		if(ans4=="y")
		{
			#getting the file from which we will remove the blind RE sites
			if(ans2 == "y")
			{
				temp3 <- read.table(paste("~/Analyze4C/sgr_orgs/",removed_name,sep=""))
				file.name2 <- paste("~/Analyze4C/sgr_orgs/",removed_name,sep="")
			}
			else
			{
				temp3 <- read.table(paste("~/Analyze4C/temp/",file.name,sep=""))
				#temp3 <- read.delim(paste("~/Analyze4C/temp/",file.name,sep=""),header=FALSE,quote="")
				file.name2 <- paste("~/Analyze4C/sgr_orgs/",file.name,sep="")
			}
			
			RE_seq2 <- readline(prompt="\nplease provide the sequence of the first RE site cut in 4C/3C (e.g. HindIII sequence is AAGCTT):\n\n")
							
			#align the sequence, get the indices where the RE site sits
			system(paste("bowtie -a -v 0 ",link1," --suppress 1,6,7 -c ",RE_seq2," | awk '$1==\"+\"' | sort -k2,2 -k3,3n > ~/Analyze4C/temp/firstRE_temp.txt",sep=""))

			#remove 'Mt' and 'Pt' if user wants
			ans1 <- readline(prompt="\nwould you like to remove 'Mt' and 'Pt'?\ny/n\n")
			if(ans1=="y")
			{
				temp4 <- read.table("~/Analyze4C/temp/firstRE_temp.txt")
				system("rm ~/Analyze4C/temp/firstRE_temp.txt")
				write.table(temp4[temp4[,2]!='Mt' & temp4[,2]!='Pt',],"~/Analyze4C/temp/firstRE_temp.txt",sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
			}
			
			#getting the firstRE_temp data
			firstRE_temp <- read.table(paste("~/Analyze4C/temp/firstRE_temp.txt",sep=""))
			#firstRE_temp <- read.delim(paste("~/Analyze4C/temp/firstRE_temp.txt",sep=""),header=FALSE,quote="")
			
			#adding a fourth column with the same number for every row (here we are adding 1)
			firstRE_temp2 <- unname(cbind(firstRE_temp,rep(1,nrow(firstRE_temp))))
			
			#taking the previous made Dpn file, creating a new file with that data by adding the same number in each row (here we are adding 2)
			secondRE_temp <- unname(cbind(temp3,rep(2,nrow(temp3))))
			
			#combining both files to one, sorting the indices (first by chromosome then by index)
			all_firstANDsecondRE <- rbind(data.frame(firstRE_temp2),data.frame(secondRE_temp))
			all_firstANDsecondRE <- all_firstANDsecondRE[order(all_firstANDsecondRE[,2],all_firstANDsecondRE[,3]),]
			rownames(all_firstANDsecondRE) <- seq(length=nrow(all_firstANDsecondRE))
			
			#go over the fourth column, everywhere that we see more than two '2' in a row (without a '1' in between), we will mark all the middle dpn sites (e.g. by changing their number to '3')			
			all2 <- c()
			for(z in unique(all_firstANDsecondRE[,2]))
			{
				all_firstANDsecondRE_chr <- all_firstANDsecondRE[all_firstANDsecondRE[,2]==z,]
				rownames(all_firstANDsecondRE_chr) <- seq(length=nrow(all_firstANDsecondRE_chr))
				first_inds <- which(all_firstANDsecondRE_chr[,5] %in% 1)
				for(a in first_inds)
				{
					if(a==1 & length(first_inds)>=2)
					{
						if(all_firstANDsecondRE_chr[a+1,5]==2)#if the next index points to a secondary RE site (Dpn)
						{
							all_firstANDsecondRE_chr[a+1,5] <- 3
						}
					}
					else if(a==nrow(all_firstANDsecondRE_chr) & length(first_inds)>=2)
					{
						if(all_firstANDsecondRE_chr[a-1,5]==2)#if the previous index points to a secondary RE site (Dpn)
						{
							all_firstANDsecondRE_chr[a-1,5] <- 3
						}						
					}
					else if(length(first_inds)<2 & length(first_inds)>0)
					{
						cat("\nthis chromosome contains only 1 RE site\n\n")
					}
					else
					{
						if(all_firstANDsecondRE_chr[a+1,5]==2)#if the next index points to a secondary RE site (Dpn)
						{
							all_firstANDsecondRE_chr[a+1,5] <- 3
						}

						if(all_firstANDsecondRE_chr[a-1,5]==2)#if the previous index points to a secondary RE site (Dpn)
						{
							all_firstANDsecondRE_chr[a-1,5] <- 3
						}						
					}
				}
				all2 <- rbind(all2,all_firstANDsecondRE_chr)
			}	
			
			rownames(all2) <- seq(length=nrow(all2)) #giving row names in order
			
			#we only take the rows with '3' in the fourth column, meaning we are left with all the dpn sites that can potentially be used in the 4C
			final <- all2[all2[,5]==3,]
			final <- final[,-5]
			
			#creating the file with the final data of RE sites
			spl <- strsplit(file.name2,".txt")
			write.table(final,paste(spl[1],".noBlindREs.txt",sep=""),sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
			
			system("rm ~/Analyze4C/temp/firstRE_temp.txt")
		}
		
		#moving the original txt file of RE sites to the folder 'sgr_orgs'
		system(paste("mv ~/Analyze4C/temp/",file.name," ~/Analyze4C/sgr_orgs",sep=""))
	}
	else
	{
		cat("\nthe RE sites file for these specific sites and organism already exists\n\n")
	}
}
