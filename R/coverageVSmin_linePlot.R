#' @export

#this function creates a graph of all the chromosomes, trans, or by chromosome. Showing with lines how the minimum number of reads influences the coverage percentage.

coverageVSmin_linePlot <- function(dat,file.name,cis)
{
	flag1 <- 1
	while(flag1 == 1)
	{
		#choosing the number of examples to be used in the plot:
		flag2 <- 0
		while(flag2 == 0)
		{
			flag2 <- 1
			num <- as.integer(readline(prompt=cat("\nchoose how many examples of the minimum number of reads or p-scores would you like to view in the plot (the minimum number appear in order from smallest to biggest on plot and are taken in that order from the data):\n\n")))
			for(m in unique(dat[,1]))
			{
				if(num > length(unique(dat[dat[,1]==m,3])))
				{
					flag2 <- 0
					cat("\nthe number chosen is too large\nchoose a different one\n\n")	
					break
				}
			}	
		}
		
		#asking if the plot should include the 100% coverage, meaning that the minimum reads/p-score will be 0
		ans <- readline(prompt=cat("\nshould the plot include 100% coverage?\ny/n\n\n"))
		if(ans == "y")
		{
			fi <- 1
		}
		else
		{
			fi <- 2
		}

		cat("\ncreating the plot\nplease wait...\n\n")
		#all
		rds <- sort(unique(dat[,3]))
		vec <- rep(0,length(rds))
		ind <- 1 #this will be the index counter for vec
		for(x in rds)
		{
			vec[ind] <- ((nrow(dat[dat[,3]>=x,]))/(nrow(dat)))*100 #getting the coverage percentage for all the chromosomes
			ind <- ind + 1
		}
		chromosomes <- rep("all",length(rds))
		df_all <- data.frame(chromosomes,rds,vec)
		df_all <- df_all[fi:num,]

		#trans
		rds <- sort(unique(dat[dat[,1]!=cis,3]))
		vec <- rep(0,length(rds))
		ind <- 1 #this will be the index counter for vec
		for(x in rds)
		{
			vec[ind] <- ((nrow(dat[dat[,1]!=cis & dat[,3]>=x,]))/(nrow(dat[dat[,1]!=cis,])))*100 #getting the coverage percentage for trans
			ind <- ind + 1
		}
		chromosomes <- rep("trans",length(rds))
		df_trans <- data.frame(chromosomes,rds,vec)
		df_trans <- df_trans[fi:num,]
			
		#by chromosome
		df_chrs <- c()
		for(m in unique(dat[,1]))
		{
			rds <- sort(unique(dat[dat[,1]==m,3]))
			vec <- rep(0,length(rds))
			ind <- 1 #this will be the index counter for vec
			for(x in rds)
			{
				vec[ind] <- ((nrow(dat[dat[,1]==m & dat[,3]>=x,]))/(nrow(dat[dat[,1]==m,])))*100 #getting the coverage percentage per chromosome
				ind <- ind + 1
			}
			
			if(m == cis)
			{
				chromosomes <- rep(paste("cis_chr_",m,sep=""),length(rds))
			}
			else
			{
				chromosomes <- rep(paste("chr_",m,sep=""),length(rds))
			}
			df_chrs_temp <- data.frame(chromosomes,rds,vec)
			df_chrs_temp <- df_chrs_temp[fi:num,]
			df_chrs <- rbind(df_chrs,df_chrs_temp)
		}

		df <- rbind(df_all,df_trans,df_chrs)
		rownames(df) <- seq(1,nrow(df),1)
		#creating the line plot
		print(ggplot2::ggplot(data=df, ggplot2::aes(x=rds, y=vec,group=chromosomes,color=chromosomes)) + ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::ggtitle(paste("minimum number reads/p-score VS. coverage line-plot for",file.name)) + ggplot2::labs(x="minimum reads/p-score",y="coverage(%)") + ggplot2::theme(plot.title = ggplot2::element_text(size=9)) + ggplot2::scale_x_continuous(breaks = seq(0,(max(rds)+1))))
				
		#save plots
		ans1 <- readline(prompt=cat("\nwould you like to save the plot?\ny/n\n\n"))
		if(ans1 == "y")
		{
			wd <- as.numeric(readline(prompt=cat("\nenter the width size of the plot:\n\n")))
			ht <- as.numeric(readline(prompt=cat("\nenter the height size of the plot:\n\n")))
			file.name2 <- strsplit(file.name,"[.]")
			ggplot2::ggsave(paste("~/Analyze4C/plots/",file.name2[[1]][1],"_coverageVSmin_linePlot_",num,"examples.png",sep=""),width=wd,height=ht)			
		}
		
		#try another number of examples
		ans2 <- readline(prompt=cat("\nwould you like to choose a different number of examples to use in the plot?\ny/n\n\n"))
		if(ans2 == "n")
		{
			flag1 <- 0
		}
	}	
}
