#' @export

#venn_DF recieves a data frame of ones and zeros representing the number intersected for each file
#num_of_files represents the number of files (number of columns in the data frame)
#extra tells us if there is an extra file that should be added to num_of_files, this is usually the reference file which will all be ones, this is when extra equals 1

venn_creator <- function(venn_DF,num_of_files,extra=0,vennName = -1)
{
	#check how many files there are
	if(num_of_files + extra > 5) #there can be up to 5 sets (due to VennDiagram package constraints)
	{
		cat("\nthere are too many files in order to create a venn diagram\n\n")
	}
	else
	{
		#choosing names for each file
		cat("\nchoose names for each one of the files:\n\n")
		nms <- c()
		for(i in 1:(num_of_files + extra))
		{
			tmp <- as.character(readline(prompt=cat("enter a name for file no.",i,":\n\n",sep="")))
			nms <- c(nms,tmp)
		}
		
		cat("\ncreating the venn diagram\nplease wait...\n\n")
		
		grid::grid.newpage()
		#creating the venn diagrams according to number of files	
		switch((num_of_files + extra),
		#one file
		{
			venn.plot <- VennDiagram::draw.single.venn(area1=sum(venn_DF[,1]),
				category = nms, lty = "blank", fill = "skyblue")
		},
		#two files
		{
			venn.plot <- VennDiagram::draw.pairwise.venn(area1=sum(venn_DF[,1]),area2=sum(venn_DF[,2]),
				n12=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1),
				category = nms, lty = "blank", fill = c("skyblue", "pink1"))
		},
		#three files (7 combinations)
		{
			venn.plot <- VennDiagram::draw.triple.venn(area1=sum(venn_DF[,1]),area2=sum(venn_DF[,2]),area3=sum(venn_DF[,3]),
				n12=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1),
				n13=sum(venn_DF[,1] == 1 & venn_DF[,3] == 1),
				n23=sum(venn_DF[,2] == 1 & venn_DF[,3] == 1),				
				n123=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,3] == 1),
				category = nms, lty = "blank", fill = c("skyblue", "pink1", "mediumorchid"))	
		},
		#four files (15 combinations)
		{
			venn.plot <- VennDiagram::draw.quad.venn(area1=sum(venn_DF[,1]),area2=sum(venn_DF[,2]),area3=sum(venn_DF[,3]),area4=sum(venn_DF[,4]),
				n12=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1),
				n13=sum(venn_DF[,1] == 1 & venn_DF[,3] == 1),
				n14=sum(venn_DF[,1] == 1 & venn_DF[,4] == 1),
				n23=sum(venn_DF[,2] == 1 & venn_DF[,3] == 1),
				n24=sum(venn_DF[,2] == 1 & venn_DF[,4] == 1),
				n34=sum(venn_DF[,3] == 1 & venn_DF[,4] == 1),				
				n123=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,3] == 1),
				n124=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,4] == 1),
				n134=sum(venn_DF[,1] == 1 & venn_DF[,3] == 1 & venn_DF[,4] == 1),
				n234=sum(venn_DF[,2] == 1 & venn_DF[,3] == 1 & venn_DF[,4] == 1),
				n1234=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,3] == 1 & venn_DF[,4] == 1),
				category = nms, lty = "blank",fill = c("skyblue", "pink1", "mediumorchid","orange"))
		},
		#five files	(31 combinations)
		{
			venn.plot <- VennDiagram::draw.quintuple.venn(area1=sum(venn_DF[,1]),area2=sum(venn_DF[,2]),area3=sum(venn_DF[,3]),area4=sum(venn_DF[,4]),area5=sum(venn_DF[,5]),
				n12=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1),
				n13=sum(venn_DF[,1] == 1 & venn_DF[,3] == 1),
				n14=sum(venn_DF[,1] == 1 & venn_DF[,4] == 1),
				n15=sum(venn_DF[,1] == 1 & venn_DF[,5] == 1),
				n23=sum(venn_DF[,2] == 1 & venn_DF[,3] == 1),
				n24=sum(venn_DF[,2] == 1 & venn_DF[,4] == 1),
				n25=sum(venn_DF[,2] == 1 & venn_DF[,5] == 1),
				n34=sum(venn_DF[,3] == 1 & venn_DF[,4] == 1),
				n35=sum(venn_DF[,3] == 1 & venn_DF[,5] == 1),
				n45=sum(venn_DF[,4] == 1 & venn_DF[,5] == 1),		
				n123=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,3] == 1),
				n124=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,4] == 1),
				n125=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,5] == 1),
				n134=sum(venn_DF[,1] == 1 & venn_DF[,3] == 1 & venn_DF[,4] == 1),
				n135=sum(venn_DF[,1] == 1 & venn_DF[,3] == 1 & venn_DF[,5] == 1),
				n145=sum(venn_DF[,1] == 1 & venn_DF[,4] == 1 & venn_DF[,5] == 1),
				n234=sum(venn_DF[,2] == 1 & venn_DF[,3] == 1 & venn_DF[,4] == 1),
				n235=sum(venn_DF[,2] == 1 & venn_DF[,3] == 1 & venn_DF[,5] == 1),
				n245=sum(venn_DF[,2] == 1 & venn_DF[,4] == 1 & venn_DF[,5] == 1),
				n345=sum(venn_DF[,3] == 1 & venn_DF[,4] == 1 & venn_DF[,5] == 1),
				n1234=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,3] == 1 & venn_DF[,4] == 1),
				n1235=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,3] == 1 & venn_DF[,5] == 1),
				n1245=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,4] == 1 & venn_DF[,5] == 1),
				n1345=sum(venn_DF[,1] == 1 & venn_DF[,3] == 1 & venn_DF[,4] == 1 & venn_DF[,5] == 1),
				n2345=sum(venn_DF[,2] == 1 & venn_DF[,3] == 1 & venn_DF[,4] == 1 & venn_DF[,5] == 1),
				n12345=sum(venn_DF[,1] == 1 & venn_DF[,2] == 1 & venn_DF[,3] == 1 & venn_DF[,4] == 1 & venn_DF[,5] == 1),
				category = nms, lty = "blank",fill = c("skyblue", "pink1", "mediumorchid","orange","green"))
		}
		)
		
		print(venn.plot)
		
		#asking if to save
		if(readline(prompt=cat("\nwould you like to save this venn diagram?\ny/n\n\n")) == "y")
		{
			#save the venn diagram
			if(vennName == -1)
			{
				#getting the date and time
				DandT1 <- toString(Sys.time())
				DandT2 <- gsub(" ","_",DandT1)
				DandT2 <- gsub(":","",DandT2)
				vennName <- paste("venn_",DandT2,".jpg",sep="")
			}
			
			jpeg(paste("~/Analyze4C/plots/",vennName,sep=""))
			grid.draw(venn.plot)
			dev.off()

			return(1) #states that the diagram was saved
		}	
	}
	return(0)
}
