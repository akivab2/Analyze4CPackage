#' @export

contactBands <- function(filename,RE_gap,bp_gap,cis,co,genome_sizes)
{
    #first create bands of RE sites, they need to be without a break of more than RE_gap (the user decides how many) RE sites in between them and no more than bp_gap (the user decides) kb distance between them.
    #these bands are then compared between L1 and L2 using bedtools
    #since i am using distances in kb i will need to rearrange any of the files that the data has moved around in when moving the translocation
    #co represents the cutoff that above that number it will be considered positive
    
    dat <- filename[filename[,1]!=cis,] #here cis is removed before creating contact bands, if there is no need to remove cis then just remove this line
    #creating bands:
    
    #first gets the beginning of each band
    first <- 0
    
    #counts number of RE sites with 0 in a row
    zero_REs <- 0
    
    #has the previous site with positive RE site
    prev <- 0
    
    #contains all the contacts
    conts <- c()
    
    #gets the size of the data
    len <- nrow(dat)
    
    for(i in 1:len)
    {
		#if this is the first positive RE site dealt with in the genome
        if(dat[i,3]>co && first == 0)
        {
            first <- dat[i,]
            prev <- first
            last <- first
            zero_REs <- 0
            cur_chr <- first[1]
        }
        else if(all(first != 0)) #the 'all' makes it easier to test if 'first' is true
        {
            if(dat[i,3]>co && dat[i,1]!=cur_chr) #if we reached a different chromosome
            {
                band <- c(first[1],first[2],last[2])
                conts <- rbind(conts,band)
                first <- dat[i,]
                prev <- first
                last <- first
                zero_REs <- 0
                cur_chr <- first[1]
            }
            else #if we stay in the same chromosome
            {
				#if we reached a positive RE (above cutoff) and we are less than RE_gap from the last one and less than bp_gap bp from the last one
                if((dat[i,3]>co) && (zero_REs <= RE_gap) && ((dat[i,2]-prev[2]) <= bp_gap))
                {
                    last <- dat[i,]
                    prev <- last
                }
				#if we reached a positive RE (above cutoff) and we are more than RE_gap from the last one or more than bp_gap bp from the last one
                else if((dat[i,3]>co) && ((zero_REs > RE_gap) || ((dat[i,2]-prev[2]) > bp_gap)))
                {
					if((first[2] == last[2]) && ((first[2] + 1) <= genome_sizes[as.numeric(cur_chr),2]))
					{
						last[2] <- first[2] + 1
					}
                    band <- c(first[1],first[2],last[2])
                    conts <- rbind(conts,band)
                    first <- dat[i,]
                    prev <- first
                    last <- first
                    zero_REs <- 0
                }
				#if we reached a band that isn't above cutoff
                else if(dat[i,3] <= co)
                {
                    zero_REs <-  zero_REs + 1
                }
            }
        }
    }
    band <- c(first[1],first[2],last[2])
    conts <- rbind(conts,band)
    return(conts)
}
