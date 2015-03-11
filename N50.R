# v0.1 - CA Hamm - University of Kansas


library("Biostrings")


contigStats <- function(N, reflength, style = "data", pch = 19, 
	xlab="Percentage of Assembly Covered by Contigs of Size >= Y", 
	ylab = "Contig Size [bp]", main = "Cumulative Length of Contigs", 
	sizetitle = 14, sizex = 12, sizey = 12, sizelegend = 9, 
	xlim, ylim) {
        Nl <- lapply(names(N), function(x) rev(sort(N[[x]])))
        names(Nl) <- names(N)
        Nlcum <- lapply(names(Nl), function(x) cumsum(Nl[[x]]))
        names(Nlcum) <- names(Nl)
        N50 <- sapply(seq(along = N), function(x) Nl[[x]]
        [which(Nlcum[[x]] - reflength[x]/2 >= 0)[1]])
        names(N50) <- names(N)
    if(style == "data") {
        N75 <- sapply(seq(along = N), function(x) Nl[[x]]
        	[which(Nlcum[[x]] - reflength[x] * 0.75 >= 0)[1]])
        names(N50) <- names(N)
        N25 <- sapply(seq(along = N), function(x) Nl[[x]]
        	[which(Nlcum[[x]] - reflength[x] * 0.25 >= 0)[1]])
        names(N50) <- names(N)
    stats <- cbind(N25, N50, N75, Longest = sapply(N, max), 
        Mean = sapply(N, mean), Median = sapply(N, median), Shortest = sapply(N, min), N_Contigs = sapply(N, length))
                return(c(Nlcum, Contig_Stats = list(stats)))
        }
	if(style == "plot") {
            if(missing(xlim)) xlim <- c(0, 100)
            if(missing(ylim)) ylim <- c(0, max(unlist(N))) 
            split.screen(c(1,1))
            for(i in seq(along = Nl)) {
                    if(i == 1) {
            plot(Nlcum[[i]]/reflength[[i]] * 100, Nl[[i]], col = i, 
            pch = pch, xlim = xlim, ylim = ylim, xlab = xlab, 
            ylab = ylab, main = main)  
                    } 
            screen(1, new = FALSE)
            plot(Nlcum[[i]]/reflength[[i]] * 100, Nl[[i]], col = i, 
            pch = pch, xlim = xlim, ylim = ylim, xaxt = "n", 
            yaxt = "n", ylab = "", xlab = "", main = "", bty = "n")  
            }
            legend("topright", legend = paste(names(N50), ": N50 = ", 
            N50, sep = ""), cex = 1.2, bty = "n", pch = 19, 
            pt.cex = 1.2, col = seq(along = Nl)) 
            close.screen(all = TRUE)
        } 
}


#use as follows
genome <- readDNAStringSet("", format="fasta", use.names = FALSE)
N1 <- list(assembly1 = width(genome))
reflength1 <- sapply(N1, sum)

stats1 <- contigStats(N = N1, reflength = reflength1, style = "data") 
stats1[["Contig_Stats"]]

stats2 <- contigStats(N = N1, reflength = reflength1, style = "plot") 

