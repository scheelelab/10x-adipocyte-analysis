# R
# Nagendra Palani - Scheele lab NNF CBMR
library(dplyr)
library(ggplot2)


setwd("/home/projects/ku_00021/data/preadipo_scrnaseq")

df_branch <- readRDS("df_branch.RDS")
secgenes <- readRDS("secgenes.RDS")

xx <- length(secgenes)




plot_pseudotime <- function(dataf, gene, span, diffvalue){
  p <- ggplot(dataf, aes_string(x='Pseudotime', y=gene)) +
    geom_point(shape = 21, colour = "black", size = 1, stroke = 0, aes(fill=State.labels), show.legend=F) +
    geom_smooth(se = FALSE, aes(color=Branch), span=span, method='loess', size=2, show.legend=F)
    #set the correct aes for geom_smooth

  
  g <- ggplot_build(p) #gets the plotting data from the ggplot object
  subvalue <- g$data[[2]]$y[g$data[[2]]$group == 1] - g$data[[2]]$y[g$data[[2]]$group == 2]
  differences <- abs(subvalue)
  #g$data[[1]] contains the data point, g$data[[2]] the smoothed means
  pseudotime = g$data[[2]]$x[g$data[[2]]$group == 1]
  gene_diff <- data.frame(pseudotime, differences,subvalue)
  
  
  if (is.finite(min(gene_diff$pseudotime[gene_diff$differences > diffvalue]))) {
  diverge_point_ptmin <- min(gene_diff$pseudotime[gene_diff$differences > diffvalue])
  diverge_point_ptmax <- max(gene_diff$pseudotime[gene_diff$differences > diffvalue])
  diverge_point_sv  <- min(gene_diff$subvalue[gene_diff$differences > diffvalue])
     
  } else
  {
     diverge_point_ptmin <- 0
  diverge_point_ptmax <- 0
  diverge_point_sv  <- 0
  }
  
  xx <- list(gene,diverge_point_ptmin,diverge_point_ptmax,diverge_point_sv)
  return(xx)
}


outlist <- list()

 df_depot <- df_branch %>% dplyr::filter(tissue == "brown") # filter by tissue 
 
 # for white, change the filter value to "white", 
            # change name of the output file below, and rerun script.

 for(i in 1:xx)  {
   #set.seed(2020)

outlist[[i]] <- plot_pseudotime(df_depot, secgenes[i], span=0.9, 0.15) # diffvalue is the minimum difference in expression value

  }



outdf  <-  as.data.frame(matrix(unlist(outlist), nrow=length(unlist(outlist[1])))) %>% t
rownames(outdf) = NULL
outdf <- as.data.frame(outdf)

 outdf_net <- outdf %>% filter(V2 > 0 & V3 >98) # stretched pseudotime value cutoffs.

data.table::fwrite(outdf_net,file = "secreted_genes_brown_diverge.csv")

quit(save= "no")
