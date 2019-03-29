#' Plot Spatial Transcriptomics expression profile for one or more genes
#' @details This ggplot2-style function takes a count_matrix and a coordinate_matrix as well as a list of genes whose expression to visualize in 2D using one color per gene
#' @param counts_matrix input matrix of counts with genes in columns and spot identifier in rows
#' @param coordinate_matrix coordinate table which converts the spot identifier into an x and y coordinate on a grid
#' @param gene_set names of genes to plot, one or more
#' @param only_max TURE/FALSE wether to only plot the gene with the highest relative expression for each spot or all gene-expressions on top of each other
#' @param z_score TRUE/FALSE wether to normalize expression of genes withint the gene_set by z_score or by dividing by the max expression across al spots
#' @return ggplot2 plot object
#' @export

genewise_ggplot <- function(counts_matrix,gene_set,coordinate_matrix,only_max=TRUE,z_score=TRUE){
    require(ggplot2)
    hues = seq(15, 375, length = length(gene_set) + 1)
    gg_colours <- hcl(h = hues, l = 65, c = 100)[1:length(gene_set)]
    if (length(gg_colours)==1)  gg_colours <- "#000000"
    if (length(gene_set)>1){
        if (z_score) {
            normalized_counts <- apply(counts_matrix[,gene_set], MARGIN = 2,function(x) x/sd(x))
        }else{
            normalized_counts <- apply(counts_matrix[,gene_set], MARGIN = 2,function(x) x/max(x))
        }
        max_expression <- apply(normalized_counts[,gene_set], MARGIN = 1, max)
        which_expression <- apply(normalized_counts[,gene_set], MARGIN = 1, which.max)
        ggplot_data <- data.frame(cbind(x_coord=coordinate_matrix[,1],y_coord=coordinate_matrix[,2],expression=max_expression,which_max = which_expression))
        ggplot_data$which_max <- as.factor(ggplot_data$which_max)
    }else{
        if (z_score) {
            normalized_counts <- counts_matrix[,gene_set]/sd(counts_matrix[,gene_set])
        }else{
            normalized_counts <- counts_matrix[,gene_set]/max(counts_matrix[,gene_set])
        }
        ggplot_data <- data.frame(cbind(x_coord=coordinate_matrix[,1],y_coord=coordinate_matrix[,2],expression=normalized_counts,which_max = 1))
        ggplot_data$which_max <- as.factor(ggplot_data$which_max)
    }

    plot_1 <- ggplot(ggplot_data[ggplot_data$expression!=0,],aes(x=x_coord,y=y_coord,size=expression,colour=which_max)) + geom_point() + scale_color_manual(values = gg_colours,labels = gene_set,name='Genes:') + theme_bw() + guides(size=FALSE) + theme(plot.margin = unit(c(0,0,0,0),'cm'),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y = element_blank(),legend.key.width = unit(.1,'line'))
    return(plot_1)
}
