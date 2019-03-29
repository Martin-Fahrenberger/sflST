#' Building sparse fused lasso models for given genes from their count data
#' @details This function will build a sparse fused lasso for each input gene from it's expression vector and the associated spatial coordinates and save the resulting models to an .Rdata object
#' @param gene_list a set of genen names for which to build the models
#' @param counts_matrix a matrix with gene expression levels across genes (columns) and associated spots (rows)
#' @param ids_table table  associating each rowname of the counts_matrix with an x and y coordinate
#' @param name name of the object to save to
#' @param output_folder folder to save object in
#' @param include_diag logical value wether to include diagonal spots as neigbours when computing the lasso, Warning slow if active (default FALSE)
#' @param gamma ratio of penalization between sparsity and continuity (default 1, high = sparse)
fused_lasso_complete_fixed_gamma <- function(gene_list,counts_matrix,ids_table,name,output_folder,include_diag=FALSE,gamma=1){
    require(genlasso)
    source('./Chunks/fusedlasso2d_diag.R')
    lasso_set <- c()
    for (i in gene_list) {
        tmp_matrix <- matrix(data = 0,nrow = max(ids_table[rownames(counts_matrix),1]),ncol = max(ids_table[rownames(counts_matrix),2]))
        for (j in 1:length(counts_matrix[,i])){
            tmp_matrix[ids_table[names(counts_matrix[,i][j]),1],ids_table[names(counts_matrix[,i][j]),2]] <- counts_matrix[,i][j]
        }
        if (include_diag == TRUE) {
            tmp_lasso <- list(fusedlasso2d_diag(tmp_matrix , gamma = gamma))
        } else {
            tmp_lasso <- list(fusedlasso2d(tmp_matrix,gamma = gamma))
        }
        lasso_set <- append(lasso_set,tmp_lasso)
        cat('X')
    }
    cat(paste('\n Saving ',output_folder,'/',name,'.RData to disk \n',sep = ''))
    names(lasso_set) <- gene_list
    save(lasso_set,file = paste(output_folder,'/',name,'.RData',sep = ''))
}
