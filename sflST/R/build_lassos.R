#' Calculate 2D fused-lasso solutionpath for a countmatrix and save as object
#' @details This function builds a solution path for each gene (column) in a input matrix and saves it as an .Rdata object
#' @param counts_matrix input matrix of counts with genes in columns and spot identifier in rows
#' @param coords_table coordinate table which converts the spot identifier into an x and y coordinate on a grid
#' @param name Name of the object to save to.
#' @param ncores number of parallel processes to start (default 8)
#' @param include_diag wether or not to include diagonal neighbours on the coordinate grid (default FALSE) Warning slow if activated.
#' @param output_folder folder in which to save the generated lasso models
#' @param gamma ratio of penalization between sparsity and continuity (default 1, high = sparse)
#' @return does not return output but saves to files
#' @export
build_lassos <- function(counts_matrix,coords_table,name,ncores=8,include_diag=FALSE,output_folder,gamma=1){
    require(genlasso)
    require(parallel)
    number_of_sets <- ceiling(length(colnames(counts_matrix))/100)
    list_of_sets <- c()
    for (i in 1:(number_of_sets-1)){
        list_of_sets <- append(list_of_sets,list(colnames(counts_matrix)[(1+(i-1)*100):(i*100)]))
    }
    list_of_sets <- append(list_of_sets,list(colnames(counts_matrix)[(1+(number_of_sets-1)*100):length(colnames(counts_matrix))]))
    #dim1 <- max(coords_table[rownames(counts_matrix),1])
    #dim2 <- max(coords_table[rownames(counts_matrix),2])
    #D_matrix <- build_expanded_neighbours(dim1, dim2)
    #D_matrix <- getD2dSparse(dim1,dim2)
    #mclapply(1:length(list_of_sets),function(x) fused_lasso_complete_fixed_gamma(list_of_sets[[x]],counts_matrix,D_matrix = D_matrix,coords_table,paste(name,x,sep='_'),include_diag,output_folder),mc.cores = ncores)
    #mclapply(1:2,function(x) fused_lasso_complete_fixed_gamma(list_of_sets[[x]],counts_matrix,D_matrix = D_matrix,coords_table,paste(name,x,sep='_'),include_diag,output_folder),mc.cores = ncores)
    mclapply(1:length(list_of_sets),function(x) fused_lasso_complete_fixed_gamma(list_of_sets[[x]],counts_matrix,coords_table,name=paste(name,x,sep='_'),output_folder=output_folder,include_diag,gamma=gamma),mc.cores = ncores)
}
