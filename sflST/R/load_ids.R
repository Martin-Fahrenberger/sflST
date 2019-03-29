#' Building Spatial Transcriptomics coordinate object
#' @details This function will build a sparse fused lasso for each input gene from it's expression vector and the associated spatial coordinates and save the resulting models to an .Rdata object
#' @param file_name path to the Spatial Transcriptomics id file
#' @param output_folder where to save the .RDS object
#' @return does not return output but saves to files
#' @export
load_ids <- function(file_name,output_folder){

    ids_table <-read.table(file_name)
    rownames(ids_table) <- ids_table$V1
    ids_table <- ids_table[,c(2,3)]
    saveRDS(ids_table,file=paste(output_folder,'ids_table.RDS',sep='/'))

}
