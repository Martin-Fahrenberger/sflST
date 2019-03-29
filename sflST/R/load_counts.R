#' load all count-files from a Spatial Transcriptomics data-set that was processed using simpleST
#' @details This function will load the count-files from the 1088 separate sequencing experiments as processed by simpleST and save the resulting count matrix as a .RDS file in the output-folder
#' @param alignments_folder path to the parent folder containing one sub-folder per alignment named by their spatial barcodes
#' @param output_folder path where to save the R object containing the count_matrix
#' @return does not return output but saves to files
#' @export
load_counts <- function(alignments_folder,output_folder) {
    require(data.table)
    start_time <- Sys.time()
    spot_dirs<- list.dirs(alignments_folder,recursive = FALSE)[-1]
    cat('Found',length(spot_dirs),'spot folders \n',sep=' ')
    succesfully_aligned<- unlist(lapply(spot_dirs,function (x) length(list.files(x)[!(grepl(list.files(x),pattern='_STARtmp'))])==16))
    cat('Of which',sum(succesfully_aligned),'contained successful alignments \n',sep=' ')
    spot_dirs<- spot_dirs[succesfully_aligned]
    spot_names<- list.dirs(alignments_folder,full.names = FALSE,recursive = FALSE)[-1]
    spot_names<- spot_names[succesfully_aligned]
    spot_files<- paste(spot_dirs,'/counts.txt',sep = '')

    counts_UMI<- c()
    for (i in spot_files){
        #cat(i,'\n')
        current_file <- as.data.frame(fread(i,header = TRUE))
        current_spot <- current_file[,7]
        names(current_spot) <- current_file$Geneid
        counts_UMI<- rbind(counts_UMI,current_spot)
    }
    rownames(counts_UMI) <- spot_names

    counts_final <- counts_UMI[,colSums(counts_UMI)!=0]
    counts_final <- counts_final[,which(!(colnames(counts_final) %in% c('Malat1','mt-Rnr1','mt-Rnr2','Rn18s-rs5')))]
    cat('Saving count matrix to',paste(output_folder,'counts_final.RDS \n',sep='_'))
    system(paste('mkdir -p',output_folder,sep=' '))
    saveRDS(counts_final,file=paste(output_folder,'counts_final.RDS',sep='/'))
    end_time <- Sys.time()
    cat('Build count matrix containing',dim(counts_final)[1],'spots and',dim(counts_final)[2],'genes in',end_time - start_time,'minutes \n')
}
