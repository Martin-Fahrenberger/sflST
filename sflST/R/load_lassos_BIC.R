#' load all lasso models optioimized by BIC
#' @details This function loads all models from a folder as built by the build_lassos function optimizes them by BIC and saves a counts object containing a count_matrix, a fits object all the fitted models, a BICs object containing the BIC value for each model and a lls object containing the loglikelihood value of each model.
#' @param lasso_folder path to the directory containing the lasso R objects from the build_lassos function
#' @param name to append to the save files
#' @return does not return output but saves to files in working directory
#' @export
load_lassos_BIC <- function(lasso_folder,name){
    library(genlasso)
    source('Chunks/coef_error.R')
    source('Chunks/calc_BICs.R')
    all_lasso_files <- list.files(lasso_folder,full.names=T)
    counts_mob_gamma_1 <- c()
    fits_mob_gamma_1 <- c()
    BICs_mob_gamma_1 <- c()
    lls_mob_gamma_1 <- c()
    for (i in all_lasso_files){
        load(i)
        tmp_fits_set <- c()
        fits_mob_set <- c()
        BICs_mob_set <- c()
        lls_mob_set <- c()
        for (j in names(lasso_set)){
            BIC_list <- calc_BICs(lasso_set[[j]])
            lambda_BIC <- BIC_list[which.min(BIC_list[,4]),2]
            fits_mob <- coef_error(lasso_set[[j]], lambda=lambda_BIC)
            tmp_fit <- fits_mob$beta
            tmp_fits_set <- cbind(tmp_fits_set,tmp_fit)
            fits_mob_set <- append(fits_mob_set,list(fits_mob))
            BICs_mob_set <- append(BICs_mob_set,min(BIC_list[,4]))
            lls_mob_set <- append(lls_mob_set,min(BIC_list[,5]))
        }
        colnames(tmp_fits_set) <- names(lasso_set)
        counts_mob_gamma_1 <- cbind(counts_mob_gamma_1,tmp_fits_set)
        fits_mob_gamma_1 <- append(fits_mob_gamma_1,fits_mob_set)
        BICs_mob_gamma_1 <- append(BICs_mob_gamma_1,BICs_mob_set)
        lls_mob_gamma_1 <-  append(lls_mob_gamma_1,lls_mob_set)
    }
    saveRDS(Matrix(counts_mob_gamma_1,sparse=TRUE),file=paste('./',name,'/','counts_lasso_BIC.RDS',sep = ''))
    saveRDS(fits_mob_gamma_1,file=paste('./',name,'/','fits_lasso_BIC.RDS',sep = ''))
    saveRDS(BICs_mob_gamma_1,file=paste('./',name,'/','BICs_lasso_BIC.RDS',sep = ''))
    saveRDS(lls_mob_gamma_1,file=paste('./',name,'/','lls_lasso_BIC.RDS',sep = ''))

}
