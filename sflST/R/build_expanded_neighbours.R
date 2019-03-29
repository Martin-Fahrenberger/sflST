#' Build 2D neighbourhood matrix including diagonal neigbours
#' @details Builds a 2D neighbourhood matrix as used as D in the genlasso function
#' @param dim1 grid size in x direction
#' @param dim2 grid size in y direction
#' @return returns a sparse neighbourhood matrix
build_expanded_neighbours <- function(dim1,dim2){
    D12 <- getD2dSparse(dim1,dim2)
    s2 <- 71

    D3=Matrix(0,nrow = dim1*dim2,ncol=dim1*dim2)
    for (i in 1:(dim1*dim2)){ #east neigbours for all rows
        D3[i,i] <- -s2
        if (i+1+dim1 <= dim1*dim2) D3[i,i+1+dim1] <- s2
    }
    D3 = D3[1:((dim1-1)*(dim2-1)), ]

    D4=Matrix(0,nrow = dim1*dim2,ncol=dim1*dim2)
    for (i in 1:(dim1*dim2)){ #east neigbours for all rows
        D4[i,i] <- -s2
        if (i+dim1-1 <= dim1*dim2 ) D4[i,i+dim1-1] <- s2
    }
    D4 = D4[(seq(0, dim1 * dim2-1)%%dim1) != 0, ]
    D4 = D4[1:((dim1-1)*(dim2-1)), ]
    #return(Matrix(rbind(D1,D2,D3,D4),sparse=TRUE))
    return(Matrix(rbind(D12,D3,D4),sparse=TRUE))
}
