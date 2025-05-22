#' @title Extract loadings
#' @description  Extract loadings (species scores) from a gllvm object.
#' 
#' @param object an object of class 'gllvm'.
#' @param rotate whether to rotate the loadings based on the singular value decomposition of the site scores
#' @param scale whether to scale the loadings in accordance with the scaling rules of the "ordiplot" function
#' @param ... not used
#' 
#' @details
#' Function retrieves the loadings a.k.a. species scores for a GLLVM. For the optima of a quadratic response model, see \code{\link{optima.gllvm}}
#'@aliases getLoadings getLoadings.gllvm
#'@method getLoadings gllvm
#'@export
#'@export getLoadings.gllvm

getLoadings.gllvm <- function(object, scale = FALSE, rotate = FALSE, ...)
{
  if(inherits(object, "gllvm.quadratic"))stop("For the quadratic response model, please use 'optima.gllvm' instead.")
  sploads <- object$params$theta
  
  if(object$num.lv>0)sploads[,(ncol(sploads)-object$num.lv+1):ncol(sploads)] <- sploads[,(ncol(sploads)-object$num.lv+1):ncol(sploads)]%*%diag(tail(object$params$sigma.lv, object$num.lv))
  # if(object$num.lv.c>0 & scale)sploads[,1:object$num.lv.c] <- sploads[,1:object$num.lv.c]%*%diag(head(object$params$sigma.lv, object$num.lv.c))
  
  lvs <- getLV.gllvm(object) # for rotating/scaling loadings
  
  if (rotate) {
    svd_rotmat <- svd(lvs)$v # rotation matrix for sites
    sploads <- sploads %*% svd_rotmat # rotate loadings
  }
  
  if (scale) { # scale using Bert's mystery formula
    bothnorms <- vector("numeric",ncol(sploads))
    for(i in 1:ncol(sploads)){
      bothnorms[i] <- sqrt(sum(lvs[,i]^2)) * sqrt(sum(sploads[,i]^2))
    }
    
    ## Standardize to unit norm then scale using bothnorms. Note alpha = 0.5 so both have same norm. Otherwise "significance" becomes scale dependent
    scaled_cw_species <- sploads
    for(i in 1:ncol(scaled_cw_species)){
      scaled_cw_species[,i] <- sploads[,i] / sqrt(sum(sploads[,i]^2)) * (bothnorms[i]^(1-0.5)) 
    }
    
    sploads <- scaled_cw_species
  }
  
  return(sploads)
}

#'@export getLoadings
getLoadings <- function(object, ...)
{
  UseMethod(generic = "getLoadings")
}