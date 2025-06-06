#' @title Extract latent variables
#' @description  Extract latent variables from gllvm object.
#' 
#' @param object an object of class 'gllvm'.
#' @param type type of latent variable scores to retrieve from a gllvm object. For models with unconstrained latent variables, defaults to "residual". For models with constrained latent variables, defaults to conditional. Alternatively, "marginal" returns linear combination scores without residual error.
#' @param rotate whether to rotate the latent variable scores based on their singular value decomposition 
#' @param scale whether to scale the latent variable scores in accordance with the scaling rules of the "ordiplot" function
#' @param ... not used
#' 
#' @details
#' Function retrieves the site scores for a GLLVM. Each type corresponds to a separate term of the model. For a GLLVM with unconstrained latent variables
#' the default is "residual". "Residual" scores represent the error term in concurrent ordination, and are not available for constrained ordination.
#' 
#' For GLLVMs with informed latent variables, "conditional" returns the complete site scores, due to both fixed- and latent effects,
#' where the latent effect is always scaled by the diagonal of the species loadings so that it can be small relative to the fixed-effects. "Conditional"
#' here means conditional on the random-effect i.e. the residual.
#' 
#' Type "marginal"  returns linear combination scores, i.e. the site scores only due to fixed-effects. These are available for constrained and concurrent ordination.
#' 
#' If both unconstrained and constrained latent variables are included in the model, type "marginal" returns linear combination scores for constrained latent variables
#' but "residual" scores for unconstrained latent variables.
#'@aliases getLV getLV.gllvm
#'@method getLV gllvm
#'@export
#'@export getLV.gllvm

getLV.gllvm <- function(object, rotate = FALSE, scale = FALSE, type = NULL, ...)
{
  n <- nrow(object$y)
  
  if((object$num.RR+object$num.lv.c+object$num.lv)==0){
    stop("No latent variables in model.")
  }
  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X #for backward compatibility
  if(!is.null(type)){
    if(!type%in%c("residual","conditional","marginal")){
    stop("Type should be one of: residual, conditional, or marginal.")
    }
    if(type=="conditional"&object$num.lv.c==0){
      stop("Can only retrieve conditional site scores for concurrent ordination.")
    }else if(type=="residual"&object$num.lv.c==0&object$num.lv==0){
      stop("Can only retrieve residual site scores for for a random-effects ordination.")
    }else if(type=="marginal"&object$num.lv.c==0&object$num.RR==0){
      stop("Can only retrieve marginal site scores for an ordination that includes predictors.")
    }
  }
  
  if(is.null(type)){
  if(object$num.lv.c==0&object$num.RR==0){
    type = "residual"
  }else{
    if(object$num.lv.c>0){
      type = "conditional"
    }else if(object$num.RR>0){
      type = "marginal"
    }
    
  }
  } 

  if(type == "residual"){
    lvs <- object$lvs
  }else if(type=="conditional"&object$num.lv.c>0){
      lvs <- object$lvs
      lvs[,1:object$num.lv.c] <- t(t(lvs[,1:object$num.lv.c,drop=F])*object$params$sigma.lv[1:object$num.lv.c])
      lvs[,1:object$num.lv.c] <- lvs[,1:object$num.lv.c]+object$lv.X.design%*%object$params$LvXcoef[,1:object$num.lv.c,drop=F]
      if(object$num.RR>0){
        lvs <- cbind(lvs[,1:object$num.lv.c,drop=F],object$lv.X.design%*%object$params$LvXcoef[,-c(1:object$num.lv.c),drop=F], lvs[,-c(1:object$num.lv.c),drop=F])
      }
  }else if(type=="marginal"){
    lvs <- object$lv.X.design%*%object$params$LvXcoef
  }
  
  #only allow residual or marginal with RR + num.lv, since num.lv is not uncorrelated with predictors
  if(type=="marginal"&object$num.lv>0)object$num.lv <- 0
  if(type=="residual"&object$num.RR>0)object$num.RR <- 0
  
  
  if((object$num.lv.c+object$num.RR)>0){
    if((object$num.lv.c+object$num.RR)>0&object$num.lv>0)colnames(lvs)<-c(paste("CLV",1:(object$num.lv.c+object$num.RR),sep=""),paste("LV",1:object$num.lv,sep=""))
    if((object$num.lv.c+object$num.RR)>0&object$num.lv==0)colnames(lvs)<-paste("CLV",1:(object$num.lv.c+object$num.RR),sep="")
  }
  
  if (rotate) {
    svd_rotmat_sites <- svd(lvs)$v # SVD rotation matrix
    lvs <- lvs %*% svd_rotmat_sites # rotate lvs NB: does the scaling affect the svd??
  }
  
  if (scale) { # scale using Bert's mystery formula
    # get species loadings scaled by sigma.lv
    choose.lv.coefs <- getLoadings(object)
    
    bothnorms <- vector("numeric",ncol(choose.lv.coefs))
    for(i in 1:ncol(choose.lv.coefs)){
      bothnorms[i] <- sqrt(sum(lvs[,i]^2)) * sqrt(sum(choose.lv.coefs[,i]^2))
    }
    
    ## Standardize both to unit norm then scale using bothnorms. Note alpha = 0.5 so both have same norm. Otherwise "significance" becomes scale dependent
    scaled_cw_sites <- t(t(lvs) / sqrt(colSums(lvs^2)) * (bothnorms^0.5)) 
    
    lvs <- scaled_cw_sites
  }
  
  
  return(lvs)
}

#'@export getLV
getLV <- function(object, ...)
{
  UseMethod(generic = "getLV")
}