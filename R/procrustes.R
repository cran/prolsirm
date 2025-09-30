##########################################################################
## function that performs procrustes matching on posterior samples of latent
## positions from latent space item response model
##
## returns the matched samples of latent positions after procrustes transformation
##
## The matching is based on Jeon, M., Jin, I. H., Schweinberger, M., & Baugh, S. (2021).
##    Mapping Unobserved Item–Respondent Interactions: A Latent Space Item Response Model with Interaction Map.
##    Psychometrika, 86(2), 378–403. https://doi.org/10.1007/s11336-021-09762-5
##
## The procrustes transformation is based on Borg and Groenen 1997. Modern Multidimensional
##   Scaling. New York: Springer. pp. 340-342.
##
##
## Jinwen Luo
## UCLA
## 8/28/2023
##
## Copyright (C) 2023-present Jinwen Luo, Minjeong Jeon.
##########################################################################

#' Procrustes matching of posterior samples of latent positions of persons ( \code{z}) and items (\code{w}) from LSIRM with \code{q}-dimensional interaction map
#'
#' The function performs Procrustes matching by aligning the \code{M}-1 posterior samples of \code{z} and \code{w} to the reference configuration (e.g., iteration 1). Users can select any configuration as the reference with the \code{ref} argument. The function returns \code{M} matched lists of \code{z} and \code{w} after the Procrustes matching process.
#'
#' @title procrustes
#' @param z A list of length \code{M}. \code{M} posterior samples of person latent positions \code{N} x \code{q}, where \code{N} is the number of respondents and \code{q} is the dimension of the interaction map.
#' @param w A list of length \code{M}. \code{M} posterior samples of item latent positions  \code{I} x \code{q}, where \code{I} is the number of items  and \code{q} is the dimension of the interaction map.
#' @param ref Reference configuration (i.e., iteration) index. Default is 1 (i.e., posterior samples at iteration 1) ).
#' @param dilation Logical; allow *uniform scaling* during alignment
#'   (default \code{FALSE}). Set \code{TRUE} only for plotting/overlay use cases.
#' @return A list of \code{M} matched posterior samples of \code{z}  and a list of \code{M} matched posterior samples of \code{w}.
#' @importFrom MCMCpack procrustes
#' @examples
#' # Load package
#' library(prolsirm)
#'
#' # Generate example posterior samples
#' # M=3 samples of person latent positions (N=50, q=2)
#' z <- list(matrix(rnorm(100), ncol = 2),
#'           matrix(rnorm(100), ncol = 2),
#'           matrix(rnorm(100), ncol = 2))
#'
#' # M=3 samples of item latent positions (I=5, q=2)
#' w <- list(matrix(rnorm(10), ncol = 2),
#'           matrix(rnorm(10), ncol = 2),
#'           matrix(rnorm(10), ncol = 2))
#'
#' # Perform Procrustes matching
#' matched_data <- procrustes(z = z, w = w)
#'
#' @export
#'
#' @references Borg and Groenen. 1997. \emph{Modern Multidimensional Scaling}.
#' New York: Springer. pp. 340-342.
#' @references Jeon, M., Jin, I. H., Schweinberger, M., & Baugh, S. 2021. \emph{Mapping Unobserved Item–Respondent Interactions: A Latent Space Item Response Model with Interaction Map}. Psychometrika, 86(2), 378–403.
#' @references Andrew D. Martin, Kevin M. Quinn, Jong Hee Park. 2011. MCMCpack: Markov Chain Monte Carlo in R. Journal of Statistical Software. 42(9): 1-21.
#'
procrustes  <-function(z,w,ref=1,dilation=FALSE) {

        # take the z, w as input
        # z is a M * N * q list object containing the latent coordinates of row objects.
        # M is the length of list, or the number of latent space configurations to be matched.
        # N is the number of row objects, and two coordinates for two dimensional latent space
        # w is a M * I * q list object containing the latent coordinates of column objects.
        # I is the number of columns, , and two coordinates for two dimensional latent space
        # the default reference for procrustes matching will be the first configurations in the list
        M.w <- length(w)
        M.z <- length(z)
        if (M.w != M.z) {stop(sprintf("Person and item posterior sample lengths do not match!"))}else{
              M <- M.z
        }
        zz0 <- z[[ref]]
        ww0 <- w[[ref]]
        # define sample sizes
        nz = nrow(zz0)
        nw = nrow(ww0)
        if (ncol(zz0) != ncol(ww0)) {stop(sprintf("Person and item dimensions do not match!"))}
        # center the reference wz matrix
        wz_ref=rbind(z[[ref]],w[[ref]])
        mean_pos=apply(wz_ref,2,mean)
        wz0=sweep(wz_ref,2,mean_pos)
        # procrustes rotations through MCMCpack
        matched_zz = z
        matched_ww = w

        # R package
        # require(MCMCpack)
        for (ii in 1: M) {
                # target matrix: wz0 (centered)
                #center first
                wz=rbind(z[[ii]],w[[ii]])
                mean_pos=apply(wz,2,mean)
                wz_centered=sweep(wz,2,mean_pos)

                proc=MCMCpack::procrustes(wz_centered, wz0, translation=TRUE, dilation=dilation)

                wz_matched0  <- proc$X.new

                # results
                matched_zz[[ii]]=wz_matched0[1:nz,]
                matched_ww[[ii]]=wz_matched0[(nz+1):(nz+nw),]
        }

        return( list(zz=matched_zz, ww=matched_ww))
}
