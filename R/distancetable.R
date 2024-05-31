#' G0I Distances Table
#'
#' \code{distancesg0itable} G0I Distances Table
#'
#'
#' @param image an array image with dimensions X,Y,Z
#' @param marks1 vector 1 of marks with the follow order c(left, right, up, down)
#' @param marks2 vector 2 of marks with the follow order c(left, right, up, down)
#' @param p     distance Renyi parameter
#'
#' @return Distances table between samples obtained of marks.
#'
#' @references
#' \itemize{
#'   \item Nascimento, Abraao DC, Renato J. Cintra, and Alejandro C. Frery. "Hypothesis testing in speckled data with stochastic distances." IEEE Transactions on geoscience and remote sensing 48.1 (2009): 373-385.
#' }
#'
#' @examples
#'
#' set.seed(1) # determining a seed
#'
#' data("SanFrancisco150")
#' image <- SanFrancisco150
#'
#' # left,right, up,down)
#' marks1 <- c(10,30,50,100)
#' marks2 <- c(20,50,80,130)
#'
#' teste <- distancesg0itable(image,marks1, marks2, p = 1/2)
#' teste
#'
#'
#' @importFrom stats integrate
#'
#' @export
distancesg0itable = function(image, marks1, marks2, p = 1/2){
  UseMethod("distancesg0itable")
}

#' @export
distancesg0itable.default = function(image, marks1, marks2, p = 1/2){


  # --------------------------------------------------------
  if(length(marks1) != 4 | length(marks2) != 4 | any(marks1 == 0 ) | any(marks2 == 0 ) ){
    stop("Your need to redefine the marks1 and/or marks2")
  }
  # --------------------------------------------------------

  # Selecting the samples
  sample1 <- sapply(1:3, function(i) Re(image[marks1[3]:marks1[4], marks1[1]:marks1[2],i]))
  sample2 <- sapply(1:3, function(i) Re(image[marks2[3]:marks2[4], marks2[1]:marks2[2],i]))

  # Estimating Parameters
  estimatingsample1 <- lapply(1:3, function(i) g0iestimation(c(-1.1,1,1), sample1[,i]))
  estimatingsample2 <- lapply(1:3, function(i) g0iestimation(c(-1.1,1,1), sample2[,i]))

  #--------------
  # Calculating the distances
  vec <- c("KL", "Renyi", "Bhattacharyya", "Hellinger",
           "Arithmetic-Geometric", "Triangular", "Harmonic-Mean")
  #--------------
  cat("Calculating the distances")

  res <- t(sapply(1:3,function(j){
      sapply(1:length(vec), function(i){
        # --
        #cat(".")
        cat(paste(j,"-", i,vec[i],"|"))
        # --
        distancesg0i.default(estimatingsample1[[j]],estimatingsample2[[j]],
                             type = vec[i])$distance_value})
    }))
  #--------------
  colnames(res) <- vec
  rownames(res) <- c("HH", "HV", "VV")

  #-------------------------------------------------------
  structure(list(
    marks1 = marks1,
    marks2 = marks2,
    distancestable = res
  ),
  class = "distancesg0itable"
  )

}


#' @export
print.distancesg0itable <- function(x, ...) {
  # -----------------
  cat("\nDistances Table\n\n")
  # -----------------
  print(x$distancestable)
  # -----------------
}

