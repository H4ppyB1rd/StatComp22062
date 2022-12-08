#' @title l2 Norm Sq
#' @description we use this function to calculate the square of l2 norm of a vector.
#' @param b a vector
#' @examples
#' \dontrun{
#' x = c(1:3)
#' print(l2NormSq(x))
#' }
#' @export
l2NormSq = function(b){
  return( crossprod(b,b)[1,1] )
}
