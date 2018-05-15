#' AnimalUni that extends \code{"stanfit"} class.
#'
#' AnimalUni provides some resonable print and summary methods for animal model fits
#'
#' @name animalUni-class
#' @rdname animalUni-class
#' @export
amimalUni = setClass("animalUni", contains = "stanfit")

#'@export
setMethod("show", "animalUni", 
          function(object) print(object, pars = c("sigma_G", "sigma_E", "beta")))
#'@export
print.animalUni <- function(x, pars = c("sigma_G", "sigma_E", "beta"), 
                          probs = c(0.025, 0.25, 0.5, 0.75, 0.975), 
                          digits_summary = 2, include = TRUE, ...) { 
  rstan:::print.stanfit(x, pars, probs, digits_summary, include, ...)
}

#'@export
setMethod("summary", signature = "animalUni",
          function(object, pars = c("sigma_G", "sigma_E", "beta"),
                   probs = c(0.025, 0.25, 0.50, 0.75, 0.975), use_cache = TRUE, ...) {
            callNextMethod(object, pars, probs, use_cache, ...)
            })