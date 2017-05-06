#' Plot a simulated model
#'
#' Basic plotting functionality for an object of class \code{azrsim}, generated using
#' the simulation function \code{AZRsim::simulate.azrmod()}.
#'
#' @param object An object of class \code{azrsim}
#' @param ... Additional arguments to \code{plot()}.
#' @return A lattice plot of the variable simulations over time.
#' @examples
#' sho <- create_model(system.file("examples/sho.txt", package="AZRsim"))
#' sho_sim <- simulate(sho, 100)
#' plot(sho_sim, lwd = 2)
#' @export
plot.azrsim <- function(object, mfrow_val = NULL, pars = NULL, ...) {
  call <- match.call()
  if (!("azrsim" %in% class(object))) {
    stop(paste0(call$object, " must be of class 'azrsim'."))
  }
  param_names <- names(object)[-which(names(object) == "TIME")]
  if(is.null(mfrow_val) & is.null(pars)) {
    param_count <- length(param_names)
    num_plots <- c(1:min(param_count, 9))

    if(param_count %in% c(1,2,3)) {
      par(mfrow = c(1, param_count))
      for(name in param_names[num_plots]) {
        plot(object[["TIME"]], object[[name]], main = name, ylab = "", xlab = "Time", type = "l", ...)
      }
    }
    else if(param_count %in% c(4,5,6)) {
      par(mfrow = c(2, param_count))
      for(name in param_names[num_plots]) {
        plot(object[["TIME"]], object[[name]], main = name, ylab = "", xlab = "Time", type = "l", ...)
      }
    }
    else if(param_count %in% c(7,8,9)) {
      par(mfrow = c(3, 3))
      for(name in param_names[num_plots]) {
        plot(object[["TIME"]], object[[name]], main = name, ylab = "", xlab = "Time", type = "l", ...)
      }
    }
    else {
      warning(paste0("Plotting 9 out of ", param_count, " variables"), call. = FALSE)
      par(mfrow = c(3, 3))
      for(name in param_names[num_plots]) {
        plot(object[["TIME"]], object[[name]], main = name, ylab = "", xlab = "Time", type = "l", ...)
      }
    }
  }
  else {
    if(!is.null(mfrow_val) & is.null(pars)) {
      par(mfrow = mfrow_val)
      for(name in param_names[1:prod(mfrow_val)]) {
        plot(object[["TIME"]], object[[name]], main = name, ylab = "", xlab = "Time", type = "l", ...)
      }
    }
    else if(!is.null(mfrow_val) & !is.null(pars)) {
      if(!(prod(mfrow_val) == length(pars))) {
        stop("Mismatch between number of plots and number of parameters declared.")
      }
      par(mfrow = mfrow_val)
      for(name in pars) {
        plot(object[["TIME"]], object[[name]], main = name, ylab = "", xlab = "Time", type = "l", ...)
      }
    }
    else if (is.null(mfrow_val) & !is.null(pars)) {
      if(length(pars) <= 3) {
        par(mfrow = c(1,length(pars)))
        for(name in pars) {
          plot(object[["TIME"]], object[[name]], main = name, ylab = "", xlab = "Time", type = "l", ...)
        }
      }
      else if(length(pars) <= 6) {
        par(mfrow = c(1,length(pars)))
        for(name in pars) {
          plot(object[["TIME"]], object[[name]], main = name, ylab = "", xlab = "Time", type = "l", ...)
        }
      }
      else if(length(pars) <= 9) {
        par(mfrow = c(1,length(pars)))
        for(name in pars) {
          plot(object[["TIME"]], object[[name]], main = name, ylab = "", xlab = "Time", type = "l", ...)
        }
      }
    }
  }
}
