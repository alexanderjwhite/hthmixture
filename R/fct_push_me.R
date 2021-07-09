#' Pushover Capability for Large Simulation
#'
#' @param message a message about script progress
#'
#' @return nothing
#' @export
#'
#' @import pushoverr
#' @examples fct_push_me("test")
fct_push_me <- function(message = "complete"){
  pushoverr::pushover(message = message,
                      user = "umdyjtpwtyn4eg41b1zztc9kp4197o",
                      app = "aihd2o1bzds63ezxz8567dqv9fsv22")
}