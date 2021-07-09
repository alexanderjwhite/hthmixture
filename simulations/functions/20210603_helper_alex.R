boltzmann_pert <- function(lik){
  if(is.null(lik)|length(lik)==1){
    return(1)
  }
  # const <- 100000
  # lik <- c(1,1.5,10,10)*const
  rel <- ((lik %>% lead())/lik) %>% .[!is.na(.)]
  # kt <- norm(lik, type = "2")
  # kt <- median(lik)
  kt <- 1
  
  pert <- 1-exp(rel[length(rel)]/kt)/sum(exp(rel/kt))
  # print(paste(lik,"pert:",pert))
  return(pert)
  
}
