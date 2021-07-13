result %>% 
  purrr::map_dfr(.f = function(.x){
    result <- .x %>% purrr::pluck("res")
    true <- result$true
    est <- result$est
    shuffle <- clue::solve_LSAP(table(true, est), maximum = TRUE)
    
    acctbl <- (table(true, est)[,shuffle])
    acc <- (acctbl %>% diag() %>% sum()) / .x$N
    return(tibble(id = .x$id, N = .x$N, k = .x$k, sigma = .x$sigma, dim = .x$dim, s = .x$s, r = .x$r, rep = .x$rep, iter = result$iter, time = result$time, acc = acc))
    
  })
