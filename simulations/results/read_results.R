sim_res <- readRDS("G:/My Drive/Dissertation/HTH Mixture/hthmixture/simulations/results/20210709_results.rds")

results_summ <- sim_res %>% 
  purrr::map_dfr(.f = function(.x){
    result <- .x %>% purrr::pluck("res")
    
    if(class(result) == "list"){
    
    true <- factor(result$true, levels = 1:(.x$k))
    est <- factor(result$est, levels = 1:(.x$k))
    shuffle <- clue::solve_LSAP(table(true, est), maximum = TRUE)

    acctbl <- (table(true, est)[,shuffle])
    acc <- (acctbl %>% diag() %>% sum()) / .x$N
    return(tibble(id = .x$id, N = .x$N, k = .x$k, sigma = .x$sigma, dim = .x$dim, s = .x$s, r = .x$r, rep = .x$rep, iter = result$iter, time = result$time, acc = acc))
    # return(tibble(id = .x$id, N = .x$N, k = .x$k, sigma = .x$sigma, dim = .x$dim, s = .x$s, r = .x$r, rep = .x$rep, iter = result$iter, time = result$time, acc = "test"))
    } else {
      return(tibble(id = .x$id, N = .x$N, k = .x$k, sigma = .x$sigma, dim = .x$dim, s = .x$s, r = .x$r, rep = .x$rep, iter = NULL, time = NULL, acc = NULL))
    }
    print("hello")
  }) 

results_summ %>% 
  group_by(N,k,sigma,dim,s,r) %>% 
  summarize(iter = mean(iter),
            time=mean(time),
            acc = mean(acc)) %>% 
  tidyr::drop_na() %>% 
  filter(k==2) %>% 
  View()

