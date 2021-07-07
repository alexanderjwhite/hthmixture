library(dplyr)
library(ggplot2)
# library(valse)
# source("./functions/20210205_sarrs_alex.R")
# set.seed(19921124)
num_pert <- 20
chains <- 100
lam <- 1
maxiter <- 100
N <- 100
prob <- c(0.5,0.5)
k <- prob %>% length()
nvld <- 1e4
rho <- c(0,0)
sigma <- 1
p <- 200
m <- 200
s <- 10
r <- c(1,1)
b <- c(5,10)
int <- prob %>% cumsum()
rand_assign <- runif(N)
names <- paste0("c_",1:k)
# if(k*s > p){print("FOR SEPARATION, VERIFY THAT K*S < P")}

# make sure k*s is less than p

clust_assign_true <- (rand_assign) %>% 
  purrr::map_int(.f = function(.x){
    clust <- (.x <= int) %>% 
      which() %>% 
      min()
    return(clust)
  }) %>% 
  sort()

clust_assign_true_key <- clust_assign_true %>% 
  tibble() %>% 
  mutate(order = 1:N) %>% 
  arrange((.)) 

clust_assign_true_vec <- clust_assign_true_key %>% 
  pull(order)


n <- clust_assign_true %>% 
  as_tibble() %>% 
  group_by(value) %>% 
  summarize(n = n()) %>% 
  pull(n)

set.seed(19921124)
models_with_lam <- list()
models_cv_lam <- list()
models_both <- list()
for (i in 1:1){
  clust_iter <<- 1
  clust_min <<- 1
  clust_max <<- s
  data_k <- n %>% 
    list(r,rho,b) %>% 
    purrr::pmap(.f = function(.n,.r,.rho,.b){
      a_rows <- clust_min:clust_max
      clust_iter <<- clust_iter + 1
      clust_min <<- clust_max+1
      clust_max <<- clust_iter*s
      
      sim <- fct_sim_mixrrr(.n,a_rows,p,m,.r,.rho,sigma,.b)
      A <- sim %>% 
        purrr::pluck("A")
      S <- sim %>% 
        purrr::pluck("sig")
      X <- sim %>% 
        purrr::pluck("X")
      Y <- sim %>% 
        purrr::pluck("Y")
      return(list(X=X,Y=Y,A=A,S=S))
    })
  
  x <- data_k %>% 
    purrr::map_dfr(.f = function(.x){
      X <- .x %>% 
        purrr::pluck("X")
      return(as_tibble(X))
    }) %>% 
    slice(clust_assign_true_vec) %>% 
    as.matrix()
  
  y <- data_k %>% 
    purrr::map_dfr(.f = function(.x){
      Y <- .x %>% 
        purrr::pluck("Y")
      return(as_tibble(Y))
    }) %>% 
    slice(clust_assign_true_vec) %>% 
    as.matrix()
  
  n_iter <- 436
  model <- hthmix(x, y, chains = 1, maxiter = 7)
  final_assign <- model %>% purrr::pluck("result","assign","final_assign")
  
  model %>% 
    purrr::pluck("result","assign","assign_store") %>% 
    mutate(obs = rep(1:N, (n_iter+1)),
           truth = rep(clust_assign_true, (n_iter+1)),
           correct = assign==truth) %>% 
    ggplot() +
    geom_raster(aes(x=iter,y=obs, fill = correct))
  
  model %>% 
    purrr::pluck("result","assign","lik_store") %>% 
    ggplot() +
    geom_line(aes(x=iter,y=ll))
  
  shuffle <- clue::solve_LSAP(table(clust_assign_true, final_assign), maximum = TRUE)
  acctbl <- (table(clust_assign_true, final_assign)[,shuffle])
  acc <- (acctbl %>% diag() %>% sum()) / N
  acc
  models_with_lam <- append(models_with_lam,list(hthmix(x, y, chains = 5, lam = 1)))
  models_cv_lam <- append(models_cv_lam, list(hthmix(x, y, chains = 1)))
  # models_both <- append(models_both, list(hthmix(x, y, chains = 10)))
  
  best_model <- model %>%
    arrange(-llik) %>%
    slice(1) %>%
    pull(assign) %>%
    purrr::pluck(1)

  shuffle <- clue::solve_LSAP(table(clust_assign_true, best_model), maximum = TRUE)

  acctbl <- (table(clust_assign_true, best_model)[,shuffle])
  acc <- (acctbl %>% diag() %>% sum()) / N
  acc
  print("##############")
  print(i)
  print("##############")
}

# models_with_lam %>% readr::write_rds("./results/no_cv.rds")
# models_cv_lam %>% readr::write_rds("./results/cv.rds")
# models_both %>% readr::write_rds("./results/both.rds")

1 %>% 
  purrr::map_dfr(.f = function(.i){
    # time_1 <- models_with_lam %>% 
    #   purrr::pluck(.i,"time")
    
    time_2 <- models_cv_lam %>% 
      purrr::pluck(.i,"time")
    
    # time_3 <- models_both %>% 
      # purrr::pluck(.i,"time")
    # 
    # best_model_1 <- models_with_lam %>% 
    #   purrr::pluck(.i, "result") %>% 
    #     arrange(-llik) %>%
    #     slice(1) %>%
    #     pull(assign) %>%
    #     purrr::pluck(1)
    # 
    #   shuffle_1 <- clue::solve_LSAP(table(clust_assign_true, best_model_1), maximum = TRUE)
    # 
    #   acctbl_1 <- (table(clust_assign_true, best_model_1)[,shuffle_1])
    #   acc_1 <- (acctbl_1 %>% diag() %>% sum()) / N
      
      best_model_2 <- models_cv_lam %>% 
        purrr::pluck(.i, "result") %>% 
        arrange(-llik) %>%
        slice(1) %>%
        pull(assign) %>%
        purrr::pluck(1)
      
      shuffle_2 <- clue::solve_LSAP(table(clust_assign_true, best_model_2), maximum = TRUE)
      
      acctbl_2 <- (table(clust_assign_true, best_model_2)[,shuffle_2])
      acc_2 <- (acctbl_2 %>% diag() %>% sum()) / N
      
      # best_model_3 <- models_both %>% 
      #   purrr::pluck(.i, "result") %>% 
      #   arrange(-llik) %>%
      #   slice(1) %>%
      #   pull(assign) %>%
      #   purrr::pluck(1)
      # 
      # shuffle_3 <- clue::solve_LSAP(table(clust_assign_true, best_model_3), maximum = TRUE)
      # 
      # acctbl_3 <- (table(clust_assign_true, best_model_3)[,shuffle_3])
      # acc_3 <- (acctbl_3 %>% diag() %>% sum()) / N
    
    
    return(tibble(iter = .i, 
                  time_lam = time_1, 
                  time_cv = time_2,
                  # time_both = time_3,
                  acc_lam = acc_1,
                  acc_cv = acc_2#,
                  # acc_both = acc_3
                  ))
  }) %>% 
  summarize_if(is.numeric,mean)
































