library(dplyr)
library(ggplot2)
source("./functions/20210205_sarrs_alex.R")
# set.seed(19921124)
lam <- 1
maxiter <- 100
N <- 100
prob <- c(0.5,0.5)
k <- prob %>% length()
nvld <- 1e4
rho <- c(0,0)
sigma <- 1
p <- 50
m <- 20
s <- 5
r <- c(1,1)
b <- c(5,10)
int <- prob %>% cumsum()
rand_assign <- runif(N)
names <- paste0("c_",1:k)
if(k*s > p){print("FOR SEPARATION, VERIFY THAT K*S < P")}

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

gamma_store <- as_tibble(matrix(rep(0,k),nrow=1)) %>% 
  rename_if(is.numeric,~names) %>% 
  mutate(w_ll = 0, iter = 0)
  

  
# clust_assign_true <- c(rep(1,39),rep(2,52), rep(3,109))
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
    
    sim <- sim_sarrs(.n,nvld,.n,a_rows,p,m,.r,.rho,sigma,.b)
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

stack_X <- data_k %>% 
  purrr::map_dfr(.f = function(.x){
    X <- .x %>% 
      purrr::pluck("X")
    return(as_tibble(X))
  }) %>% 
  slice(clust_assign_true_vec) %>% 
  as.matrix()

stack_Y <- data_k %>% 
  purrr::map_dfr(.f = function(.x){
    Y <- .x %>% 
      purrr::pluck("Y")
    return(as_tibble(Y))
  }) %>% 
  slice(clust_assign_true_vec) %>% 
  as.matrix()


init_int <- rep(1/k,k) %>% cumsum()
init_rand_assign <- runif(N)

clust_assign <- (init_rand_assign) %>% 
  purrr::map_int(.f = function(.x){
    clust <- (.x <= init_int) %>% 
      which() %>% 
      min()
    return(clust)
  })

# clust_assign <- clust_assign_true

conv <- Inf
main_clust <- Inf
iter <- 0
while(conv>0 & iter < maxiter){
# while(iter < maxiter){
  iter <- iter+1
  
  pi_vec <- 1:k %>% 
    purrr::map_dbl(
      .f = function(.x){
        sum(clust_assign==.x)/N
      }
    )
  
  
  gamma <- 1:k %>% 
    list() %>% 
    purrr::pmap_dfc(.f = function(.x){
      
      select <- (clust_assign==.x) %>% 
        which()
      
      nk <- select %>% length()
      
      X_k <- stack_X %>% 
        as_tibble() %>% 
        slice(select) %>% 
        as.matrix()
      
      Y_k <- stack_Y %>% 
        as_tibble() %>% 
        slice(select) %>% 
        as.matrix()
      
      # P=X%*% ginv(t(X)%*%X) %*%t(X)
      # rhat = sum(svd(P%*%Y)$d > sigmahat*(sqrt(2*m)+sqrt(2*min(n,p))))
      
      rhat <- r[.x]
      estimate <- SARRS(Y_k,X_k,rhat, lam, "grLasso")
      
      A_k <- estimate %>% 
        purrr::pluck("Ahat")
      
      sig_vec <- estimate %>% 
        purrr::pluck("sigvec")
      
      mu_mat <- (stack_X %>% 
                   bind_cols(int = rep(1,N)) %>% 
                   as.matrix()) %*% A_k
      
      gam <- array(dim = N)
      for(i in 1:N){
        mu_i <- mu_mat[i,]
        # sig_i <- ((stack_Y[i,]-mu_i) %>% matrix()) %*% t((stack_Y[i,]-mu_i) %>% matrix() %>% diag()) %>% abs() %>% sqrt()
        # gam[i, .x] <- 1:m %>% 
        gam[i] <- 1:m %>% 
          purrr::map_dbl(
            .f = function(.y){
              # dnorm(stack_Y[i,.y], mean = mu_i[.y], sd = sig_i[.y])
              dnorm(stack_Y[i,.y], mean = mu_i[.y], sd = sig_vec[.y])
            }
          ) %>%
          log() %>% 
          sum() #%>% 
          # exp()
      } 
      # gam[is.infinite(gam)] <- min(gam[is.finite(gam)], na.rm = TRUE)
      return(tibble(gam) %>% setNames(names[.x]))
      
    }) 
  
  weighted_ll <- (gamma*pi_vec) %>% 
    tibble() %>% 
    mutate(total = rowSums(across(where(is.numeric)))) %>% 
    pull(total)
  
  gamma_store <- gamma_store %>% 
    bind_rows(tibble(gamma, w_ll = weighted_ll, iter = iter))
  
  # gamma_trans <- 1:N %>% 
  #   purrr::map_dfr(
  #     .f = function(.z){
  #       obs <- gamma %>% 
  #         slice(.z) %>%
  #         as.matrix(dimnames = NULL) %>% 
  #         as.numeric()
  #       
  #       obs[is.infinite(obs)] <- NA
  #       obs <- obs-median(obs,na.rm=TRUE)
  #       obs[is.na(obs)] <- -Inf
  #       obs <- obs-max(obs)
  #       obs <- as.matrix(obs) %>% t()
  #       return(as_tibble(obs) %>% setNames(names))
  #     }
  #   ) %>% 
  #   as.matrix()
  # 
  # %>% 
  #   rowwise() %>% 
  #   mutate(med = median(c_across(1:k), na.rm = TRUE)) %>% 
  #   mutate(med = function(.x){median(c_across(.x))})
  #   mutate_all(~(.x-med)) %>% 
  #   select(-med) %>% 
  #   as.matrix()
  # 

  # 
  # gamma %>% 
  #   rowwise() %>% 
  #   mutate(med = median(c_across(1:k))) %>% 
  #   mutate_all(~(.x-med)) %>% 
  #   select(-med)
  # 
  # gamma_trans <- list(gamma[,1], gamma[,2]) %>% 
  #   purrr::pmap_dfr(
  #     .f = function(.x,.y){
  #       res <- c(.x,.y)-mean(c(.x,.y))
  #       return(tibble(val1 = res[1], val2 = res[2]))
  #     }
  #   ) %>%
  #   as.matrix()
  
  # W <- array(dim=c(N,k))
  # for(i in 1:k){
  #   for(j in 1:N){
  #     #W[j,i] <- pi_vec[i]*gamma[j,i]/((pi_vec*gamma[j,]) %>% sum())
  #     part_1 <- log(pi_vec[i])+gamma_trans[j,i]
  #     part_2 <- (log(pi_vec)+gamma_trans[j,]) %>% exp() %>% sum() %>% log()
  #     W[j,i] <- exp(part_1 - part_2)
  #   }
  # }
  # W[is.na(W)] <- 0
  
  clust_assign_old <- clust_assign
  clust_assign <- 1:N %>% 
    purrr::map_int(
      .f = function(.x){
       # res <- W[.x,] %>% which.max()
        res <- (gamma %>% as.matrix())[.x,] %>% which.max()
        # if(is.na(res)){res <- which.max(runif(k))}
      }
    )
  
  # shuffle <- clue::solve_LSAP(table(clust_assign, clust_assign_true), maximum = TRUE)
  # clust_assign_temp <- clust_assign*10
  # for (i in 1:k){
  #   clust_assign_temp[clust_assign_temp==(shuffle[i]*10)]
  # }
  
  conv <- (clust_assign!=clust_assign_old) %>% sum()
  #shuffle <- clue::solve_LSAP(table(clust_assign_true, clust_assign), maximum = TRUE)
  print((table(clust_assign_true, clust_assign)))
  #print(clust_assign_old)
  #print(clust_assign)
  
  
  # 
  # shuffle <- clue::solve_LSAP(table(clust_assign_true, clust_assign), maximum = TRUE)
  # 
  # 
  # main_clust_old <- main_clust
  # main_clust <- (table(clust_assign_true, clust_assign)[,shuffle]) %>% diag() %>% sum()
  # conv <- (main_clust==main_clust_old)
  print(paste("Iter:",iter,"| conv:",conv))
  # print(main_clust)
}

shuffle <- clue::solve_LSAP(table(clust_assign_true, clust_assign), maximum = TRUE)
(table(clust_assign_true, clust_assign)[,shuffle])

p1 <- gamma_store %>% 
  filter(iter > 0) %>% 
  group_by(iter) %>% 
  summarize_all(~sum(.)) %>% 
  tidyr::pivot_longer(cols = tidyselect::starts_with("c"), names_to = "cluster") %>% 
  ggplot() +
  geom_path(aes(x = iter, y = value, colour = cluster))
# plotly::ggplotly(p1)
p1

p2 <- gamma_store %>% 
  filter(iter > 0) %>% 
  group_by(iter) %>% 
  summarize_all(~sum(.)) %>% 
  ggplot() +
  geom_path(aes(x = iter, y = w_ll)) 
# plotly::ggplotly(p2)
p2
