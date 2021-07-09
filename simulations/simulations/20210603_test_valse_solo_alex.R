library(dplyr)
library(ggplot2)
library(valse)
source("./functions/20210205_sarrs_alex.R")
set.seed(19921124)
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

gamma_store <<- as_tibble(matrix(rep(0,k),nrow=1)) %>% 
  rename_if(is.numeric,~names) %>% 
  mutate(w_ll = 0, iter = 0)

# chain_store <- tibble(obs = 0, assign = 0, llik = 0, chain = 0)
tictoc::tic()
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

runValse(stack_X, stack_Y, procedure = "LassoRank", fast = FALSE, verbose = TRUE, kmin = 2, kmax = 2)
