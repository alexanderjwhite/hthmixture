library(dplyr)
set.seed(19921124)
lam <- 1
maxiter <- 50
N <- 200
prob <- c(0.2, 0.3, 0.5)
k <- prob %>% length()
nvld <- 1e4
rho <- c(0,0,0)
sigma <- 0.01
p <- 3
m <- 3
s <- 1
r <- c(2,2,2)
b <- c(5,10,15)
int <- prob %>% cumsum()
rand_assign <- runif(N)

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
    E <- sim %>% 
      purrr::pluck("E")
    return(list(X=X,Y=Y,A=A,S=S,E=E))
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

clust_assign <- clust_assign_true



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
    
    rhat <- 2
    estimate <- SARRS(Y_k,X_k,rhat, lam, "grLasso")
    
    A_k <- estimate %>% 
      purrr::pluck("Ahat")
    
    mu_mat <- (stack_X %>% 
                 bind_cols(int = rep(1,N)) %>% 
                 as.matrix()) %*% A_k
    
    gam <- array(dim = N)
    for(i in 1:N){
      mu_i <- mu_mat[i,]
      sig_i <- ((stack_Y[i,]-mu_i) %>% matrix()) %*% t((stack_Y[i,]-mu_i) %>% matrix() %>% diag()) %>% abs() %>% sqrt()
      # gam[i, .x] <- 1:m %>% 
      gam[i] <- 1:m %>% 
        purrr::map_dbl(
          .f = function(.y){
            dnorm(stack_Y[i,.y], mean = mu_i[.y], sd = sig_i[.y])
            # dnorm(stack_Y[i,.y], mean = mu_i[.y], sd = 1)
          }
        ) %>%
        log() %>% 
        sum() %>% 
        exp()
    } 
    return(gam)
    
  }) %>% 
  as.matrix()

W <- array(dim=c(N,k))
for(i in 1:k){
  for(j in 1:N){
    W[j,i] <- pi_vec[i]*gamma[j,i]/((pi_vec*gamma[j,]) %>% sum())
  }
}
W[is.na(W)] <- 0

clust_assign_old <- clust_assign
clust_assign <- 1:N %>% 
  purrr::map_int(
    .f = function(.x){
      res <- W[.x,] %>% which.max()
      # if(is.na(res)){res <- which.max(runif(k))}
    }
  )


iter1 <- which(clust_assign_true!=clust_assign)
iter1

case_1 <- 46 #Should be 2, clustered as 1
case_2 <- 49

X_case <- stack_X %>% 
  tibble() %>% 
  slice(case_1) %>% 
  pull() %>% 
  as.matrix()

Y_case <- stack_Y %>% 
  tibble() %>% 
  slice(case_1) %>% 
  pull() %>% 
  as.matrix()

A_1 <- data_k[[1]] %>% 
  purrr::pluck("A")

A_2 <- data_k[[2]] %>% 
  purrr::pluck("A")

Ahat_1 <- A_k

Ahat_2 <- A_k

XA_1 <- X_case %*% A_1

XA_2 <- X_case %*% A_2

sig_1_case <- ((Y_case-XA_1) %>% matrix()) %*% t((Y_case-XA_1) %>% matrix() %>% diag()) %>% abs() %>% sqrt()

sig_2_case <- ((Y_case-XA_2) %>% matrix()) %*% t((Y_case-XA_2) %>% matrix() %>% diag()) %>% abs() %>% sqrt()

XAhat_1 <- (X_case %>% 
              tibble() %>% 
              bind_cols(int = 1) %>% 
              as.matrix()) %*% Ahat_1

XAhat_2 <- (X_case %>% 
              tibble() %>% 
              bind_cols(int = 1) %>% 
              as.matrix()) %*% Ahat_2

sighat_1_case <- ((Y_case-XAhat_1) %>% matrix()) %*% t((Y_case-XAhat_1) %>% matrix() %>% diag()) %>% abs() %>% sqrt()

sighat_2_case <- ((Y_case-XAhat_2) %>% matrix()) %*% t((Y_case-XAhat_2) %>% matrix() %>% diag()) %>% abs() %>% sqrt()

lik_1_i <- rep(0,m)
# for(i in 1:m){lik_1_i[i] <- dnorm(Y_case[i], mean = XA_1[i], sd = sig_1_case[i])}
for(i in 1:m){lik_1_i[i] <- dnorm(Y_case[i], mean = XA_1[i], sd = 1)}
lik_1 <- lik_1_i %>% log() %>% sum() %>% exp()

lik_2_i <- rep(0,m)
# for(i in 1:m){lik_2_i[i] <- dnorm(Y_case[i], mean = XA_2[i], sd = sig_2_case[i])}
for(i in 1:m){lik_2_i[i] <- dnorm(Y_case[i], mean = XA_2[i], sd = 1)}
lik_2 <- lik_2_i %>% log() %>% sum() %>% exp()

likhat_1_i <- rep(0,m)
# for(i in 1:m){likhat_1_i[i] <- dnorm(Y_case[i], mean = XAhat_1[i], sd = sighat_1_case[i])}
for(i in 1:m){likhat_1_i[i] <- dnorm(Y_case[i], mean = XAhat_1[i], sd = 1)}
likhat_1 <- likhat_1_i %>% log() %>% sum() %>% exp()

likhat_2_i <- rep(0,m)
# for(i in 1:m){likhat_2_i[i] <- dnorm(Y_case[i], mean = XAhat_2[i], sd = sighat_2_case[i])}
for(i in 1:m){likhat_2_i[i] <- dnorm(Y_case[i], mean = XAhat_2[i], sd = 1)}
likhat_2 <- likhat_2_i %>% log() %>% sum() %>% exp()
