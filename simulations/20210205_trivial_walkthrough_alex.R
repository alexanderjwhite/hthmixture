# Regular Multivariate regression
tictoc::tic()
p <- 5000
m <- 5000
s <- 30
a_rows <- 1:s
r <- 2
b <- 10
rho <- 0
n <- 200
lam <- 0.5

A <- matrix(0,p,m)
nrow <- length(a_rows)
B0 <- matrix(rnorm(nrow*r),nrow,r)
B1 <- matrix(rnorm(r*m),r,m)
tictoc::tic()
A[a_rows,] <- b*B0%*%B1
tictoc::toc()

Sigma <- matrix(1,p,p)
for(j in 1:p) for (k in 1:p) Sigma[j,k] <- rho^abs(j-k)

tictoc::tic()
X <- MASS::mvrnorm(n,rep(0,p),Sigma)	
tictoc::toc()
tictoc::tic()
E <- matrix(rnorm(n*m),n,m)
tictoc::toc()

tictoc::tic()
Y <- X%*%A + E
tictoc::toc()

# A_reg <- array(0,dim = c(p,m))
# for(i in 1:m){
#   model <- lm(Y[,i]~X[,1:s]) 
#   A_reg[1:(s+1),i] <- model$coefficients
# }
# 
# par(mfrow=c(1,2))
# NMF::aheatmap(A, Rowv = NA, Colv = NA)
# NMF::aheatmap(A_reg[2:p,], Rowv = NA, Colv = NA)

# SARRS

rhat <- 2
lam <- 1
tictoc::tic()
estimate <- SARRS(Y,X,rhat, lam, "grLasso")
tictoc::toc()
A_sarrs <- estimate %>% 
  purrr::pluck("Ahat")

tictoc::tic()
((A_sarrs[1:p,]-A) %>% norm(type = "2"))/(A %>% norm(type="2"))
# A %>% norm(type="2")
tictoc::toc()
tictoc::toc()

par(mfrow=c(1,2))
NMF::aheatmap(A[1:200,1:200], Rowv = NA, Colv = NA)
NMF::aheatmap(A_sarrs[1:200,1:200], Rowv = NA, Colv = NA)

# Single iteration

lam <- 1
maxiter <- 50
N <- 200
prob <- c(0.2, 0.3, 0.5)
k <- prob %>% length()
nvld <- 1e4
rho <- c(0,0,0)
sigma <- 0.01
p <- 50
m <- 20
s <- 15
r <- c(2,2,2)
b <- c(5e3,10e3,15e3)
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
iter2 <- which(clust_assign_true!=clust_assign)

true_A <- data_k[[3]]$A
par(mfrow=c(1,2))
NMF::aheatmap(true_A, Rowv = NA, Colv = NA)
NMF::aheatmap(A_k, Rowv = NA, Colv = NA)

