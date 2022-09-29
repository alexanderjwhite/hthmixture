N <- 100
k <- 2
sigma <- 1
p <- 200
m <- 250
rank <- 1
b <- 1
d <- 20
h <- 0.2
case <- "independent"

set.seed(1)
sim <- fct_simulate(N, k, sigma, p, m, rank, b, d, h, case)

model <- hthmix(sim$x, sim$y, k, nstart = 1, init_assign = NULL, selection = "universal", 
                alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, rank = NULL, max_rank = 3, 
                maxiter = 1e3, verbose = TRUE, temp = 1000, p1 = 0.9, p2 = 0.95, 
                sim_N = 200, checks = 3, true = sim$true)

# model <- hthmix(sim$x, sim$y, k, selection = "universal", maxiter = 10)

clust <- model$assign
true <- sim$true
clust_reorder <- clue::solve_LSAP(table(true, clust), maximum = TRUE)[clust]

ari <- mclust::adjustedRandIndex(true, clust_reorder)
ari

# fct_j_lik(sim$x, sim$y, k, sim$true)
# fct_j_lik(sim$x, sim$y, k, model$assign)
# fct_j_lik(sim$x, sim$y, k, fct_initialize(k,N))

library(ggplot2)
p1 <- model$ll_store %>% 
  ggplot(aes(x = iter, y = ll, color = type )) +
  geom_point()
# 30851.8341042094


results <- read.table("C:/Users/whitealj/OneDrive - Indiana University/FromGoogle/Dissertation/HTH Mixture/iu_2022_09_25_01/summary.csv")
colnames(results) <- c("N", "k", "sigma", "p", "m", "seed", "rank", "b", "ari", "time")
results %>% mutate(d = 20) %>% View()
results %>% View()

results %>% 
  group_by(N,k,sigma,p,m,rank,b) %>% 
  summarize(n = n(),
            ari_mean = mean(ari),
            ari_sd = sd(ari),
            max_time = max(time)) %>% 
  View()


model <- readRDS("C:/Users/whitealj/OneDrive - Indiana University/FromGoogle/Dissertation/HTH Mixture/iu_2022_09_25_01//output/iu_hthmix_100_3_0.1_dependent_200_250_4_4_1_3.rds")
model$models[[1]]$ll_store %>% 
  ggplot(aes(x = iter, y = ll, color = type )) +
  geom_point()













N <- 100
k <- 3
sigma <- 1
p <- 200
m <- 250
rank <- 1
b <- 1
d <- 20
h <- 0.2
case <- "independent"
max_rank <- 3
aris <- rep(0,10)
# for (seed in 2001:2010){
  # print(seed)
  set.seed(1)
  sim <- fct_simulate(N, k, sigma, p, m, rank, b, d, h, case)
  
  model <- hthmix(sim$x, sim$y, k, nstart = 1, init_assign = NULL, selection = "universal", 
                  alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, rank = NULL, max_rank = 3, 
                  maxiter = 1e3, verbose = TRUE, temp = 1000, p1 = 0.9, p2 = 0.95, 
                  sim_N = 200, checks = 3, true = sim$true)
  
  # model <- hthmix(sim$x, sim$y, k, selection = "universal", maxiter = 10)
  
  clust <- model$assign
  true <- sim$true
  clust_reorder <- clue::solve_LSAP(table(true, clust), maximum = TRUE)[clust]
  
  aris <- mclust::adjustedRandIndex(true, clust_reorder)
# } 0.5652265
aris

######################################
model$llik
fct_j_lik(sim$x,sim$y,k,sim$true,lambda = c(8.63996131521922,2.5211124681541,8.63996131521922))

ind <- which(model$assign==1)
x_k <- sim$x[ind,]
y_k <- sim$y[ind,]

sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))


safe_rank <- purrr::safely(fct_rank)
eta_k <- 3
rank_sest <- safe_rank(x_k, y_k, sigma_hat, eta_k)
rank_hat <- ifelse(is.null(rank_sest$result),1,rank_sest$result)
rank_hat <- min(rank_hat, max_rank)
fct_sarrs(y_k, x_k, r = 3, lam = NULL, 
          alpha = 2*sqrt(3), beta = 1, sigma = sigma_hat,
          ptype = "grLasso", y_sparse = TRUE)$lambda_store










#############################
N <- 100
k <- 3
sigma <- 1
p <- 200
m <- 250
rank <- 1
b <- 1
d <- 20
h <- 0.2
case <- "independent"
set.seed(1)
sim <- fct_simulate(N, k, sigma, p, m, rank, b, d, h, case)

model <- hthmix(x = sim$x, y = sim$y, k = k, nstart = 1, init_assign = NULL, init_lambda = lambda, alt_iter = 3, anneal_iter = 1e3, 
                em_iter = 1e3, temp = 1e3, mu = 0.95, eps = 1e-6, accept_prop = 0.95, sim_N = 200)



clust <- model$assign
true <- sim$true
clust_reorder <- clue::solve_LSAP(table(true, clust), maximum = TRUE)[clust]

mclust::adjustedRandIndex(true, clust_reorder)
# 0.5652265
