
N <- 100
k <- 2
rho <- 0.8
sigma <- 1
p <- 200
m <- 250
s_x <- 10
s_y <- 10
rank <- 1
b <- 1
M <- 5
tp <- tn <- fp <- fn <- pos <- neg <- matrix(0, nrow = M, ncol = k)
aris <- rep(0,M)
for(seed in 1:M){
  print(seed)
  set.seed(seed)
  sim <- fct_simulate(N, k, rho, sigma, p, m, s_x, s_y, rank, b, type = "2")

  for(i in 1:k){
    x_k <- sim$x[sim$true==i,]
    y_k <- sim$y[sim$true==i,]
    sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))
    # # lam_univ <- fct_lambda(sigma_hat, ncol(sim$x), nrow(sim$x))
    # lam_univ <- 1
    model <- fct_sarrs(y_k,x_k, 1, lam_univ, sigma = sigma_hat)$Ahat[-(ncol(x_k)+1),]
    op <- which(sim$a[[i]]!=0)
    on <- which(sim$a[[i]]==0)
    ep <- which(model!=0)
    en <- which(model==0)
    tp[seed,i] <- sum(ep %in% op)
    tn[seed,i] <- sum(en %in% on)
    fp[seed,i] <- sum(!(ep %in% op))
    fn[seed,i] <- sum(!(en %in% on))
    pos[seed,i] <- length(op)
    neg[seed,i] <- length(on)
  }


  model <- hthmix(sim$x, sim$y, k, nstart  = 1)
  clust_shift <- clue::solve_LSAP(table(sim$true, model$assign), maximum = TRUE)
  clust_reorder <- clust_shift[model$assign]
  aris[seed] <- mclust::adjustedRandIndex(sim$true, clust_reorder)

}
acc <- (tp+tn)/(pos+neg)
tpr <- tp/pos
tnr <- tn/neg
f1 <- 2*tp/(2*tp+fp+fn)

mean(acc)
mean(tpr)
mean(tnr)
mean(f1)
mean(aris)

N <- 100
k <- 2
rho <- 0.8
sigma <- 1
p <- 200
m <- 250
s_x <- 10
s_y <- 10
rank <- 1
b <- 1
M <- 5

set.seed(1)
sim <- fct_simulate(N, k, rho, sigma, p, m, s_x, s_y, rank, b, type = "2")
model <- hthmix(sim$x, sim$y, k, nstart  = 1)
clust_shift <- clue::solve_LSAP(table(sim$true, model$assign), maximum = TRUE)
clust_reorder <- clust_shift[model$assign]
mclust::adjustedRandIndex(sim$true, clust_reorder)
