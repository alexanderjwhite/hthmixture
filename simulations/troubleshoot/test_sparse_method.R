N <- 100
k <- 1
rho <- 0.8
sigma <- 1
p <- 500
m <- 250
s_x <- 10
s_y <- 10
rank <- 1
b <- 1
M <- 1
tp <- tn <- fp <- fn <- pos <- neg <- cv_error <- rep(0,M)
set.seed(1)
for(i in 1:M){
  print(i)
  sim <- fct_simulate(N, k, rho, sigma, p, m, s_x, s_y, rank, b, type = "2")
  sigma_hat <- fct_sigma(sim$y, nrow(sim$y), ncol(sim$y))
  
  model <- fct_sarrs(sim$y, sim$x, 1, median(store), sigma = sigma_hat)$Ahat[-(ncol(sim$x)+1),]
  op <- which(sim$a[[1]]!=0)
  on <- which(sim$a[[1]]==0)
  ep <- which(model!=0)
  en <- which(model==0)
  tp[i] <- sum(ep %in% op)
  tn[i] <- sum(en %in% on)
  fp[i] <- sum(!(ep %in% op))
  fn[i] <- sum(!(en %in% on))
  pos[i] <- length(op)
  neg[i] <- length(on)

}


acc <- (tp+tn)/(pos+neg)
tpr <- tp/pos
tnr <- tn/neg
f1 <- 2*tp/(2*tp+fp+fn)

mean(acc)
mean(tpr)
mean(tnr)
mean(f1)

set.seed(1)
sim <- fct_simulate(N, k, rho, sigma, p, m, s_x, s_y, rank, b, type = "2")
store <- NULL
for(j in 1:100){
  print(j)
  clust_assign <- fct_initialize(k,N)
  for(i in 1:k){
    x_k <- sim$x[which(clust_assign==i),]
    y_k <- sim$y[which(clust_assign==i),]
    sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))
    store <- c(store,fct_sarrs(y_k, x_k, 1, NULL, sigma = sigma_hat)$lambda_store)
  }
}
min(store)
