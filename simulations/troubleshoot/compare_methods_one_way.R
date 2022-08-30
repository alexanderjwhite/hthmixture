library(devtools)
library(dplyr)
library(ggplot2)
load_all()


load("G:/My Drive/Dissertation/HTH Mixture/hthmixture/simulations/real data/Data.RData")

x <- as.matrix(dat.mrm)
y <- as.matrix(dat.vol.std)
N <- x %>% nrow()
k <- 3

N <- x %>% 
  nrow()
p <- x %>% 
  ncol()
m <- y %>% 
  ncol()

set.seed(1)
clust_assign <- fct_initialize(k, N)

# model_1 <- model_k
A_1 <- model_1$Ahat
A_1[A_1!=0] <- 1
NMF::aheatmap(A_1, Rowv = NA, Colv = NA)
which(rowSums(A_1)!=0)

mu_mat1 <- cbind(x,1) %*% model_1$Ahat
sigvec1 <- model_1$sigvec
gam1 <- fct_log_lik(mu_mat1, sigvec1, y, N, m)


# model_2 <- model_k
A_2 <- model_2$Ahat
A_2[A_2!=0] <- 1
NMF::aheatmap(A_2, Rowv = NA, Colv = NA)
which(rowSums(A_2)!=0)

mu_mat2 <- cbind(x,1) %*% model_2$Ahat
sigvec2 <- model_2$sigvec
gam2 <- fct_log_lik(mu_mat2, sigvec2, y, N, m)

# model_3 <- model_k
A_3 <- model_3$Ahat
A_3[A_3!=0] <- 1
NMF::aheatmap(A_3, Rowv = NA, Colv = NA)
which(rowSums(A_3)!=0)

mu_mat3 <- cbind(x,1) %*% model_3$Ahat
sigvec3 <- model_3$sigvec
gam3 <- fct_log_lik(mu_mat3, sigvec3, y, N, m)

gam <- cbind(gam1,gam2,gam3)
apply(gam,1,which.max)
