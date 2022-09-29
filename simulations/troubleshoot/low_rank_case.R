library(devtools)
library(ggplot2)
library(gridExtra)
load_all()

# Set up simulation
N <- 100
k <- 2
rho <- 0.9
sigma <- 1
p <- 200
m <- 250
s_x <- 10
s_y <- 10
rank <- 1
b <- 0.5

for(seed in 1:5){
  # reproduce
  set.seed(seed)
  sim <- fct_simulate(N, k, rho, sigma, p, m, s_x, s_y, rank, b, type = "2")
  
  # setup
  
  x <- sim$x
  y <- sim$y
  N <- x %>% nrow()
  
  assignments <- NULL
  likelihood <- NULL
  A <- NULL
  
  clust_assign <- fct_initialize(k, N)
  new_change_ll <- -Inf
  skip_check <- FALSE
  old_ll <- -Inf
  clust_reorder <- clue::solve_LSAP(table(sim$true, clust_assign), maximum = TRUE)[clust_assign]
  clust_data <- tibble(iter = factor("0", levels = c(0:100,"true")), row_id = 1:N, assign = as.factor(clust_reorder), ari = 0)
  
  N <- nrow(x)
  p <- ncol(x)
  m <- ncol(y)
  
  iter <- 0
  conv <- Inf
  clust_store <- tibble(iter=rep(iter,N),assign=clust_assign)
  ll_store <- tibble(iter = 0, ll = -Inf)
  
  # open pdf
  pdf(file = paste0("./simulations/troubleshoot/low_rank_case_method_",seed,".pdf"), width = 12, height = 12)
  while(conv > 0){
    # iteration 1
    iter <- iter + 1
    # print(paste("Iteration",i,"..."))
    pi_vec <- fct_pi_vec(clust_assign, k, N)
    
    gamma_model <- fct_gamma(x, y, k, N, clust_assign, selection = "universal", alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, max_rank = 3, rank = NULL)
    gamma <- gamma_model$gamma
    # print(gamma)
    A <- gamma_model$A
    sig_vec <- gamma_model$sig_vec
    weighted_ll <- fct_weighted_ll(gamma)
    if(weighted_ll > old_ll){
      ll_store <- ll_store %>% bind_rows(tibble(iter = iter, ll = weighted_ll))
      clust_assign_old <- clust_assign
      clust_assign <- fct_update_clust(gamma, N)
    } else {
      clust_assign_old <- clust_assign
      skip_check <- TRUE
    }
    old_ll <- weighted_ll
    clust_shift <- clue::solve_LSAP(table(sim$true, clust_assign), maximum = TRUE)
    clust_reorder <- clust_shift[clust_assign]
    ari <- mclust::adjustedRandIndex(sim$true, clust_reorder)
    clust_data <- rbind(clust_data, tibble(iter = paste0(iter), row_id = 1:N, assign = as.factor(clust_reorder), ari = ari))
    
    # evolution of assignment
    p1 <- clust_data %>%
      rbind(tibble(iter = "true", row_id = 1:N, assign = as.factor(sim$true), ari = 1)) %>%
      # mutate(iter = factor(iter, levels = c(0:iter,"true"))) %>%
      ggplot(aes(x = iter, y = row_id, fill = assign)) +
      geom_tile() +
      ggtitle(paste("Iteration:",iter,"| ARI:",round(ari,3)))
    
    p2 <- clust_data %>%
      dplyr::select(iter, ari) %>%
      distinct() %>%
      ungroup() %>%
      mutate(iter = as.numeric(iter)) %>%
      ggplot(aes(x = iter, y = ari)) +
      geom_line() +
      geom_point() +
      ggtitle("ARI")
    
    # evolution of likelihood
    p3 <- ll_store %>%
      ggplot(aes(x = iter, y = ll)) +
      geom_line() +
      geom_point() +
      ggtitle("Likelihood")
    
    # Accuracy of A1
    A1 <- A[[clust_shift[1]]]
    A1[A1 != 0] <- 1
    p4 <- as_tibble(A1) %>%
      tidyr::pivot_longer(cols = tidyselect::everything(),names_to = "col", values_to = "val") %>%
      mutate(row = rep(1:(nrow(A1)), each = ncol(A1))) %>%
      mutate(val = as.factor(val)) %>%
      ggplot(aes(x = col, y = row, fill = val)) +
      geom_tile() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ggtitle("A1hat")
    
    A1_true <- sim$a[[1]]
    A1_true[A1_true != 0] <- 1
    p5 <- as_tibble(A1_true) %>%
      tidyr::pivot_longer(cols = tidyselect::everything(),names_to = "col", values_to = "val") %>%
      mutate(row = rep(1:(nrow(A1_true)), each = ncol(A1_true))) %>%
      mutate(val = as.factor(val)) %>%
      ggplot(aes(x = col, y = row, fill = val)) +
      geom_tile() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ggtitle("A1 True")
    
    # Accuracy of A2
    A2 <- A[[clust_shift[2]]]
    A2[A2 != 0] <- 1
    p6 <- as_tibble(A2) %>%
      tidyr::pivot_longer(cols = tidyselect::everything(),names_to = "col", values_to = "val") %>%
      mutate(row = rep(1:(nrow(A2)), each = ncol(A2))) %>%
      mutate(val = as.factor(val)) %>%
      ggplot(aes(x = col, y = row, fill = val)) +
      geom_tile() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ggtitle("A2hat")
    
    A2_true <- sim$a[[2]]
    A2_true[A2_true != 0] <- 1
    p7 <- as_tibble(A2_true) %>%
      tidyr::pivot_longer(cols = tidyselect::everything(),names_to = "col", values_to = "val") %>%
      mutate(row = rep(1:(nrow(A2_true)), each = ncol(A2_true))) %>%
      mutate(val = as.factor(val)) %>%
      ggplot(aes(x = col, y = row, fill = val)) +
      geom_tile() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ggtitle("A2 True")
    
    plots <- list(p1,p2,p3,p4,p5,p6,p7)
    lay <- rbind(c(1,1,2,2,3,3),c(4,4,4,5,5,5),c(6,6,6,7,7,7))
    do.call(grid.arrange, list(grobs = plots, layout_matrix = lay))
    
    conv <- (clust_assign != clust_assign_old) %>% sum()
    print(paste("i: ",iter, "| conv: ", conv))
    
    # if(conv==0 & !skip_check){
    #
    #   print("checking for final convergence...")
    #   clust_assign_old <- clust_assign
    #   # new_change <- fct_conv_check2(x, y, k, N, clust_assign, selection = "universal", alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, max_rank = 3, rank = NULL,ll = weighted_ll)
    #   clust_assing <- fct_sim_anneal(sim$x, sim$y, k, clust_assign, t_1 = 1000, mu = 0.9, eps = 1e-6, p = 0.9, N = 200)
    #   # clust_assign <- new_change$clust_assign
    #   conv <- (clust_assign != clust_assign_old) %>% sum()
    #   print(paste("i: ",iter, "| conv: ", conv))
    # }
    
  }
  
  perfect <- fct_gamma(x, y, k, N, clust_assign = sim$true, selection = "universal", alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, max_rank = 3, rank = NULL)
  # Accuracy of A2
  # Accuracy of A1
  A1 <- perfect$A[[1]]
  A1[A1 != 0] <- 1
  p1 <- as_tibble(A1) %>%
    tidyr::pivot_longer(cols = tidyselect::everything(),names_to = "col", values_to = "val") %>%
    mutate(row = rep(1:(nrow(A1)), each = ncol(A1))) %>%
    mutate(val = as.factor(val)) %>%
    ggplot(aes(x = col, y = row, fill = val)) +
    geom_tile() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggtitle("A1hat")
  
  A1_true <- sim$a[[1]]
  A1_true[A1_true != 0] <- 1
  p2 <- as_tibble(A1_true) %>%
    tidyr::pivot_longer(cols = tidyselect::everything(),names_to = "col", values_to = "val") %>%
    mutate(row = rep(1:(nrow(A1_true)), each = ncol(A1_true))) %>%
    mutate(val = as.factor(val)) %>%
    ggplot(aes(x = col, y = row, fill = val)) +
    geom_tile() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggtitle("A1 True")
  
  # Accuracy of A2
  A2 <-perfect$A[[2]]
  A2[A2 != 0] <- 1
  p3 <- as_tibble(A2) %>%
    tidyr::pivot_longer(cols = tidyselect::everything(),names_to = "col", values_to = "val") %>%
    mutate(row = rep(1:(nrow(A2)), each = ncol(A2))) %>%
    mutate(val = as.factor(val)) %>%
    ggplot(aes(x = col, y = row, fill = val)) +
    geom_tile() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggtitle("A2hat")
  
  A2_true <- sim$a[[2]]
  A2_true[A2_true != 0] <- 1
  p4 <- as_tibble(A2_true) %>%
    tidyr::pivot_longer(cols = tidyselect::everything(),names_to = "col", values_to = "val") %>%
    mutate(row = rep(1:(nrow(A2_true)), each = ncol(A2_true))) %>%
    mutate(val = as.factor(val)) %>%
    ggplot(aes(x = col, y = row, fill = val)) +
    geom_tile() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggtitle("A2 True")
  
  plots <- list(p1,p2,p3,p4)
  do.call(grid.arrange, list(grobs = plots, ncol=2))
  
  # close pdf
  dev.off()
}


# Set up simulation
N <- 100
k <- 2
rho <- 0.9
sigma <- 1
p <- 200
m <- 250
s_x <- 10
s_y <- 10
rank <- 1
b <- 1


M <- 5
times <- rep(0,M)
aris <- rep(0,M)
best_lik <- rep(0,M)
true_lik <- rep(0,M)
set.seed(1)
for (i in 1:M){
  set.seed(i)
  print(paste("M =",i))
  sim <- fct_simulate(N, k, rho, sigma, p, m, s_x, s_x, rank, b, type = "2")
  start_time <- Sys.time()
  model <- hthmix(sim$x, sim$y, k = k, nstart = 25, selection = "universal", y_sparse = TRUE, maxiter = 100)
  comp_time <- Sys.time()
  clust <- model$assign
  true <- sim$true
  clust_reorder <- clue::solve_LSAP(table(true, clust), maximum = TRUE)[clust]
  best_lik[i] <- fct_j_lik(sim$x, sim$y, k, model$assign)
  true_lik[i] <- fct_j_lik(sim$x, sim$y, k, sim$true)
  
  # out_string <- paste0("/N/project/zhangclab/alex/hthmix/iu_2022_09_16_01/output/iu_hthmix_",N, "_",k,"_", rho,"_",sigma,"_",p,"_",m,"_",s_x,"_",s_x,"_",rank,"_",b,".rds")
  times[i] <- difftime(comp_time, start_time, units = "secs")
  aris[i] <- mclust::adjustedRandIndex(true, clust_reorder)
  print(aris[i])
}
mean(aris)

sim$true
set.seed(1)
# sim <- fct_simulate(N, k, rho, sigma, p, m, s_x, s_x, rank, b, type = "2")
sim <- fct_simulate(N, k, rho, sigma, p, m, s_x, s_y, rank, b, type = "2")
sigma_hat <- fct_sigma(sim$y, nrow(sim$y), ncol(sim$y))
# lam_univ <- fct_lambda(sigma_hat, ncol(sim$x), nrow(sim$x))
lam_univ <- 55
model <- fct_sarrs(sim$y[sim$true==2,], sim$x[sim$true==2,], 1, lam_univ, sigma = sigma_hat)$Ahat[-(ncol(sim$x[sim$true==2,])+1),]
which(sim$a[[2]]!=0)
which(model!=0)
op <- which(sim$a[[1]]!=0)
on <- which(sim$a[[1]]==0)
ep <- which(model!=0)
en <- which(model==0)
en <- which(model==0)
tp <- sum(ep %in% op)
tn <- sum(en %in% on)
fp <- sum(!(ep %in% op))
fn <- sum(!(en %in% on))
pos <- length(op)
neg <- length(on)

acc <- (tp+tn)/(pos+neg)
tpr <- tp/pos
tnr <- tn/neg
f1 <- 2*tp/(2*tp+fp+fn)

mean(acc)
mean(tpr)
mean(tnr)
mean(f1)
