library(dplyr)
source("./functions/20210205_sarrs_alex.R")
future::plan("multiprocess")
dims <- c(100, 200, 500, 1000, 5000)
ss <- c(10, 30, 50, 100, 500)
ns <- c(50, 100, 200, 500)
bs <- c(5,10,20,50)
r <- 2
rho <- 0
lam <- 1

results <- expand.grid(.dim = dims, .s = ss, .n = ns, .b = bs) %>% 
  as.list() %>% 
  furrr::future_pmap_dfr(
    .progress = TRUE,
    .f = function(.dim,.s,.n,.b){
      if(.s > .dim){return(tibble(dim = .dim, s = .s, n = .n, b = .b, ratio = NA))}
      m <- .dim
      p <- .dim
      a_rows <- 1:.s
      A <- matrix(0,p,m)
      nrow <- length(a_rows)
      B0 <- matrix(rnorm(nrow*r),nrow,r)
      B1 <- matrix(rnorm(r*m),r,m)
      A[a_rows,] <- .b*B0%*%B1
      
      Sigma <- matrix(1,p,p)
      for(j in 1:p) for (k in 1:p) Sigma[j,k] <- rho^abs(j-k)
      
      X <- MASS::mvrnorm(.n,rep(0,p),Sigma)
      
      E <- matrix(rnorm(.n*m),.n,m)
      
      Y <- X%*%A + E
      
      
      rhat <- 2
      estimate <- SARRS(Y,X,r, lam, "grLasso")
      
      A_sarrs <- estimate %>% 
        purrr::pluck("Ahat")
      
      
      ratio <- ((A_sarrs[1:p,]-A) %>% norm(type = "2"))/(A %>% norm(type="2"))
      
      return(tibble(dim = .dim, s = .s, n = .n, b = .b, ratio = ratio))
    }
  )

results %>% 
  tidyr::drop_na() %>% 
  readr::write_rds(file = "./results/20210205_high_dim_sarrs_sim_alex.rds")

results %>% 
  tidyr::drop_na() %>%
  filter(dim == 5000) %>% 
  arrange(ratio) %>% 
  View()


# Play with plots
# Remove NAs
# plot_list <- expand.grid(.n = ns, .b = bs) %>% 
#   as.list() %>% 
#   purrr::pmap(
#     .f = function(.n, .b){
#       results %>% 
#         tidyr::drop_na() %>% 
#         filter(n == .n & b == .b) %>% 
#         plotly::plot_ly(x = ~dim, y = ~ratio, type = "scatter", mode = "markers", color = ~s,
#                         text = ~paste("dim:",dim,"s:",s,"ratio",ratio)) %>% 
#         plotly::layout(title = paste("n =",.n,"b =",.b))
#     }
#   )
# 
# plotly::subplot(plot_list[1:8], nrows = 2)
# 
# library(ggplot2)
# results %>% 
#   tidyr::drop_na() %>% 
#   filter(n == 500 & b == 20) %>% 
#   mutate(dim_s = dim/s) %>% 
#   ggplot(aes(x=dim_s, y = ratio)) +
#   geom_point()
#   
# results %>% 
#   tidyr::drop_na() %>% 
#   filter(n == 500 & b == 20) %>% 
#   plotly::plot_ly(x = ~dim, y = ~ratio, type = "scatter", mode = "markers", color = ~s)
# 
# 
# plotList <- function(nplots) {
#   lapply(seq_len(nplots), function(x) plotly::plot_ly())
# }
