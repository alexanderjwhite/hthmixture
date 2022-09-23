#' Title
#'
#' @param N doc
#' @param k doc
#' @param rho doc
#' @param sigma doc
#' @param p doc
#' @param m doc
#' @param s_x doc
#' @param s_y doc
#' @param rank doc
#' @param b doc
#'
#' @import dplyr stats purrr
#'
#' @return doc
#' @export
#'
fct_simulate <- function(N = 200, k = 3, rho = 0, sigma = 1, p = 100, m = 50, s_x = 10, s_y = 15, rank = 2, b = 1, type = "2", h = 0.2, case = "independent"){
  prob <- rep(1/k,k)
  int <- prob %>% cumsum()
  rand_assign <- stats::runif(N)

  clust_assign_true <- (rand_assign) %>%
    purrr::map_int(.f = function(.x){
      clust <- (.x <= int) %>%
        which() %>%
        min()
      return(clust)
    }) %>%
    sort()

  clust_assign_true_key <- clust_assign_true %>%
    dplyr::tibble() %>%
    dplyr::mutate(order = 1:N) %>%
    dplyr::arrange((.))

  clust_assign_true_vec <- clust_assign_true_key %>%
    dplyr::pull(order)


  n <- clust_assign_true %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(value) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::pull(n)

  clust_iter <- 1
  clust_min_x <- 1
  clust_max_x <- s_x
  clust_min_y <- 1
  clust_max_y <- s_y

  x <- NULL
  y <- NULL
  a <- NULL
  for(i in 1:k){
    a_rows <- clust_min_x:clust_max_x
    a_cols <- clust_min_y:clust_max_y
    clust_iter <- clust_iter + 1
    clust_min_x <- clust_max_x+1
    clust_max_x <- clust_iter*s_x
    clust_min_y <- clust_max_y+1
    clust_max_y <- clust_iter*s_y

    if(type == "2"){
      sim <- fct_sim_tsvd(n = n[i], p = p, m = m, b = b, rank = rank, h = h, case = case)
    } else {
      sim <- fct_sim_mixrrr(n[i],a_rows,a_cols,p,m,rank,rho,sigma,b)
    }

    a <- c(a,list(sim$A))
    x <- x %>% rbind(sim$X)
    y <- y %>% rbind(sim$Y)
  }
  x <- x[clust_assign_true_vec,]
  y <- y[clust_assign_true_vec,]

  list(x = x, y = y, a = a, true = clust_assign_true)
}
