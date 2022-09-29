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
fct_simulate <- function(N = 200, k = 3, sigma = 1, p = 100, m = 50, rank = 2, b = 1, d = 20, h = 0.2, case = "independent"){
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


  x <- NULL
  y <- NULL
  a <- NULL
  for(i in 1:k){
    
    sim <- fct_sim_tsvd(n = n[i], p = p, m = m, b = b, d = d, rank = rank, h = h, case = case)

    a <- c(a,list(sim$A))
    x <- x %>% rbind(sim$X)
    y <- y %>% rbind(sim$Y)
  }
  x <- x[clust_assign_true_vec,]
  y <- y[clust_assign_true_vec,]

  list(x = x, y = y, a = a, true = clust_assign_true)
}
