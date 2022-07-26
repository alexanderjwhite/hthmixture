#' Split data for Cross Validation
#'
#' @param x matrix
#' @param y matrix
#' @param val_frac fraction for validation
#'
#' @return training and testing sets
#' @export
#' 
#' @import dplyr
#'
fct_data_split <- function(x, y, val_frac){
  
  n <- x %>% 
    dplyr::tibble() %>% 
    nrow()
  
  test_rows <- (1:n) %>% 
    dplyr::tibble() %>% 
    dplyr::slice_sample(prop = val_frac) %>% 
    dplyr::pull()
  
  train_rows <- (1:n) %>% 
    dplyr::tibble() %>% 
    dplyr::filter(!((.) %in% test_rows)) %>% 
    dplyr::pull()
  
  n_train <- train_rows %>% 
    length()
  
  n_test <- test_rows %>% 
    length()
  
  x_train <- x %>% 
    dplyr::tibble() %>% 
    dplyr::slice(train_rows) %>% 
    as.matrix()
  
  x_test <- x %>% 
    dplyr::tibble() %>% 
    dplyr::slice(test_rows) %>% 
    as.matrix()
  
  y_train <- y %>% 
    dplyr::tibble() %>% 
    dplyr::slice(train_rows) %>% 
    as.matrix()
  
  y_test <- y %>% 
    dplyr::tibble() %>% 
    dplyr::slice(test_rows) %>% 
    as.matrix()
  
  return(list(x_train = x_train,
              x_test = x_test,
              y_train = y_train,
              y_test = y_test,
              n_train = n_train,
              n_test = n_test))
  
}
