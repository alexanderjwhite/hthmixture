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
#' @examples
fct_data_split <- function(x, y, val_frac){
  
  n <- x %>% 
    tibble() %>% 
    nrow()
  
  test_rows <- (1:n) %>% 
    tibble() %>% 
    slice_sample(prop = val_frac) %>% 
    pull()
  
  train_rows <- (1:n) %>% 
    tibble() %>% 
    filter(!((.) %in% test_rows)) %>% 
    pull()
  
  n_train <- train_rows %>% 
    length()
  
  n_test <- test_rows %>% 
    length()
  
  x_train <- x %>% 
    tibble() %>% 
    slice(train_rows) %>% 
    as.matrix()
  
  x_test <- x %>% 
    tibble() %>% 
    slice(test_rows) %>% 
    as.matrix()
  
  y_train <- y %>% 
    tibble() %>% 
    slice(train_rows) %>% 
    as.matrix()
  
  y_test <- y %>% 
    tibble() %>% 
    slice(test_rows) %>% 
    as.matrix()
  
  return(list(x_train = x_train,
              x_test = x_test,
              y_train = y_train,
              y_test = y_test,
              n_train = n_train,
              n_test = n_test))
  
}
