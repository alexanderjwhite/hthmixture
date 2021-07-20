library(hthmixture)
load("G:/My Drive/Dissertation/HTH Mixture/hthmixture/simulations/real data/Data.RData")
model <- hthmix_cv(dat.mrm, dat.vol.std, k_min = 2, k_max = 8, r_min = 1, r_max = 2, chains = 1, maxiter = 1e3, penal_search = 1:100/100)
model %>% saveRDS("G:/My Drive/Dissertation/HTH Mixture/hthmixture/simulations/real data/20210713_model.rds")

model %>% 
  purrr::map_dfr(.f = function(.x){
    return(tibble(k = .x$k, r = .x$r, error = .x$error))
  })
