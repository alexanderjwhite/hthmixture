library(hthmixture)
library(ggplot2)
load("G:/My Drive/Dissertation/HTH Mixture/hthmixture/simulations/real data/Data.RData")
model <- hthmix_cv(dat.mrm, dat.vol.std, k_min = 2, k_max = 3, r_min = 1, r_max = 1, chains = 1, maxiter = 1e3, penal_search = 1:100/100)
model %>% saveRDS("G:/My Drive/Dissertation/HTH Mixture/hthmixture/simulations/real data/20210713_model.rds")

model %>% 
  purrr::map_dfr(.f = function(.x){
    return(tibble(k = .x$k, r = .x$r, error = .x$error))
  })

best_model_test <- hthmix(dat.mrm, dat.vol.std, k = 2, rank = 1, chains = 5, maxiter = 1e3, penal_search = 1:100/100)

# best_model_test %>% readr::write_rds("G:/My Drive/Dissertation/HTH Mixture/hthmixture/simulations/real data/bestk2.rds")
# best_model %>% readr::write_rds("G:/My Drive/Dissertation/HTH Mixture/hthmixture/simulations/real data/bestk3.rds")
best_model_test$result %>% 
  select(chain, llik, iter) %>% 
  distinct() %>% 
  arrange(-llik)
best_model <- best_model_test$result %>% filter(chain==3)
test <- best_model_test
test$result <- test$result %>% filter(chain==3)
best_assign <- test %>% purrr::pluck("result", "assign", "final_assign")
best_model %>% purrr::pluck("result","assign", "assign_store") %>% 
  group_by(iter) %>% 
  summarize(n=n())

p1 <- dat.demo %>% 
  bind_cols(tibble(clust = as.factor(best_assign))) %>% 
  ggplot() +
  geom_boxplot(aes(x = clust, y = CDRSB_bl))
  
p2 <- dat.demo %>% 
  bind_cols(tibble(clust = as.factor(best_assign))) %>% 
  ggplot() +
  geom_boxplot(aes(x = clust, y = ADAS11_bl))

p3 <- dat.demo %>% 
  bind_cols(tibble(clust = as.factor(best_assign))) %>% 
  ggplot() +
  geom_boxplot(aes(x = clust, y = ADAS13_bl))

p4 <- dat.demo %>% 
  bind_cols(tibble(clust = as.factor(best_assign))) %>% 
  ggplot() +
  geom_boxplot(aes(x = clust, y = FAQ_bl))

p5 <- dat.demo %>% 
  bind_cols(tibble(clust = as.factor(best_assign))) %>% 
  ggplot() +
  geom_boxplot(aes(x = clust, y = MMSE_bl))

p6 <- dat.demo %>% 
  bind_cols(tibble(clust = as.factor(best_assign))) %>% 
  ggplot() +
  geom_bar(aes(x = clust, fill = DX))
p6
