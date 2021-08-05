paths <- .libPaths()
paths <- c("/geode2/home/u100/whitealj/BigRed3/r_4_0_4_library/",paths)
.libPaths(paths)
library(hthmixture)
load("/geode2/home/u100/whitealj/BigRed3/scripts/COAD_tumor_logFPKM_Methylation.RData")

x <- COAD_Methylation450 %>% 
  tidyr::drop_na() %>% 
  t() %>% 
  dplyr::as_tibble()

y <- COAD_logFPKM %>% 
  t() %>% 
  dplyr::as_tibble()

safe_hthmix <- purrr::safely(hthmix_cv)
model <- safe_hthmix(x, y, k_min = 2, k_max = 8, r_min = 1, r_max = 2, chains = 1, penal_search = 1:100/100)

model %>% saveRDS("rna_search_results_20210726.rds")