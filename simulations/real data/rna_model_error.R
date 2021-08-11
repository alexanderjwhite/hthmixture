paths <- .libPaths()
paths <- c("/geode2/home/u100/whitealj/Carbonate/r_4_0_4_library/",paths)
.libPaths(paths)
library(hthmixture)
options(error=traceback)
load("/geode2/home/u100/whitealj/Carbonate/scripts/COAD_tumor_logFPKM_Methylation.RData")

x <- COAD_Methylation450 %>% 
  tidyr::drop_na() %>% 
  t() %>% 
  dplyr::as_tibble()

y <- COAD_logFPKM %>% 
  t() %>% 
  dplyr::as_tibble()

model <- hthmix_cv(x, y, k_min = 2, k_max = 6, r_min = 1, r_max = 2, chains = 1, penal_search = 1:100/100)

model %>% saveRDS("rna_results_20210802.rds")