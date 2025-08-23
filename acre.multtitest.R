
source("~/SCR_func.R")

# D_vals <- 6:11
# n_sims <- 500
 sim_results <- lapply(D_vals, function(D_val) {
  cat("Running simulation for D =", D_val, "\n")
  acre_model.sim(n_sims,
                 D = D_val,
                 g0 = 0.75,
                 sigma = 60,
                 xlim = c(-500, 900),
                 ylim = c(-500, 900),
                 traps = traps.1,
                 buffer = 300)
})


processed_list <- list()
tidylist <- list()
for (i in seq_along(sim_results)) {
  # Extract and transform
  wide <- sim_results[[i]]$wide %>%
    mutate(
      p.g0 = ((g0 - 0.75)/0.75) * 100, 
      p.D = ((D - (i + 5))/(i + 5)) * 100, 
      p.sig = ((sigma - 60)/60) * 100
    ) %>%
    pivot_longer(
      cols = c("p.g0", "p.D", "p.sig"),
      names_to = "parameter",
      values_to = "estimate"
    ) %>%
    mutate(sim_id = i+5) %>%
    select(-c("D", "g0", "sigma"))
  tidy <- sim_results[[i]]$tidy %>% 
    mutate(sim_id = i + 5)
  tidylist[[i]] = tidy
  processed_list[[i]] <- wide
}

df <- bind_rows(tidylist)
df 
write_csv(df, "multitidy_scr.csv")

# final_df <- bind_rows(processed_list)
# write_csv(final_df, "pctdiff_multi.csv")


