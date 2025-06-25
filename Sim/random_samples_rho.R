


sim_df_MUM <- readRDS("data/MUM_sim_py3.Rds")
n_MAX <- sim_df_MUM$sim_no %>% max()

sim_df_MUM %>% head()

sim_df_MUM %>% filter(Age_Focal == 0)

sim_df_MUM %>% filter(Age_Focal == 0) %>% group_by(Focal_ID, sim_no, Age_Kin) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number, Focal_ID, sim_no, Age_Kin) %>%
  dplyr::summarise(rel_number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number, Age_Kin) %>%
  dplyr::summarise(occurance = sum(rel_number)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number) %>%
  dplyr::summarise(total_occurance = sum(occurance),
                   prob_occurance = occurance/total_occurance,
                   Age_Kin = Age_Kin) %>%
  dplyr::ungroup() %>% ggplot(aes(x = Age_Kin, y = prob_occurance)) + geom_point()

sim_df_MUM %>% filter(Age_Focal == 0) %>% group_by(Focal_ID, sim_no, Age_Kin) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number, Focal_ID, sim_no, Age_Kin) %>%
  dplyr::summarise(rel_number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number, Age_Kin) %>%
  dplyr::summarise(occurance = sum(rel_number)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number) %>%
  dplyr::summarise(total_occurance = sum(occurance),
                   prob_occurance = occurance/total_occurance,
                   Age_Kin = Age_Kin) %>%
  dplyr::ungroup()

sample_vec <- c(5, 10, 15, 20, 30, 50, 100)

rho_df <- data.frame()
for(j in sample_vec){
  for(i in seq(1, 6)){
  test_sm <- sim_df_MUM %>% filter(sim_no %in% sample(seq(1, n_MAX) , j)) %>%
    filter(Age_Focal == 0) %>% group_by(Focal_ID, sim_no, Age_Kin) %>%
    dplyr::summarise(number = sum(Alive)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(number, Focal_ID, sim_no, Age_Kin) %>%
    dplyr::summarise(rel_number = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(number, Age_Kin) %>%
    dplyr::summarise(occurance = sum(rel_number)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(number) %>%
    dplyr::summarise(total_occurance = sum(occurance),
                     prob_occurance = occurance/total_occurance,
                     Age_Kin = Age_Kin) %>%
    dplyr::ungroup()

  test_sm$replication <- i
  test_sm$sample_size <- j
  rho_df <- rbind(rho_df, test_sm)
  }
}


rho_plot <- rho_df %>% ggplot(aes(x = Age_Kin, prob_occurance)) + facet_grid(sample_size ~ replication) +
   geom_point() + ylab("Probability of mother's age") + xlab("Mother's age") +
  ggtitle("Rows = sample sizes, Cols  = independent realisations")

rho_plot
file_put <- here::here("Sim/")
ggsave(paste0(file_put, "/rho_smaples.png"), rho_plot)

