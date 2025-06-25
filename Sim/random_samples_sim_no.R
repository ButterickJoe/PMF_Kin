f_mat <- readr::read_rds("data/f_mat.Rds")
u_mat <- readr::read_rds("data/u_mat.Rds")
########### Accumulated kin total
f_int <- which(mothers_age(u_mat,f_mat)>0)%>%length()

########## Younger sisters over Focal's life #################
alt_rho <- readRDS("data/alt_rhob1")
alt_rho <- rbind(alt_rho, data.frame(Age = 51, rhob1 = 0, method = "MS"))
alt_rho$Age
alt_rho <- alt_rho[order(alt_rho$Age, decreasing = FALSE),]
alt_rho
rhob1 <- alt_rho$rhob1
age_b1 <- alt_rho$Age
rhob1 %>% length()
age_b1 %>% length()
acc_ys_list <- list()
for(age_focal in 1:11){
  ys_conv <- list()
  for(age_ys in max(0,age_focal-f_int-1):(age_focal-1)){

    YS <- Matrix_func_age_YS_alt(age_focal, age_ys, u_mat, f_mat, 8, rhob1, age_b1)[[1]] ## note both give same results
    ys_conv[[(1+length(ys_conv))]] <- YS
  }
  #pp <- convoluion_nth(length(ys_conv), ys_conv)
  pp <- convoluion_nth_Joe(length(ys_conv), ys_conv) ## fft quicker than nth convolution function
  ys_df <- data.frame(number = seq(0,7),
                      prob = pp,
                      age_focal = age_focal)
  acc_ys_list[[(1+length(acc_ys_list))]] <- ys_df
}
acc_ys_list <- do.call("rbind", acc_ys_list)
## check sum to one
acc_ys_list %>% dplyr::group_by(age_focal) %>%
  dplyr::summarise(check = sum(prob)) %>%
  dplyr::ungroup()


sim_df_YS <- readRDS("data/YS_sim_py3.Rds")
n_MAX <- sim_df_YS$sim_no %>% max()

sample_vec <- c(3, 5, 20, 50, 100, 200)

data_sim_random <- data.frame()
for(j in sample_vec){
for(i in seq(1, 6)){
test <- sim_df_YS %>% filter(sim_no %in% sample(seq(1, n_MAX) , j)) %>%
  dplyr::group_by(Focal_ID, Age_Focal, sim_no) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number, sim_no, Focal_ID, Age_Focal) %>%
  dplyr::summarise(number = mean(number),
                   freq_number = n()) %>%
  dplyr::group_by(number, Age_Focal) %>%
  dplyr::summarise(rel_number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Age_Focal) %>%
  dplyr::summarise(normaliser = sum(rel_number),
                   prob = rel_number/normaliser,
                   number = number) %>%
  dplyr::ungroup() %>%
  group_by(Age_Focal) %>%
  dplyr::summarise(X = number*prob,
                   X2 = number^2*prob,
                   expectation = sum(X),
                   variance = sum(X2) - expectation^2,
                   prob = prob) %>%
  ungroup()
test$indend_sample <- i
test$sample_size <- j
data_sim_random <- rbind(data_sim_random, test)
}
}

sis_samples <- data_sim_random %>% ggplot(aes(x = Age_Focal)) +
  facet_grid(sample_size~indend_sample) +
  geom_point(aes(y = expectation), color = "blue") + geom_point(aes(y = variance^(0.5)), color = "red") +
  geom_line(data = mat_df_YS1%>%transmute(Age_Focal = age_Focal, variance = variance),
            aes(y = variance^(0.5)), color = "black")

sis_samples
file_put <- here::here("Sim/")
ggsave(paste0(file_put, "/YS_smaples.png"), sis_samples)


mat_df_YS <- readRDS("data/YS_pmf_accum.Rds")
mat_df_YS %>% head()
mat_df_YS$X <- mat_df_YS$number*mat_df_YS$prob
mat_df_YS$X2 <- mat_df_YS$number^2*mat_df_YS$prob

mat_df_YS1 <- mat_df_YS %>%
  dplyr::group_by(age_focal) %>%
  dplyr::summarise(expectation = sum(X),
                   s2 = sum(X2),
                   variance = s2 - expectation^2,
                   cum_p = cumsum(prob),
                   number = number) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_Focal = age_focal,
                   expectation = expectation,
                   variance = variance,
                   kin = "Younger Sisters",
                   method = "PMF model",
                   inv_cum_p_95 = ifelse(cum_p >= 0.95 , number, 20 ),
                   inv_cum_p_05 = ifelse(cum_p >= 0.05 , number, 20 )) %>%
  group_by(age_Focal) %>%
  summarise(number_q95 = min(inv_cum_p_95),
            number_q05 = min(inv_cum_p_05),
            age_Focal = age_Focal,
            expectation = expectation,
            variance = variance,
            kin = "Younger Sisters",
            method = "PMF model")
