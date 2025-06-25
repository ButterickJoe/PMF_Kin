library(readxl)
library(readr)
library(dplyr)
library(MASS)
library(tidyverse)
library(ggplot2)

full_simulation_var1 <- readRDS("data/YS_sim_R.Rds")
full_simulation_var1$sim_no %>% unique()
rm(full_simulation_var1)
gc()


full_simulation_var1 %>% dplyr::filter(sim_no < 10) %>%
  dplyr::group_by(age_Focal, sim_no) %>%
  dplyr::summarise(number = sum(alive)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number, age_Focal) %>%
  dplyr::summarise(rel_number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(age_Focal) %>%
  dplyr::summarise(normaliser = sum(rel_number),
                   prob = rel_number/normaliser,
                   number = number) %>%
  dplyr::ungroup()


check <- full_simulation_var1 %>%
  dplyr::group_by(age_Focal, sim_no) %>%
  dplyr::summarise(number = sum(alive)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number, age_Focal) %>%
  dplyr::summarise(rel_number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(age_Focal) %>%
  dplyr::summarise(normaliser = sum(rel_number),
                   prob = rel_number/normaliser,
                   number = number) %>%
  dplyr::ungroup() %>%
  group_by(age_Focal) %>%
  dplyr::summarise(X = number*prob,
                   X2 = number^2*prob,
                   expectation = sum(X),
                   variance = sum(X2) - expectation^2,
                   prob = prob,
                   cum_p = cumsum(prob),
                   number = number) %>%
  ungroup()
check$method <- "Mircrosimulation (R)"
check %>% head()
check %>% dplyr::select(age_Focal, expectation, variance) %>%
  ggplot(aes(x = age_Focal)) + geom_point(aes(y = expectation), color = "blue") +
  geom_point(aes(y = variance^(0.5)), color = "red")

check2 <- readRDS("data/YS_pmf_accum.Rds")

check2 %>% head()
check2$X <- check2$prob*check2$number
check2$X2 <- check2$prob*check2$number^2
check3 <- check2 %>%
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

check3 %>% head()
check %>% head()
rbind(check %>% dplyr::select(age_Focal, expectation, variance, method),
      check3 %>% dplyr::select(age_Focal, expectation, variance, method)) %>%
  ggplot(aes(x = age_Focal)) + geom_point(aes(y = expectation, color = method)) +
  geom_point(aes(y = variance^(0.5), color = method)) + ylim(c(0,2))

### check matrix

sim_df_YS <- readRDS("data/YS_sim_py2.Rds")
#mat_df_YS <- acc_ys_list
sim_df_YS %>% head()

sim_df_YS2 <- sim_df_YS %>%
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
                   prob = prob,
                   cum_p = cumsum(prob),
                   number = number) %>%
  ungroup()
sim_df_YS2 %>% head()
sim_df_YS3 <- sim_df_YS2 %>%
  dplyr::transmute(age_Focal = Age_Focal,
                   expectation = expectation,
                   variance = variance,
                   kin = "Younger Sisters",
                   method = "Microsimulation (py)",
                   inv_cum_p_95 = ifelse(cum_p >= 0.95 , number, 20 ),
                   inv_cum_p_05 = ifelse(cum_p >= 0.05 , number, 20 )) %>%
  group_by(age_Focal) %>%
  summarise(number_q95 = min(inv_cum_p_95),
            number_q05 = min(inv_cum_p_05),
            age_Focal = age_Focal,
            expectation = expectation,
            variance = variance,
            kin = "Younger Sisters",
            method = "Microsimulation")


head(sim_df_YS3)
head(check)
head(check3)


mat_res <- check3 %>% transmute(age_Focal = age_Focal, expectation = expectation, variance = variance, method = method)
py_sim_res <- sim_df_YS3 %>% transmute(age_Focal = age_Focal, expectation = expectation, variance = variance, method = method)
r_sim_res <- check %>% transmute(age_Focal = age_Focal, expectation = expectation, variance = variance, method = method)

library(reshape2)
fdf <- rbind(mat_res,py_sim_res,r_sim_res) %>% melt(id = c("age_Focal","method")) %>%
  transmute(age_Focal = age_Focal,
            moment = variable,
            value = value,
            method = method)


fdf %>% filter(age_Focal != 0) %>% ggplot(aes(x = age_Focal, color = method, linetype = moment)) +
  geom_line(aes(y = value), size = 1) + theme(legend.position = "top")
