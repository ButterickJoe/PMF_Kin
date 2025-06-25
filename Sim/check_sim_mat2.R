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


sim_df_YS <- readRDS("data/YS_sim_py2.Rds")
sim_df_YS2 <- sim_df_YS %>% filter(Age_Focal < 100) %>%
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
## put into nice data frame
sim_df_YS3 <- sim_df_YS2 %>%
  dplyr::transmute(age_Focal = Age_Focal,
                   expectation = expectation,
                   variance = variance,
                   kin = "Younger Sisters",
                   method = "Microsimulation",
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

mat_df_YS <- acc_ys_list
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

mat_df_YS$number
## sim and matrix for younger sisters #####
YS_comp_df <- rbind(sim_df_YS3,mat_df_YS1) %>% dplyr::filter(age_Focal != 0)

YS_comp_df %>% ggplot(aes(x = age_Focal, color = method)) + geom_point(aes(y = expectation)) +
  geom_line(aes(y = variance^(1/2)))


ggsave("YS_age_mom_sim.png", input_)

############################




Q_matrix(21, 7, f_mat)
Q_matrix_joe(21, 7, f_mat)

for(i in 14:55){
  c <- Q_matrix_joe(i, 7, f_mat)
  print(colSums(c))
  print(paste0(which(c<0),"for age ", i))
}

Q_matrix(14, 7, f_mat)
mum_fert <- fert_dists(f_mat, 7)
mum_fert[[49]]

convolve(convolve(convolve(convolve(convolve(mum_fert[[14]],mum_fert[[14]]),mum_fert[[14]]),mum_fert[[14]]),mum_fert[[14]]),mum_fert[[14]])
Q_matrix_joe(13,7,f_mat)
m_pdf_14 <- list()
for(i in 1 : 7){
  m_pdf_14[[i]] <- mum_fert[[14]]
}

fourier_fft_nth(7, m_pdf_14)
plot(convolve(convolve(convolve(convolve(convolve(mum_fert[[14]],mum_fert[[14]]),mum_fert[[14]]),mum_fert[[14]]),mum_fert[[14]]),mum_fert[[14]]))
lines(abs(fourier_fft_nth(7, m_pdf_14)))
sum(abs(fourier_fft_nth(7, m_pdf_14)))
sum(fourier_fft_nth(7, m_pdf_14))
m_pdf_14


###################################

sim_list <- list()
for(sim_number in unique(sim_df_YS$sim_no)){
  age_list_df <- list()
  test_df <- sim_df_YS %>% filter(sim_no == sim_number)
  for(i in 1:99){

  to_create <- data.frame(Age = i,
                          number = seq(0,6),
                          freq_number = rep(0,7))

  test_1 <- filter(test_df, Age_Focal == i)
  number_sim <- test_1 %>% group_by(Focal_ID) %>%
    summarise(ts = sum(Alive)) %>%
    ungroup() %>%
    mutate(zero = sum(ts==0),
           one = sum(ts==1),
           two = sum(ts == 2),
           three = sum(ts == 3),
           four = sum(ts==4),
           five = sum(ts==5),
           six = sum(ts == 6))

  to_create$freq_number <- to_create$freq_number + c(number_sim$zero[1],number_sim$one[1],
                                                     number_sim$two[1],number_sim$three[1],
                                                     number_sim$four[1],number_sim$five[1],number_sim$six[1])

  age_list_df[[(1+length(age_list_df))]] <- to_create

  }


age_list_df_sim <- do.call("rbind", age_list_df)
age_list_df_sim$sim <- sim_number
sim_list[[(1+length(sim_list))]] <- age_list_df_sim
}

t_sims <- do.call("rbind",sim_list)
t_sims$sim%>%unique()
t_sims1 <- t_sims %>% group_by(sim, number, Age) %>%
  summarise(total_number_over_sims = sum(freq_number),
            number = number) %>%
  ungroup() %>%
  group_by(Age, number) %>%
  summarise(tot = sum(total_number_over_sims)) %>%
  ungroup() %>%
  group_by(Age) %>%
  summarise(normaliser = sum(tot),
            prob = tot/normaliser,
            number) %>%
  ungroup()

t_sims1$X <- t_sims1$number*t_sims1$prob
t_sims1$X2 <- t_sims1$number^2*t_sims1$prob



t_sims1 %>% group_by(Age) %>%
  summarise(expectation = sum(X),
            variance = sum(X2) - expectation^2) %>%
  ungroup() %>% ggplot(aes(x = Age)) + geom_line(aes(y = expectation), color = "blue") +
  geom_line(aes(y = variance^(0.5)), color = "red") +
  geom_point(data = mat_df_YS1%>%transmute(Age = age_Focal, expecation = expectation, variance = variance),
             aes(y = variance^(0.5)), color = "red") +
  geom_point(data = mat_df_YS1%>%transmute(Age = age_Focal, expecation = expectation, variance = variance),
             aes(y = expecation), color = "blue")


age_list_df
to_create
test_1
test_1 %>% group_by(Focal_ID) %>%
  summarise(ts = sum(Alive)) %>%
  ungroup() %>%
  mutate(zero = sum(ts==0),
         one = sum(ts==1),
         two = sum(ts == 2),
         three = sum(ts == 3),
         four = sum(ts==4),
         five = sum(ts==5),
         six = sum(ts == 6))
