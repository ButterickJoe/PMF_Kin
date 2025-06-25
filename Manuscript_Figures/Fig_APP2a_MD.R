
library(dplyr)
library(ggplot2)

fig_out <- here::here("Manuscipt_Figures_code/Figures")

### First: load in the ABM model results
sim_df_DAU <- readRDS("data/DAU_sim_py2.Rds")
sim_df_MUM <- readRDS("data/MUM_sim_py2.Rds")
#### Second: load in rates

### check matrix
f_mat <- readr::read_rds("data/f_mat.Rds")
u_mat <- readr::read_rds("data/u_mat.Rds")
f_int <- which(mothers_age(u_mat,f_mat)>0)%>%length()

##### Age pdfs #########

############################## Daughters #####################################################
## Focal age 50
joe_dau_50 <- list()
for(i in 0:49){
  val <- unconditional_Focals_age_Dau(50, i, u_mat, f_mat, 7)
  df <- data.frame(age = i,
                   number = seq(0,6),
                   prob = val)
  joe_dau_50[[(1+length(joe_dau_50))]] <- df
}
joe_dau_50 <- do.call("rbind", joe_dau_50) %>% as.data.frame()
joe_dau_50$X <- joe_dau_50$number*joe_dau_50$prob
joe_dau_50$X2 <- joe_dau_50$number^2*joe_dau_50$prob
joe_dau_50_means <- joe_dau_50 %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   var = s2 - mean^2) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 50 years",
                kin = "Daughters")
## Focal age 20
joe_dau_20 <- list()
for(i in 0:19){
  val <- unconditional_Focals_age_Dau(20, i, u_mat, f_mat, 7)
  df <- data.frame(age = i,
                   number = seq(0,6),
                   prob = val)
  joe_dau_20[[(1+length(joe_dau_20))]] <- df
}
joe_dau_20 <- do.call("rbind", joe_dau_20) %>% as.data.frame()
joe_dau_20$X <- joe_dau_20$number*joe_dau_20$prob
joe_dau_20$X2 <- joe_dau_20$number^2*joe_dau_20$prob
joe_dau_20_means <- joe_dau_20 %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   var = s2 - mean^2) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 20 years",
                kin = "Daughters")

########################################## Mother ##############################################
##  Focal age 50
joe_mum_50 <- list()
for(i in 51:99){
  val <- unconditional_mother_pdf(50, i, u_mat, f_mat, 7)
  joe_mum_50[[(1+length(joe_mum_50))]] <- val
}
joe_mum_50 <- do.call("rbind", joe_mum_50) %>% as.data.frame()
joe_mum_50 %>% head()
joe_mum_50$X <- joe_mum_50$number*joe_mum_50$prob
joe_mum_50$X2 <- joe_mum_50$number^2*joe_mum_50$prob
joe_mum_50_means <- joe_mum_50 %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   var = s2 - mean^2) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 50 years",
                kin = "Mothers")
## Focal age 20
joe_mum_20 <- list()
for(i in 21:99){
  val <- unconditional_mother_pdf(20, i, u_mat, f_mat, 7)
  joe_mum_20[[(1+length(joe_mum_20))]] <- val
}
joe_mum_20 <- do.call("rbind", joe_mum_20) %>% as.data.frame()
joe_mum_20 %>% head()
joe_mum_20$X <- joe_mum_20$number*joe_mum_20$prob
joe_mum_20$X2 <- joe_mum_20$number^2*joe_mum_20$prob
joe_mum_20_means <- joe_mum_20 %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   var = s2 - mean^2) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 20 years",
                kin = "Mothers")

### combined
dat_20 <- rbind(joe_mum_20_means,joe_dau_20_means) %>% dplyr::mutate(xint = 20)
dat_50 <- rbind(joe_mum_50_means,joe_dau_50_means) %>% dplyr::mutate(xint = 50)

legend_title_a <- "Age of Focal"
## Figure 5 in manuscript
md_age_plot <- rbind(dat_20,dat_50) %>% dplyr::filter(mean != 0) %>%
  ggplot(aes(x = age, y = mean, color = age_Focal)) +
  facet_wrap(~kin, scales = "free_x")+ geom_point(size = 1) +
  geom_vline(aes(xintercept = xint, color = age_Focal), linetype="dashed", size = 0.75) + theme_bw() +
  geom_errorbar(aes(ymin = ifelse(mean - var^(0.5) < 0, 0, mean - var^(0.5)), ymax = mean + var^(0.5)), size = 0.25) +
  ylab("Expectation \u00B1 standard deviation") + xlab("Age of kin") +
  scale_color_brewer(legend_title_a, palette = "Set2")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size=14)) + theme(strip.text = element_text(size = 14)) +
  theme(legend.position = "top")

#ggsave(paste0(fig_out,"/Age_dists_mum_dau.png"), md_age_plot, width = 9, height = 5)


##### accumulated kin by age 50 compared to microsimulation

### Daughter

## Focal 50
df_list_dau_50 <- joe_dau_50 %>% as.data.frame()
conv_list_dau_50 <- list()
for(aged in df_list_dau_50$age %>% unique()){
  p <- df_list_dau_50 %>% dplyr::filter(age == aged)
  p <- p$prob
  conv_list_dau_50[[(1+length(conv_list_dau_50))]] <- list("data" = p, "age" = aged)
}
conv_list_new_dau_50 <- lapply(1:length(conv_list_dau_50), function(x){ conv_list_dau_50[[x]][["data"]] })
pdf_f50_dau <- convoluion_nth(length(conv_list_new_dau_50), conv_list_new_dau_50)

acc_Dau_50 <- data.frame(age_Focal = "At 50 years",
                         number = seq(0,6),
                         prob = pdf_f50_dau,
                         kin = "Daughters")
acc_Dau_50$X <- acc_Dau_50$number*acc_Dau_50$prob
acc_Dau_50$X2 <- acc_Dau_50$number^2*acc_Dau_50$prob
mean_Dau_50 <- sum(acc_Dau_50$X)
var_Dau_50 <- sum(acc_Dau_50$X2) - (mean_Dau_50)^2
acc_Dau_50$expectation <- mean_Dau_50
acc_Dau_50$variance <- var_Dau_50

## Focal 20
df_list_dau_20 <- joe_dau_20 %>% as.data.frame()
conv_list_dau_20 <- list()
for(aged in df_list_dau_20$age %>% unique()){
  p <- df_list_dau_20 %>% dplyr::filter(age == aged)
  p <- p$prob
  conv_list_dau_20[[(1+length(conv_list_dau_20))]] <- list("data" = p, "age" = aged)
}
conv_list_new_dau_20 <- lapply(1:length(conv_list_dau_20), function(x){ conv_list_dau_20[[x]][["data"]] })
pdf_f20_dau <- convoluion_nth(length(conv_list_new_dau_20), conv_list_new_dau_20)

acc_Dau_20 <- data.frame(age_Focal = "At 20 years",
                         number = seq(0,6),
                         prob = pdf_f20_dau,
                         kin = "Daughters")
acc_Dau_20$X <- acc_Dau_20$number*acc_Dau_20$prob
acc_Dau_20$X2 <- acc_Dau_20$number^2*acc_Dau_20$prob
mean_Dau_20 <- sum(acc_Dau_20$X)
var_Dau_20 <- sum(acc_Dau_20$X2) - (mean_Dau_20)^2
acc_Dau_20$expectation <- mean_Dau_20
acc_Dau_20$variance <- var_Dau_20

daughter_df_20_50 <- rbind(acc_Dau_20, acc_Dau_50)

### Mother
## Focal 20
mother_conv_f20 <- list()
for(age_mother in (20 + 2) : 99){
  mother <- unconditional_mother_pdf(20, age_mother, u_mat, f_mat, 7)
  mother <- mother %>% dplyr::filter(number==1)
  mother <- mother$prob
  mother_conv_f20[[(1+length(mother_conv_f20))]] <- mother
}
pp_20 <- sum(unlist(lapply(1:length(mother_conv_f20), function(x){mother_conv_f20[[x]]})))
mother_df_20 <- data.frame(age_Focal = "At 20 years",
                           number = c(0,1),
                           prob = c(1-pp_20,pp_20),
                           kin = "Mother")
mother_df_20$X <- mother_df_20$number*mother_df_20$prob
mother_df_20$X2 <- mother_df_20$number^2*mother_df_20$prob
mean_Mum_20 <- sum(mother_df_20$X)
var_Mum_20 <- sum(mother_df_20$X2) - (mean_Mum_20)^2
mother_df_20$expectation <- mean_Mum_20
mother_df_20$variance <- var_Mum_20

## Focal 50
mother_conv_f50 <- list()
for(age_mother in (50 + 2) : 99){
  mother <- unconditional_mother_pdf(50, age_mother, u_mat, f_mat, 7)
  mother <- mother %>% dplyr::filter(number==1)
  mother <- mother$prob
  mother_conv_f50[[(1+length(mother_conv_f50))]] <- mother
}
pp_50 <- sum(unlist(lapply(1:length(mother_conv_f50), function(x){mother_conv_f50[[x]]})))
mother_df_50 <- data.frame(age_Focal = "At 50 years",
                           number = c(0,1),
                           prob = c(1-pp_50,pp_50),
                           kin = "Mother")
mother_df_50$X <- mother_df_50$number*mother_df_50$prob
mother_df_50$X2 <- mother_df_50$number^2*mother_df_50$prob
mean_Mum_50 <- sum(mother_df_50$X)
var_Mum_50 <- sum(mother_df_50$X2) - (mean_Mum_50)^2
mother_df_50$expectation <- mean_Mum_50
mother_df_50$variance <- var_Mum_50


mother_df_20_50 <- rbind(mother_df_20, mother_df_50)

## simulation

## Focal age 50
## Daughters
age_dist_sim_DAU_50 <- sim_df_DAU %>% dplyr::filter(Age_Focal == 50)
age_dist_sim_DAU_50 <- age_dist_sim_DAU_50 %>%
  dplyr::group_by(Focal_ID, Age_Focal, sim_no) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup()
## Mother
age_dist_sim_MUM_50 <- sim_df_MUM %>% dplyr::filter(Age_Focal == 50)
age_dist_sim_MUM_50 <- age_dist_sim_MUM_50 %>%
  dplyr::group_by(Focal_ID, Age_Focal, sim_no) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup()

sim_DAU_PDF_50 <- age_dist_sim_DAU_50%>%
  dplyr::group_by(number) %>%
  dplyr::summarise(number = mean(number),
            freq_number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(relative_freq = freq_number/sum(freq_number),
         X = relative_freq*number,
         X2 = relative_freq*number^2) %>%
  dplyr::transmute(age_Focal = "At 50 years",
            number = seq(0,max(number)),
            prob = relative_freq,
            kin = "Daughters",
            X = X,
            X2 = X2,
            expectation = sum(X),
            variance = sum(X2) - expectation^2)
sim_DAU_PDF_50$method <- "Microsimulation"

sim_MUM_PDF_50 <- age_dist_sim_MUM_50 %>%
  group_by(number) %>%
  summarise(number = mean(number),
            freq_number = n()) %>%
  ungroup() %>%
  mutate(relative_freq = freq_number/sum(freq_number),
         X = relative_freq*number,
         X2 = relative_freq*number^2) %>%
  transmute(age_Focal = "At 50 years",
            number = seq(0,max(number)),
            prob = relative_freq,
            kin = "Mother",
            X = X,
            X2 = X2,
            expectation = sum(X),
            variance = sum(X2) - expectation^2)
sim_MUM_PDF_50$method <- "Microsimulation"

## Focal age 20
## Daughters
age_dist_sim_DAU_20 <- sim_df_DAU %>% dplyr::filter(Age_Focal == 20)
age_dist_sim_DAU_20 <- age_dist_sim_DAU_20 %>%
  dplyr::group_by(Focal_ID, Age_Focal, sim_no) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup()
## Mother
age_dist_sim_MUM_20 <- sim_df_MUM %>% dplyr::filter(Age_Focal == 20)
age_dist_sim_MUM_20 <- age_dist_sim_MUM_20 %>%
  dplyr::group_by(Focal_ID, Age_Focal, sim_no) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup()

sim_DAU_PDF_20 <- age_dist_sim_DAU_20%>%
  group_by(number) %>%
  summarise(number = mean(number),
            freq_number = n()) %>%
  ungroup() %>%
  mutate(relative_freq = freq_number/sum(freq_number),
         X = relative_freq*number,
         X2 = relative_freq*number^2) %>%
  transmute(age_Focal = "At 20 years",
            number = seq(0,max(number)),
            prob = relative_freq,
            kin = "Daughters",
            X = X,
            X2 = X2,
            expectation = sum(X),
            variance = sum(X2) - expectation^2)
sim_DAU_PDF_20$method <- "Microsimulation"

sim_MUM_PDF_20 <- age_dist_sim_MUM_20 %>%
  group_by(number) %>%
  summarise(number = mean(number),
            freq_number = n()) %>%
  ungroup() %>%
  mutate(relative_freq = freq_number/sum(freq_number),
         X = relative_freq*number,
         X2 = relative_freq*number^2) %>%
  transmute(age_Focal = "At 20 years",
            number = seq(0,max(number)),
            prob = relative_freq,
            kin = "Mother",
            X = X,
            X2 = X2,
            expectation = sum(X),
            variance = sum(X2) - expectation^2)
sim_MUM_PDF_20$method <- "Microsimulation"

sim_mum_20_50 <- rbind(sim_MUM_PDF_20,sim_MUM_PDF_50)
sim_dau_20_50 <- rbind(sim_DAU_PDF_20,sim_DAU_PDF_50)

mother_df_20_50 %>% head()
sim_MUM_PDF_50 %>% head()
mother_df_20_50 <- mother_df_20_50 %>%
  transmute(age_Focal = age_Focal,
            number = number,
            prob = prob,
            kin = kin,
            X = X,
            X2 = X2,
            expectation = expectation,
            variance = variance,
            method = "PMF model")

daughter_df_20_50 %>% head()
daughter_df_20_50 <- daughter_df_20_50 %>%
  transmute(age_Focal = age_Focal,
            number = number,
            prob = prob,
            kin = kin,
            X = X,
            X2 = X2,
            expectation = expectation,
            variance = variance,
            method = "PMF model")

## data for plotting
m_d_sm_pdf <- rbind(daughter_df_20_50,mother_df_20_50,sim_dau_20_50,sim_mum_20_50)
err_dat <- m_d_sm_pdf %>% filter(number == 0)
err_dat$prob <- 0.02
## Figure 6 in manuscript
m_d_sm_pdf %>%
  ggplot(aes(x = number, y = prob, color = method, fill = method)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.25) +
  facet_grid(age_Focal~kin, scales = "free_x") + theme_bw() +
  geom_vline(aes(xintercept = expectation, color = method, linetype = method)) +
  geom_errorbarh(data = err_dat,
                 aes(xmin = ifelse(expectation - variance^(0.5) < 0, 0, expectation - variance^(0.5)),
                     xmax = expectation + variance^(0.5),
                     color = method, linetype = method),
                 height = 0.015) +
  scale_color_manual(labels = c("Microsimulation","PMF model"), values = c("black","orange"))+
  scale_fill_manual(labels = c("Microsimulation","PMF model"), values = c("black","orange"))+
  ylab("Probability") + xlab("Number of kin")  +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size=14)) +
  theme(strip.text = element_text(size = 14)) +
  theme(legend.position = "top") + scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  ylim(c(0,1))

ggsave(paste0(fig_out,"/PDF_mum_dau_APP.png"), md_pdf_plot, width = 9, height = 7)

md_pdf_plot


########### new
m_d_sm_pdf %>% head()

adds <- m_d_sm_pdf %>%
  group_by(age_Focal, kin, method) %>%
  summarise(number = number,
            X3 = number^3*prob,
            s2 = sum(X2),
            s3 = sum(X3),
            m3 = s3 + 2*expectation^3 - 3*expectation*s2,
            skew = m3/((variance^(0.5))^(3)),
            expectation = expectation,
            sd = variance^(0.5)  ) %>%
  ungroup() %>%
  transmute(x_text = ifelse( (method == "Microsimulation" & kin == "Mother" & age_Focal == "At 20 years"), -0.5,
                             ifelse(  (method == "PMF model" & kin == "Mother" & age_Focal == "At 20 years") , -0.5,
                                      ifelse( (method == "Microsimulation" & kin == "Daughters")  , 1,
                                              ifelse((method == "PMF model" & kin == "Daughters") , 4,
                                                     ifelse((method == "Microsimulation" & kin == "Mother" & age_Focal == "At 50 years"), -0.5,
                                                            ifelse((method == "PMF model" & kin == "Mother" & age_Focal == "At 50 years"), 0.5, NA)))))),

            expectation = expectation,
            sd = sd,
            skew = skew,
            number = number,
            kin = kin,
            method = method,
            age_Focal = age_Focal) %>%
  as.data.frame()


md_pdf_plot <- m_d_sm_pdf %>% filter(number < 7)  %>%
  ggplot(aes(x = number, y = prob , color = method, fill = method )) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.25) +
  facet_grid( age_Focal~ kin, scales = "free_x") + theme_bw() +
  scale_color_manual(labels = c("Microsimulation","PMF model"), values = c("black","orange"))+
  scale_fill_manual(labels = c("Microsimulation","PMF model"), values = c("black","orange")) +
  ylab("Probability") + xlab("Number of kin")  +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size=14)) +
  theme(strip.text = element_text(size = 14)) +
  theme(legend.position = "top") + scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  ylim(c(0,1)) +
  geom_text(data = adds %>% dplyr::select(-skew,-sd) ,
            aes(x = x_text , y = y_text,
                label = paste0("mean", "==", round(expectation, 2))
            ), parse = TRUE, hjust = 0, vjust = 1, size = 5, show.legend = FALSE) +
  geom_text(data = adds %>% dplyr::select(-expectation,-skew),
            aes(x = x_text, y = y_text1,
                label = paste0("sd", "==", round(sd,2))
            ), parse = TRUE, hjust = 0, vjust = 1 , size = 5, show.legend = FALSE) +
  geom_text(data = adds %>% dplyr::select(-expectation,-sd),
            aes(x = x_text, y_text2,
                label = paste0("skew", "==", round(skew,2))
            ), parse = TRUE, hjust = 0, vjust = 1 , size = 5, show.legend = FALSE)



adds$y_text <- ifelse( (adds$method == "PMF model" & adds$age_Focal == "At 20 years" & adds$kin == "Mother"), 0.7, 1)
adds$y_text1 <- ifelse( (adds$method == "PMF model" & adds$age_Focal == "At 20 years" & adds$kin == "Mother"), 0.6, 0.9)
adds$y_text2 <- ifelse( (adds$method == "PMF model" & adds$age_Focal == "At 20 years" & adds$kin == "Mother") , 0.5, 0.8)
