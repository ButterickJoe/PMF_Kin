library(ggplot2)
library(dplyr)

fig_out <- here::here("Manuscript_Figures/Figures")

# load in rates

### check matrix
f_mat <- readr::read_rds("data/f_mat.Rds")
u_mat <- readr::read_rds("data/u_mat.Rds")
f_int <- which(mothers_age(u_mat,f_mat)>0)%>%length()

############################################# Age-specific dists of kin when Focal is 20 ##############################
### Younger sisters
df_list_ys <- list()
for(i in 0:19){
  df <- ys_PMF(20, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_ys[[(1+length(df_list_ys))]] <- df1
}
df_list_ys_20 <- do.call("rbind", df_list_ys) %>% as.data.frame()
### Older sisters
df_list_os <- list()
for(i in 21:(21+f_int+1)){
  df <- os_PMF(20, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_os[[(1+length(df_list_os))]] <- df1
}
df_list_os_20 <- do.call("rbind", df_list_os) %>% as.data.frame()
joe_sis_20 <- rbind(df_list_ys_20,df_list_os_20) %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 20 years",
                kin = "Sisters")
dat_20_sis <- joe_sis_20 %>% dplyr::mutate(xint = 20)
### When Focal is 50
### Younger sisters
df_list_ys <- list()
for(i in 0:49){
  df <- ys_PMF(50, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_ys[[(1+length(df_list_ys))]] <- df1
}
df_list_ys_50 <- do.call("rbind", df_list_ys) %>% as.data.frame()
### Older sisters
df_list_os <- list()
for(i in 51:(51+f_int+1)){
  df <- os_PMF(50, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_os[[(1+length(df_list_os))]] <- df1
}
df_list_os_50 <- do.call("rbind", df_list_os) %>% as.data.frame()
joe_sis_50 <- rbind(df_list_ys_50, df_list_os_50) %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 50 years",
                kin = "Sisters")
dat_50_sis <- joe_sis_50 %>% dplyr::mutate(xint = 50)
dat_50_sis
######## Aunts #######################################################################################################
############################################# Age-specific dists of kin when Focal is 20 ##############################
### Younger cousins
df_list_cya_20 <- list()
for(i in 0:60){
  df <- cya_PMF(20, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_cya_20[[(1+length(df_list_cya_20))]] <- df1
}
df_list_cya_20 <- do.call("rbind", df_list_cya_20) %>% as.data.frame()

### Older cousins
df_list_coa_20 <- list()
for(i in 0:80){
  df <- coa_PMF(20, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_coa_20[[(1+length(df_list_coa_20))]] <- df1
}
df_list_coa_20 <- do.call("rbind", df_list_coa_20) %>% as.data.frame()

df_list_coa_20 %>%
  mutate(X = prob*number,
         X2 = prob*number^2) %>%
  group_by(age) %>%
  summarise(mean = sum(X),
            s2 = sum(X2),
            sd = (s2 - mean^2)^(0.5)) %>%
  ungroup() %>% ggplot(aes(x = age)) +
  geom_point(aes(y = mean)) + geom_errorbar(aes(ymin = ifelse(mean - sd<0,0,mean-sd), ymax = mean+sd))
################################################################################

joe_cousin_20 <- data.frame()
for(age1 in 0:80){
  df <- rbind(df_list_cya_20,df_list_coa_20) %>% dplyr::filter(age == age1)
  if(nrow(df)==7){prob = df$prob}
  else{prob = convoluion_nth(2, list(df$prob[1:7],df$prob[8:14]))}
  temp <- data.frame(number = seq(0,6),
                     prob = prob,
                     age = age1)
  joe_cousin_20 <- rbind(joe_cousin_20, temp)
}
joe_cousin_20 <- joe_cousin_20 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 20 years",
                kin = "Cousins")
dat_20_cousin <- joe_cousin_20 %>% dplyr::mutate(xint = 20)
### When Focal is 50
### Cousins from younger aunts
df_list_cya_50 <- list()
for(i in 0:100){
  df <- cya_PMF(50, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_cya_50[[(1+length(df_list_cya_50))]] <- df1
}
df_list_cya_50 <- do.call("rbind", df_list_cya_50) %>% as.data.frame()
### Cousins from older aunts
df_list_coa_50 <- list()
for(i in 0:100){
  df <- coa_PMF(50, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_coa_50[[(1+length(df_list_coa_50))]] <- df1
}
df_list_coa_50 <- do.call("rbind", df_list_coa_50) %>% as.data.frame()

joe_cousin_50 <- data.frame()
for(age1 in 0:100){
  df <- rbind(df_list_cya_50,df_list_coa_50) %>% dplyr::filter(age == age1)
  if(nrow(df)==7){prob = df$prob}
  else{prob = convoluion_nth(2, list(df$prob[1:7],df$prob[8:14]))}
  temp <- data.frame(number = seq(0,6),
                     prob = prob,
                     age = age1)
  joe_cousin_50 <- rbind(joe_cousin_50, temp)
}
joe_cousin_50 <- joe_cousin_50 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 50 years",
                kin = "Cousins")
dat_50_cousin <- joe_cousin_50 %>% dplyr::mutate(xint = 50)

legend_title_a <- "Age of Focal"
scale = 300
## Figure 5 in manuscript
sa_age_plot <- rbind(dat_20_sis, dat_50_sis, dat_20_cousin, dat_50_cousin) %>% dplyr::filter(mean != 0) %>%
  ggplot(aes(x = age, y = mean, color = age_Focal)) +
  facet_wrap(~kin, scales = "free_x")+ geom_point(size = 0.75) + theme_bw() +
  geom_line(aes(y = var^(0.5)), linetype = "dashed", size = 0.75) +
  geom_line(data = rbind(dat_20_sis, dat_50_sis, dat_20_cousin, dat_50_cousin) %>% dplyr::filter(skew < 60),
            aes(y = skew/scale), linetype = "solid", size = 0.75) +
  ylab("Expectation and Std. Dev.") + xlab("Age of kin") +
  scale_color_brewer(legend_title_a, palette = "Set2") +
  scale_y_continuous() +
  scale_y_continuous(sec.axis = sec_axis(~.*scale, name = "Skew")) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size=14)) + theme(strip.text = element_text(size = 14)) +
  theme(legend.position = "top")

sa_age_plot
ggsave(paste0(fig_out,"/Age_dists_sis_cous_REV.png"), sa_age_plot, width = 9, height = 5)

############################################## PDF of accumulated kin when Focal is 20
## Younger
df_list_ys_20 <- df_list_ys_20 %>% as.data.frame()
df_list_ys_20 %>% head()
conv_list_ys <- list()
for(age1 in df_list_ys_20$age %>% unique()){
  p <- df_list_ys_20 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_ys[[(1+length(conv_list_ys))]] <- list("data" = p, "age" = age1)
}
conv_list_new_ys <- lapply(1:length(conv_list_ys),
                           function(x){ conv_list_ys[[x]][["data"]] })
pdf_f20_ys <- convoluion_nth(length(conv_list_new_ys), conv_list_new_ys)
## Older
df_list_os_20 <- df_list_os_20 %>% as.data.frame()
df_list_os_20 %>% head()
conv_list_os <- list()
for(age1 in df_list_os_20$age %>% unique()){
  p <- df_list_os_20 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_os[[(1+length(conv_list_os))]] <- list("data" = p, "age" = age1)
}
conv_list_new_os <- lapply(1:length(conv_list_os), function(x){ conv_list_os[[x]][["data"]] })
pdf_f20_os <- convoluion_nth(length(conv_list_new_os), conv_list_new_os)

## df combine
acc_20_sis <- data.frame(age_Focal = "At 20 years",
                         number = seq(0,6),
                         prob = convoluion_nth(2,list(pdf_f20_ys,pdf_f20_os)),
                         kin = "Sisters")
acc_20_sis$X <- acc_20_sis$number*acc_20_sis$prob
acc_20_sis$X2 <- acc_20_sis$number^2*acc_20_sis$prob
mean_SIS_f20 <- sum(acc_20_sis$X)
var_SIS_f20 <- sum(acc_20_sis$X2) - (mean_SIS_f20)^2
acc_20_sis$expectation <- mean_SIS_f20
acc_20_sis$variance <- var_SIS_f20

############################################## PDF of accumulated kin when Focal is 50
## Younger
df_list_ys_50 <- df_list_ys_50 %>% as.data.frame()
df_list_ys_50 %>% head()
conv_list_ys <- list()
for(age1 in df_list_ys_50$age %>% unique()){
  p <- df_list_ys_50 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_ys[[(1+length(conv_list_ys))]] <- list("data" = p, "age" = age1)
}
conv_list_new_ys <- lapply(1:length(conv_list_ys),
                           function(x){ conv_list_ys[[x]][["data"]] })
pdf_f50_ys <- convoluion_nth(length(conv_list_new_ys), conv_list_new_ys)
## Older
df_list_os_50 <- df_list_os_50 %>% as.data.frame()
df_list_os_50 %>% head()
conv_list_os <- list()
for(age1 in df_list_os_50$age %>% unique()){
  p <- df_list_os_50 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_os[[(1+length(conv_list_os))]] <- list("data" = p, "age" = age1)
}
conv_list_new_os <- lapply(1:length(conv_list_os), function(x){ conv_list_os[[x]][["data"]] })
pdf_f50_os <- convoluion_nth(length(conv_list_new_os), conv_list_new_os)

## df combine
acc_50_sis <- data.frame(age_Focal = "At 50 years",
                        number = seq(0,6),
                        prob = convoluion_nth(2,list(pdf_f50_ys,pdf_f50_os)),
                        kin = "Sisters")
acc_50_sis$X <- acc_50_sis$number*acc_50_sis$prob
acc_50_sis$X2 <- acc_50_sis$number^2*acc_50_sis$prob
mean_SIS_f50 <- sum(acc_50_sis$X)
var_SIS_f50 <- sum(acc_50_sis$X2) - (mean_SIS_f50)^2
acc_50_sis$expectation <- mean_SIS_f50
acc_50_sis$variance <- var_SIS_f50

############ Same for aunts

############################################## PDF of accumulated kin when Focal is 20
## Younger
df_list_cya_20 <- df_list_cya_20 %>% as.data.frame()
df_list_cya_20 %>% head()
conv_list_cya <- list()
for(age1 in df_list_cya_20$age %>% unique()){
  p <- df_list_cya_20 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_cya[[(1+length(conv_list_cya))]] <- list("data" = p, "age" = age1)
}
conv_list_new_cya <- lapply(1:length(conv_list_cya),
                           function(x){ conv_list_cya[[x]][["data"]] })
pdf_f20_cya <- convoluion_nth(length(conv_list_new_cya), conv_list_new_cya)
## Older
df_list_coa_20 <- df_list_coa_20 %>% as.data.frame()
df_list_coa_20 %>% head()
conv_list_coa <- list()
for(age1 in df_list_coa_20$age %>% unique()){
  p <- df_list_coa_20 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_coa[[(1+length(conv_list_coa))]] <- list("data" = p, "age" = age1)
}
conv_list_new_coa <- lapply(1:length(conv_list_coa), function(x){ conv_list_coa[[x]][["data"]] })
pdf_f20_coa <- convoluion_nth(length(conv_list_new_coa), conv_list_new_coa)

## df combine
acc_20_cousins <- data.frame(age_Focal = "At 20 years",
                         number = seq(0,6),
                         prob = convoluion_nth(2,list(pdf_f20_cya,pdf_f20_coa)),
                         kin = "Cousins")
acc_20_cousins$X <- acc_20_cousins$number*acc_20_cousins$prob
acc_20_cousins$X2 <- acc_20_cousins$number^2*acc_20_cousins$prob
mean_Cous_f20 <- sum(acc_20_cousins$X)
var_Cous_f20 <- sum(acc_20_cousins$X2) - (mean_Cous_f20)^2
acc_20_cousins$expectation <- mean_Cous_f20
acc_20_cousins$variance <- var_Cous_f20

############################################## PDF of accumulated kin when Focal is 50
## Younger
df_list_cya_50 <- df_list_ya_50 %>% as.data.frame()
df_list_cya_50 %>% head()
conv_list_cya <- list()
for(age1 in df_list_cya_50$age %>% unique()){
  p <- df_list_cya_50 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_cya[[(1+length(conv_list_cya))]] <- list("data" = p, "age" = age1)
}
conv_list_new_cya <- lapply(1:length(conv_list_cya),
                           function(x){ conv_list_cya[[x]][["data"]] })
pdf_f50_cya <- convoluion_nth(length(conv_list_new_cya), conv_list_new_cya)
sum(pdf_f50_cya)
## Older
df_list_coa_50 <- df_list_coa_50 %>% as.data.frame()
df_list_coa_50 %>% head()
conv_list_coa <- list()
for(age1 in df_list_coa_50$age %>% unique()){
  p <- df_list_coa_50 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_coa[[(1+length(conv_list_coa))]] <- list("data" = p, "age" = age1)
}
conv_list_new_coa <- lapply(1:length(conv_list_coa), function(x){ conv_list_coa[[x]][["data"]] })
pdf_f50_coa <- convoluion_nth(length(conv_list_new_coa), conv_list_new_coa)

## df combine
acc_50_cousins <- data.frame(age_Focal = "At 50 years",
                         number = seq(0,6),
                         prob = convoluion_nth(2,list(pdf_f50_cya,pdf_f50_coa)),
                         kin = "Cousins")
acc_50_cousins$X <- acc_50_cousins$number*acc_50_cousins$prob
acc_50_cousins$X2 <- acc_50_cousins$number^2*acc_50_cousins$prob
mean_Cous_f50 <- sum(acc_50_cousins$X)
var_Cous_f50 <- sum(acc_50_cousins$X2) - (mean_Cous_f50)^2
acc_50_cousins$expectation <- mean_Cous_f50
acc_50_cousins$variance <- var_Cous_f50


### Combine all data frames
pdf_sibs_sim <- rbind(acc_20_sis,acc_50_sis,acc_20_cousins,acc_50_cousins)
## overlay the expectations using vlines
err_dat <- pdf_sibs_sim %>% filter(number==0)
err_dat
err_dat$prob <- 0.02 ## for plotting overlay
err_dat$kin1 <- factor(err_dat$kin, levels = c("Sisters","Cousins"))
pdf_sibs_sim$kin1 <- factor(pdf_sibs_sim$kin, levels = c("Sisters","Cousins"))


## Figure 3 in manuscript
pdf_sis_plot <- pdf_sibs_sim %>%
  filter(age_Focal %in% c("At 50 years", "At 20 years")) %>%
  ggplot(aes(x = number, y = prob , color = age_Focal, fill = age_Focal )) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.25) +
  facet_grid(age_Focal ~ kin, scales = "free_x") + theme_bw() +
  geom_vline(aes(xintercept = expectation, color = age_Focal)) +
  geom_errorbarh(data = err_dat,
                 aes(xmin = ifelse(expectation - variance^(0.5) < 0, 0, expectation - variance^(0.5)),
                     xmax = expectation + variance^(0.5),
                     color = age_Focal),
                 height = 0.015,
                 size = 0.6)  +
  scale_fill_manual(values = c("At 50 years" = "#FC8D62", "At 20 years" = "#66C2A5"))+
  scale_color_manual(values = c("At 50 years" = "#FC8D62", "At 20 years" = "#66C2A5")) +
  ylab("Probability") + xlab("Number of kin")  +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size=14)) +
  theme(strip.text = element_text(size = 14)) +
  theme(legend.position = "none") + scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  ylim(c(0,1))

pdf_sis_plot
#ggsave(paste0(fig_out,"/PDF_sis_cousins.png"), pdf_sis_plot, width = 9, height = 7)

rbind(acc_20_sis,acc_50_sis,acc_20_aunts,acc_50_aunts) %>% filter(age_Focal == "At 20 years")


## check
pdf_sibs_sim %>% dplyr::filter(kin == "Cousins", age_Focal == "At 50 years") %>%
  dplyr::mutate(X3 = number^3*prob,
                s2 = sum(X2),
                s3 = sum(X3),
                m3 = s3 + 2*expectation^3 - 3*expectation*s2,
                skew = m3/((variance^(0.5))^(3)))

mean(a20)
sd(a20)
0.7199918^(0.5)

pdf_sibs_sim %>% filter(number == 0)
pdf_sibs_sim$sd <- pdf_sibs_sim$variance^(0.5)

pdf_sis_plot <- pdf_sibs_sim %>% dplyr::group_by(kin, age_Focal) %>%
  dplyr::summarise(X3 = number^3*prob,
                   s2 = sum(X2),
                   s3 = sum(X3),
                   m3 = s3 + 2*expectation^3 - 3*expectation*s2,
                   skew = m3/((variance^(0.5))^(3)),
                   number = number,
                   prob = prob,
                   expectation = expectation,
                   sd = variance^(0.5)) %>%
  ungroup() %>%
  filter(age_Focal %in% c("At 50 years", "At 20 years")) %>%
  ggplot(aes(x = number, y = prob , color = age_Focal, fill = age_Focal )) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.25) +
  facet_grid(age_Focal ~ kin, scales = "free_x") + theme_bw() +
  scale_fill_manual(values = c("At 50 years" = "#FC8D62", "At 20 years" = "#66C2A5"))+
  scale_color_manual(values = c("At 50 years" = "#FC8D62", "At 20 years" = "#66C2A5")) +
  ylab("Probability") + xlab("Number of kin")  +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size=14)) +
  theme(strip.text = element_text(size = 14)) +
  theme(legend.position = "none") + scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  ylim(c(0,1)) +  geom_text(
    aes(x = 4.5 ,y = 1, label = paste0("mean", "==", round(expectation,2))
    ),parse = TRUE, hjust = 0,vjust = 1, size = 5 ) +
  geom_text(
      aes(x = 4.5,y = 0.9, label = paste0("sd", "==", round(sd,2))
      ),parse = TRUE, hjust = 0,vjust = 1 , size = 5) +
  geom_text(
    aes(x = 4.5, y = 0.8, label = paste0("skew", "==", round(skew,2))
    ),parse = TRUE, hjust = 0,vjust = 1 , size = 5)


pdf_sis_plot
ggsave(paste0(fig_out,"/PDF_sis_cousins.png"), pdf_sis_plot, width = 9, height = 7)

text_df
pdf_sibs_sim %>% dplyr::group_by(kin, age_Focal) %>%
  dplyr::summarise(X3 = number^3*prob,
                   s2 = sum(X2),
                   s3 = sum(X3),
                   m3 = s3 + 2*expectation^3 - 3*expectation*s2,
                   skew = m3/((variance^(0.5))^(3)),
                   number = number,
                   prob = prob) %>%
  ungroup()
