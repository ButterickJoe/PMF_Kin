library(ggplot2)
library(dplyr)

fig_out <- here::here("Revisions_Figures/Figures")

#### Second: load in rates

### check matrix
f_mat <- readr::read_rds("data/f_mat.Rds")
u_mat <- readr::read_rds("data/u_mat.Rds")
f_int <- which(mothers_age(u_mat,f_mat)>0)%>%length()

############################################# Age-specific dists of kin when Focal is 20 ##############################
### Younger sisters
df_list_ys <- list()
for(i in 0:19){
  df <- ys_pmf(20, i, u_mat, f_mat, 7)
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
dat_20 <- joe_sis_20 %>% dplyr::mutate(xint = 20)
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
dat_50 <- joe_sis_50 %>% dplyr::mutate(xint = 50)

######## Aunts #######################################################################################################
############################################# Age-specific dists of kin when Focal is 20 ##############################
unconditional_grans_age_YA(20, 60, u_mat, f_mat, 7)
Matrix_func_age_YA(20, 60, u_mat, f_mat, 7)
### Younger aunts
df_list_ya_20 <- list()
for(i in 10:80){
  df <- ya_PMF(20, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_ya_20[[(1+length(df_list_ya_20))]] <- df1
}
df_list_ya_20 <- do.call("rbind", df_list_ya_20) %>% as.data.frame()

## plot to Caswell ###############################################################
df_list_ya_20
df_list_Cas_ya_20 <- younger_aunts_dist(u_mat, f_mat)
to_add_Cas_YA <- data.frame(age_kin = seq(0,100),
                         expected_kin = df_list_Cas_ya_20[1:101,22])
to_add_Cas_YA %>% head()

df_list_ya_20 %>%
  dplyr::mutate(X = number*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(expected_kin = sum(X)) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age, expected_kin = expected_kin) %>%
  ggplot(aes(x = age_kin, y = expected_kin)) + geom_point() +
  geom_line(data = to_add_Cas_YA, aes(y = expected_kin))
################################################################################

### Older aunts
df_list_oa_20 <- list()
for(i in 33:100){
  df <- oa_PMF(20, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_oa_20[[(1+length(df_list_oa_20))]] <- df1
}
df_list_oa_20 <- do.call("rbind", df_list_oa_20) %>% as.data.frame()

## plot to Caswell ###############################################################
df_list_oa_20
df_list_Cas_oa_20 <- older_aunts_dist(u_mat, f_mat)
to_add_Cas <- data.frame(age_kin = seq(0,100),
                         expected_kin = df_list_Cas_oa_20[1:101,23])
to_add_Cas %>% head()

df_list_oa_20 %>%
  dplyr::mutate(X = number*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(expected_kin = sum(X)) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age, expected_kin = expected_kin) %>%
  ggplot(aes(x = age_kin, y = expected_kin)) + geom_point() +
  geom_line(data = to_add_Cas, aes(y = expected_kin))
################################################################################

joe_aunt_20 <- data.frame()
for(age1 in 10:100){
  df <- rbind(df_list_ya_20,df_list_oa_20) %>% dplyr::filter(age == age1)

  if(nrow(df)==7){prob = df$prob}
  else{prob = convoluion_nth(2, list(df$prob[1:7],df$prob[8:14]))}
  temp <- data.frame(number = seq(0,6),
                     prob = prob,
                     age = age1)
  joe_aunt_20 <- rbind(joe_aunt_20, temp)
}

joe_aunt_20 <- joe_aunt_20 %>%
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
                kin = "Aunts")
dat_20_aunt <- joe_aunt_20 %>% dplyr::mutate(xint = 20)
### When Focal is 50
### Younger aunts
df_list_ya_50 <- list()
for(i in 45:100){
  df <- ya_PMF(50, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_ya_50[[(1+length(df_list_ya_50))]] <- df1
}
df_list_ya_50 <- do.call("rbind", df_list_ya_50) %>% as.data.frame()
### Older aunts
df_list_oa_50 <- list()
for(i in 65:100){
  df <- oa_PMF(50, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_oa_50[[(1+length(df_list_oa_50))]] <- df1
}
df_list_oa_50 <- do.call("rbind", df_list_oa_50) %>% as.data.frame()

joe_aunt_50 <- data.frame()
for(age1 in 40:100){
  df <- rbind(df_list_ya_50,df_list_oa_50) %>% dplyr::filter(age == age1)

  if(nrow(df)==7){prob = df$prob}
  else{prob = convoluion_nth(2, list(df$prob[1:7],df$prob[8:14]))}
  temp <- data.frame(number = seq(0,6),
                     prob = prob,
                     age = age1)
  joe_aunt_50 <- rbind(joe_aunt_50, temp)
}
df
joe_aunt_50 <- joe_aunt_50 %>%
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
                kin = "Aunts")
dat_50_aunt <- joe_aunt_50 %>% dplyr::mutate(xint = 50)



legend_title_a <- "Age of Focal"

scale = 300
## Figure 5 in manuscript
sa_age_plot <- rbind(dat_20,dat_50,dat_20_aunt,dat_50_aunt) %>% dplyr::filter(mean != 0) %>%
  ggplot(aes(x = age, y = mean, color = age_Focal)) +
  facet_wrap(~kin, scales = "free_x")+ geom_point(size = 0.75) + theme_bw() +
  geom_line(aes(y = var^(0.5)), linetype = "dashed", size = 0.75) +
  geom_line(data = rbind(dat_20,dat_50,dat_20_aunt,dat_50_aunt) %>% dplyr::filter(skew < 50),
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

print(sa_age_plot)
sa_age_plot
ggsave(paste0(fig_out,"/Age_dists_sis_aunt_REV.png"), sa_age_plot, width = 9, height = 5)

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
df_list_ya_20 <- df_list_ya_20 %>% as.data.frame()
df_list_ya_20 %>% head()
conv_list_ya <- list()
for(age1 in df_list_ya_20$age %>% unique()){
  p <- df_list_ya_20 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_ya[[(1+length(conv_list_ya))]] <- list("data" = p, "age" = age1)
}
conv_list_new_ya <- lapply(1:length(conv_list_ya),
                           function(x){ conv_list_ya[[x]][["data"]] })
pdf_f20_ya <- convoluion_nth(length(conv_list_new_ya), conv_list_new_ya)
## Older
df_list_oa_20 <- df_list_oa_20 %>% as.data.frame()
df_list_oa_20 %>% head()
conv_list_oa <- list()
for(age1 in df_list_oa_20$age %>% unique()){
  p <- df_list_oa_20 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_oa[[(1+length(conv_list_oa))]] <- list("data" = p, "age" = age1)
}
conv_list_new_oa <- lapply(1:length(conv_list_oa), function(x){ conv_list_oa[[x]][["data"]] })
pdf_f20_oa <- convoluion_nth(length(conv_list_new_oa), conv_list_new_oa)

## df combine
acc_20_aunts <- data.frame(age_Focal = "At 20 years",
                         number = seq(0,6),
                         prob = convoluion_nth(2,list(pdf_f20_ya,pdf_f20_oa)),
                         kin = "Aunts")
acc_20_aunts$X <- acc_20_aunts$number*acc_20_aunts$prob
acc_20_aunts$X2 <- acc_20_aunts$number^2*acc_20_aunts$prob
mean_AUNT_f20 <- sum(acc_20_aunts$X)
var_AUNT_f20 <- sum(acc_20_aunts$X2) - (mean_AUNT_f20)^2
acc_20_aunts$expectation <- mean_AUNT_f20
acc_20_aunts$variance <- var_AUNT_f20

############################################## PDF of accumulated kin when Focal is 50
## Younger
df_list_ya_50 <- df_list_ya_50 %>% as.data.frame()
df_list_ya_50 %>% head()
conv_list_ya <- list()
for(age1 in df_list_ya_50$age %>% unique()){
  p <- df_list_ya_50 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_ya[[(1+length(conv_list_ya))]] <- list("data" = p, "age" = age1)
}
conv_list_new_ya <- lapply(1:length(conv_list_ya),
                           function(x){ conv_list_ya[[x]][["data"]] })
pdf_f50_ya <- convoluion_nth(length(conv_list_new_ya), conv_list_new_ya)
pdf_f50_ya
## Older
df_list_oa_50 <- df_list_oa_50 %>% as.data.frame()
df_list_oa_50 %>% head()
conv_list_oa <- list()
for(age1 in df_list_oa_50$age %>% unique()){
  p <- df_list_oa_50 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_oa[[(1+length(conv_list_oa))]] <- list("data" = p, "age" = age1)
}
conv_list_new_oa <- lapply(1:length(conv_list_oa), function(x){ conv_list_oa[[x]][["data"]] })
pdf_f50_oa <- convoluion_nth(length(conv_list_new_oa), conv_list_new_oa)

## df combine
acc_50_aunts <- data.frame(age_Focal = "At 50 years",
                         number = seq(0,6),
                         prob = convoluion_nth(2,list(pdf_f50_ya,pdf_f50_oa)),
                         kin = "Aunts")
acc_50_aunts$X <- acc_50_aunts$number*acc_50_aunts$prob
acc_50_aunts$X2 <- acc_50_aunts$number^2*acc_50_aunts$prob
mean_AUNT_f50 <- sum(acc_50_aunts$X)
var_AUNT_f50 <- sum(acc_50_aunts$X2) - (mean_AUNT_f50)^2
acc_50_aunts$expectation <- mean_AUNT_f50
acc_50_aunts$variance <- var_AUNT_f50


### Combine all data frames
pdf_sibs_sim <- rbind(acc_20_sis,acc_50_sis,acc_20_aunts,acc_50_aunts)
## overlay the expectations using vlines
err_dat <- pdf_sibs_sim %>% filter(number==0)
err_dat
err_dat$prob <- 0.02 ## for plotting overlay
err_dat$kin1 <- factor(err_dat$kin, levels = c("Sisters","Aunts"))
pdf_sibs_sim$kin1 <- factor(pdf_sibs_sim$kin, levels = c("Sisters","Aunts"))



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
ggsave(paste0(fig_out,"/PDF_sis_aunts.png"), pdf_sis_plot, width = 9, height = 7)

rbind(acc_20_sis,acc_50_sis,acc_20_aunts,acc_50_aunts) %>% filter(age_Focal == "At 20 years")


## check
pdf_sibs_sim %>% dplyr::filter(kin == "Aunts", age_Focal == "At 50 years") %>%
  dplyr::mutate(X3 = number^3*prob,
                s2 = sum(X2),
                s3 = sum(X3),
                m3 = s3 + 2*expectation^3 - 3*expectation*s2,
                skew = m3/((variance^(0.5))^(3)))

mean(a20)
sd(a20)
