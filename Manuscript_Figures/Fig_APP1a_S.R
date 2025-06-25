library(ggplot2)
library(dplyr)

fig_out <- "Revisions_Figures/Figures"

### First: load in the ABM model results
sim_df_YS <- readRDS("data/YS_sim_py2.Rds")
sim_df_OS <- readRDS("data/OS_sim_py2.Rds")

#### Second: load in rates

### check matrix
f_mat <- readr::read_rds("data/f_mat.Rds")
u_mat <- readr::read_rds("data/u_mat.Rds")
f_int <- which(mothers_age(u_mat,f_mat)>0)%>%length()

############################################# Age-specific dists of kin when Focal is 40 ##############################
### Younger sisters
df_list_ys <- list()
for(i in 0:39){
  df <- Matrix_func_age_YS(40, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_ys[[(1+length(df_list_ys))]] <- df1
}
df_list_ys_40 <- do.call("rbind", df_list_ys) %>% as.data.frame()
df_list_ys_40$X <- df_list_ys_40$number*df_list_ys_40$prob
df_list_ys_40$X2 <- df_list_ys_40$number^2*df_list_ys_40$prob
df_list_ys_40$X3 <- df_list_ys_40$number^3*df_list_ys_40$prob

joe_ys_40 <- df_list_ys_40 %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 40 years",
                kin = "Younger Sisters")
joe_ys_40 %>% head()
joe_ys_40 <- joe_ys_40 %>% as.data.frame()
joe_ys_40 %>% ggplot(aes(x = age, y = mean)) + geom_point() + geom_errorbar(aes(ymin = mean-var, ymax = mean+var))
### Older sisters
df_list_os <- list()
for(i in 41:(41+f_int+1)){
  df <- Matrix_func_age_OS(40, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_os[[(1+length(df_list_os))]] <- df1
}
df_list_os_40 <- do.call("rbind", df_list_os) %>% as.data.frame()
df_list_os_40$X <- df_list_os_40$number*df_list_os_40$prob
df_list_os_40$X2 <- df_list_os_40$number^2*df_list_os_40$prob
df_list_os_40$X3 <- df_list_os_40$number^3*df_list_os_40$prob
joe_os_40 <- df_list_os_40 %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 40 years",
                kin = "Older Sisters")

joe_os_40
combined_sis <- rbind(joe_ys_40,joe_os_40)
combined_sis$kin <- factor(combined_sis$kin, levels = c("Younger Sisters","Older Sisters"))
## Figure 2 in manuscript
age_dists_sis_plot <- combined_sis %>% ggplot(aes(x = age, y = mean)) +
  facet_wrap(~kin, scales = "free_x") + geom_point(size = 1, color = "orange") +
  geom_vline(xintercept = 40, linetype = "dashed") + theme_bw() +
  geom_errorbar(aes(ymin = ifelse(mean - var^(0.5) < 0, 0, mean - var^(0.5)), ymax = mean + var^(0.5)), color = "orange", size = 0.25) +
  ylab("Expectation \u00B1 standard deviation") + xlab("Age of kin") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) + theme(strip.text = element_text(size = 14))

age_dists_sis_plot
#ggsave(paste0(fig_out,"/Age_dists_sis.png"), age_dists_sis_plot, width = 9, height = 5)

############################################## PDF of accumulated kin when Focal is 40
## Younger
df_list_ys_40 <- df_list_ys_40 %>% as.data.frame()
df_list_ys_40 %>% head()
conv_list_ys <- list()
for(age1 in df_list_ys_40$age %>% unique()){
  p <- df_list_ys_40 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_ys[[(1+length(conv_list_ys))]] <- list("data" = p, "age" = age1)
}
conv_list_new_ys <- lapply((1+conv_list_ys[[1]][["age"]]):(1+conv_list_ys[[length(conv_list_ys)]][["age"]]),
                           function(x){ conv_list_ys[[x]][["data"]] })
pdf_f40_ys <- convoluion_nth(length(conv_list_new_ys), conv_list_new_ys)

## Older
df_list_os_40 <- df_list_os_40 %>% as.data.frame()
df_list_os_40 %>% head()
conv_list_os <- list()
for(age1 in df_list_os_40$age %>% unique()){
  p <- df_list_os_40 %>% dplyr::filter(age == age1)
  p <- p$prob
  conv_list_os[[(1+length(conv_list_os))]] <- list("data" = p, "age" = age1)
}
conv_list_new_os <- lapply(conv_list_os[[1]][["age"]]:conv_list_os[[length(conv_list_os)]][["age"]], function(x){ conv_list_os[[(x+1-conv_list_os[[1]][["age"]])]][["data"]] })
pdf_f40_os <- convoluion_nth(length(conv_list_new_os), conv_list_new_os)

## df combine
acc_YS_40 <- data.frame(age_Focal = "At 40 years",
                        number = seq(0,6),
                        prob = pdf_f40_ys,
                        kin = "Younger Sisters")
acc_YS_40$X <- acc_YS_40$number*acc_YS_40$prob
acc_YS_40$X2 <- acc_YS_40$number^2*acc_YS_40$prob
#acc_YS_40$X3 <- acc_YS_40$number^3*acc_YS_40$prob
mean_YS <- sum(acc_YS_40$X)
var_YS <- sum(acc_YS_40$X2) - (mean_YS)^2
acc_YS_40$expectation <- mean_YS
acc_YS_40$variance <- var_YS

acc_OS_40 <- data.frame(age_Focal = "At 40 years",
                        number = seq(0,6),
                        prob = pdf_f40_os,
                        kin = "Older Sisters")
acc_OS_40$X <- acc_OS_40$number*acc_OS_40$prob
acc_OS_40$X2 <- acc_OS_40$number^2*acc_OS_40$prob
mean_OS <- sum(acc_OS_40$X)
var_OS <- sum(acc_OS_40$X2) - (mean_OS)^2
acc_OS_40$expectation <- mean_OS
acc_OS_40$variance <- var_OS

## compare PDFs to the simulation
sim_df_YS %>% head()
## filter ABM to compare pdfs when Focal is 40
age_dist_sim_YS <- sim_df_YS %>% dplyr::filter(Age_Focal == 40)
age_dist_sim_YS1 <- age_dist_sim_YS %>%
  dplyr::group_by(Focal_ID, Age_Focal, sim_no) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup()
age_dist_sim_OS <- sim_df_OS %>% dplyr::filter(Age_Focal == 40)
age_dist_sim_OS1 <- age_dist_sim_OS %>%
  dplyr::group_by(Focal_ID, Age_Focal, sim_no) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup()
sim_YS_PDF_40 <- age_dist_sim_YS1 %>%
  dplyr::group_by(number) %>%
  dplyr::summarise(number = mean(number),
            freq_number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(relative_freq = freq_number/sum(freq_number),
         X = relative_freq*number,
         X2 = relative_freq*number^2) %>%
  dplyr::transmute(age_Focal = "At 40 years",
            number = seq(0,6),
            prob = relative_freq,
            kin = "Younger Sisters",
            X = X,
            X2 = X2,
            expectation = sum(X),
            variance = sum(X2) - expectation^2)
sim_YS_PDF_40$method <- "Microsimulation"
sim_OS_PDF_40 <- age_dist_sim_OS1 %>%
  group_by(number) %>%
  summarise(number = mean(number),
            freq_number = n()) %>%
  ungroup() %>%
  mutate(relative_freq = freq_number/sum(freq_number),
         X = relative_freq*number,
         X2 = relative_freq*number^2) %>%
  transmute(age_Focal = "At 40 years",
            number = seq(0,max(number)),
            prob = relative_freq,
            kin = "Older Sisters",
            X = X,
            X2 = X2,
            expectation = sum(X),
            variance = sum(X2) - expectation^2)
sim_OS_PDF_40$method <- "Microsimulation"
acc_YS_40$method <- "PMF model"
acc_OS_40$method <- "PMF model"
### Combine all data frames
pdf_sibs_sim <- rbind(sim_YS_PDF_40,sim_OS_PDF_40,acc_YS_40,acc_OS_40)
## overlay the expectations using vlines
err_dat <- pdf_sibs_sim %>% filter(number==0)
err_dat$prob <- 0.02 ## for plotting overlay
err_dat$kin1 <- factor(err_dat$kin, levels = c("Younger Sisters","Older Sisters"))
pdf_sibs_sim$kin1 <- factor(pdf_sibs_sim$kin, levels = c("Younger Sisters","Older Sisters"))
## Figure 3 in manuscript
pdf_sis_plot <- pdf_sibs_sim %>%
  ggplot(aes(x = number, y =prob, fill = method)) + geom_bar(position = "dodge", stat = "identity", alpha = 0.2) +
  facet_wrap(~kin1)  + theme_bw() +
  geom_vline(aes(xintercept = expectation, color = method, linetype = method), size = 0.6) +
  geom_errorbarh(data = err_dat,
                 aes(xmin = ifelse(expectation - (variance)^(1/2)<0,0,expectation - (variance)^(1/2)),
                     xmax = expectation + (variance)^(1/2),
                     color = method, linetype = method),
                 height = 0.015, size = 0.6)+
  scale_color_manual(labels = c("Microsimulation","PMF model"), values = c("black","orange"))+
  scale_fill_manual(labels = c("Microsimulation","PMF model"), values = c("black","orange"))+
  xlab("Number of sisters") + ylab("Probability") + scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) + theme(strip.text = element_text(size = 14)) +
  theme(legend.position = "top")

#ggsave(paste0(fig_out,"/PDF_sis.png"), pdf_sis_plot, width = 9, height = 5)
pdf_sis_plot

pdf_sibs_sim %>% group_by(age_Focal, kin, method) %>%
  summarise(X3 = number^3*prob,
            s2 = sum(X2),
            s3 = sum(X3),
            m3 = s3 + 2*expectation^3 - 3*expectation*s2,
            skew = m3/((variance^(0.5))^(3)),
            number = number,
            prob = prob) %>%
  ungroup()



pdf_sis_plot <- pdf_sibs_sim %>% dplyr::select(-age_Focal) %>% filter(number < 7)  %>%
  ggplot(aes(x = number, y = prob , color = method, fill = method )) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.25) +
  facet_grid( ~ kin, scales = "free_x") + theme_bw() +
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
  geom_text(data = adds1,
            aes(x = x_text , y = 1,
                label = paste0("mean", "==", round(expectation,2))
            ), parse = TRUE, hjust = 0, vjust = 1, size = 5, show.legend = FALSE) +
  geom_text(data = adds1,
            aes(x = x_text, y = 0.9,
                label = paste0("sd", "==", round(sd,2))
            ), parse = TRUE, hjust = 0, vjust = 1 , size = 5, show.legend = FALSE) +
  geom_text(data = adds1,
            aes(x = x_text, y = 0.8,
                label = paste0("skew", "==", round(skew,2))
            ), parse = TRUE, hjust = 0, vjust = 1 , size = 5, show.legend = FALSE)

pdf_sis_plot
ggsave(paste0(fig_out,"/PDF_sis_APP.png"), pdf_sis_plot, width = 9, height = 5)

adds <- pdf_sibs_sim %>%
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
  transmute(x_text = ifelse(method == "Microsimulation", 2, 4.5),
            expectation = expectation,
            sd = sd,
            skew = skew,
            number = number,
            kin = kin,
            method = method) %>%
  as.data.frame()


adds1

adds1 <- adds %>% filter(number == 0) %>%
  dplyr::select(kin, method, expectation, sd, skew, x_text)
