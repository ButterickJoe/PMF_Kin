library(ggplot2)
library(dplyr)
library(reshape2)

fig_out <- "Revisions_Figures/Figures"

### check matrix
f_mat <- readr::read_rds("data/f_mat.Rds")
u_mat <- readr::read_rds("data/u_mat.Rds")
########### Accumulated kin total
f_int <- which(mothers_age(u_mat,f_mat)>0)%>%length()

########## Younger sisters over Focal's life #################
acc_ys_list <- list()
for(age_focal in 1:99){
  ys_conv <- list()
  for(age_ys in max(0,age_focal-f_int-1):(age_focal-1)){
    ## Notice the time differences in function choice
    #YS <- unconditional_mothers_age_YS(age_focal, age_ys, u_mat, f_mat, 7)
    #YS <- YS$prob
    ## Matrix model faster
    YS <- Matrix_func_age_YS(age_focal, age_ys, u_mat, f_mat, 7) ## note both give same results
    ys_conv[[(1+length(ys_conv))]] <- YS
  }
  #pp <- convoluion_nth(length(ys_conv), ys_conv)
  pp <- fourier_fft_nth(length(ys_conv), ys_conv) ## fft quicker than nth convolution function
  ys_df <- data.frame(number = seq(0,6),
                          prob = pp,
                          age_focal = age_focal)
  acc_ys_list[[(1+length(acc_ys_list))]] <- ys_df
}
acc_ys_list <- do.call("rbind", acc_ys_list)
## check sum to one
acc_ys_list %>% dplyr::group_by(age_focal) %>%
  dplyr::summarise(check = sum(prob)) %>%
  dplyr::ungroup()
# save output
#saveRDS(acc_ys_list, "data/YS_pmf_accum.Rds")
############### Older sisters over Focal's life ##################
acc_os_list <- list()
for(age_focal in 0:99){
  os_conv <- list()
  for(age_os in (age_focal+1):min(100,age_focal+f_int+1) ){
    #OS <- unconditional_mothers_age_OS(age_focal, age_os, u_mat, f_mat, 7)
    #OS <- OS$prob
    OS <- Matrix_func_age_OS(age_focal, age_os, u_mat, f_mat, 7)
    os_conv[[(1+length(os_conv))]] <- OS
  }
  pp <- fourier_fft_nth(length(os_conv), os_conv)
  os_df <- data.frame(number = seq(0,6),
                      prob = pp,
                      age_focal = age_focal)
  acc_os_list[[(1+length(acc_os_list))]] <- os_df
}
acc_os_list <- do.call("rbind", acc_os_list)
## check sum to one
acc_os_list %>% dplyr::group_by(age_focal) %>%
  dplyr::summarise(check = sum(prob)) %>%
  dplyr::ungroup()
## save data out
#saveRDS(acc_os_list, "data/OS_pmf_accum.Rds")

########### Read in simulation/PDF data #################################################
sim_df_YS <- readRDS("data/YS_sim_py2.Rds")
mat_df_YS <- readRDS("data/YS_pmf_accum.Rds")
## put into nice data frame
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
                   X3 = number^3*prob,
                   expectation = sum(X),
                   variance = sum(X2) - expectation^2,
                   s2 = sum(X2),
                   s3 = sum(X3),
                   m3 = s3 + 2*expectation^3 - 3*expectation*s2,
                   skew = m3/((variance^(0.5))^(3)),
                   prob = prob,
                   cum_p = cumsum(prob),
                   number = number) %>%
  ungroup()
## put into nice data frame
sim_df_YS3 <- sim_df_YS2 %>%
  dplyr::transmute(age_Focal = Age_Focal,
                   expectation = expectation,
                   variance = variance,
                   skew = skew,
                   kin = "Younger Sisters",
                   method = "Microsimulation",
                   inv_cum_p_95 = ifelse(cum_p >= 0.95 , number, 20 ),
                   inv_cum_p_05 = ifelse(cum_p >= 0.05 , number, 20 ),
                   number = number,
                   prob = prob) %>%
  group_by(age_Focal) %>%
  summarise(number_q95 = min(inv_cum_p_95),
            number_q05 = min(inv_cum_p_05),
            age_Focal = age_Focal,
            expectation = expectation,
            variance = variance,
            skew = skew,
            number = number,
            prob = prob,
            kin = "Younger Sisters",
            method = "Microsimulation")
### Same with matrix output -- put into nice data frame
mat_df_YS$X <- mat_df_YS$number*mat_df_YS$prob
mat_df_YS$X2 <- mat_df_YS$number^2*mat_df_YS$prob
mat_df_YS$X3 <- mat_df_YS$number^3*mat_df_YS$prob
mat_df_YS1 <- mat_df_YS %>%
  dplyr::group_by(age_focal) %>%
  dplyr::summarise(expectation = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   variance = s2 - expectation^2,
                   m3 = s3 + 2*expectation^3 - 3*expectation*s2,
                   skew = m3/((variance^(0.5))^(3)),
                   cum_p = cumsum(prob),
                   number = number,
                   prob = prob) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_Focal = age_focal,
                   expectation = expectation,
                   variance = variance,
                   skew = skew,
                   kin = "Younger Sisters",
                   method = "PMF model",
                   inv_cum_p_95 = ifelse(cum_p >= 0.95 , number, 20 ),
                   inv_cum_p_05 = ifelse(cum_p >= 0.05 , number, 20 ),
                   number = number,
                   prob = prob) %>%
  group_by(age_Focal) %>%
  summarise(number_q95 = min(inv_cum_p_95),
            number_q05 = min(inv_cum_p_05),
            age_Focal = age_Focal,
            expectation = expectation,
            variance = variance,
            skew = skew,
            number = number,
            prob = prob,
            kin = "Younger Sisters",
            method = "PMF model")

## sim and matrix for younger sisters #####
YS_comp_df <- rbind(sim_df_YS3,mat_df_YS1) %>% dplyr::filter(age_Focal != 0)

############################################### Older Sisters #########################
sim_df_OS <- readRDS("data/OS_sim_py2.Rds")
mat_df_OS <- readRDS("data/OS_pmf_accum.Rds")
## put into nice data frame
sim_df_OS2 <- sim_df_OS %>%
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
                   X3 = number^3*prob,
                   expectation = sum(X),
                   variance = sum(X2) - expectation^2,
                   s2 = sum(X2),
                   s3 = sum(X3),
                   m3 = s3 + 2*expectation^3 - 3*expectation*s2,
                   skew = m3/((variance^(0.5))^(3)),
                   prob = prob,
                   cum_p = cumsum(prob),
                   number = number) %>%
  ungroup()
## put into nice data frame
sim_df_OS3 <- sim_df_OS2 %>%
  dplyr::transmute(age_Focal = Age_Focal,
                   expectation = expectation,
                   variance = variance,
                   skew = skew,
                   kin = "Older Sisters",
                   method = "Microsimulation",
                   inv_cum_p_95 = ifelse(cum_p >= 0.95 , number, 20 ),
                   inv_cum_p_05 = ifelse(cum_p >= 0.05 , number, 20 ),
                   number = number,
                   prob = prob) %>%
  group_by(age_Focal) %>%
  summarise(number_q95 = min(inv_cum_p_95),
            number_q05 = min(inv_cum_p_05),
            age_Focal = age_Focal,
            expectation = expectation,
            variance = variance,
            skew = skew,
            number = number,
            prob = prob,
            kin = "Older Sisters",
            method = "Microsimulation")
mat_df_OS$X <- mat_df_OS$number*mat_df_OS$prob
mat_df_OS$X2 <- mat_df_OS$number^2*mat_df_OS$prob
mat_df_OS$X3 <- mat_df_OS$number^3*mat_df_OS$prob
## put into nice data frame
mat_df_OS1 <- mat_df_OS %>%
  dplyr::group_by(age_focal) %>%
  dplyr::summarise(expectation = sum(X),
                   s2 = sum(X2),
                   variance = s2 - expectation^2,
                   s3 = sum(X3),
                   m3 = s3 + 2*expectation^3 - 3*expectation*s2,
                   skew = m3/((variance^(0.5))^(3)),
                   cum_p = cumsum(prob),
                   number = number,
                   prob = prob) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_Focal = age_focal,
                   expectation = expectation,
                   variance = variance,
                   skew = skew,
                   kin = "Older Sisters",
                   method = "PMF model",
                   inv_cum_p_95 = ifelse(cum_p >= 0.95 , number, 20 ),
                   inv_cum_p_05 = ifelse(cum_p >= 0.05 , number, 20 ),
                   number = number,
                   prob = prob) %>%
  group_by(age_Focal) %>%
  summarise(number_q95 = min(inv_cum_p_95),
            number_q05 = min(inv_cum_p_05),
            age_Focal = age_Focal,
            expectation = expectation,
            variance = variance,
            skew = skew,
            number = number,
            prob = prob,
            kin = "Older Sisters",
            method = "PMF model")

## combine older sisters sim and matirx
OS_comp_df <- rbind(sim_df_OS3,mat_df_OS1)

####### combine sim and martrix for both younger and older ############
full_sib_comp <- rbind(YS_comp_df, OS_comp_df)
full_sib_comp
full_sib_comp$kin1 <- factor(full_sib_comp$kin, levels = c("Younger Sisters","Older Sisters"))
## create df to take quartiles
newSISDF <- full_sib_comp %>%
  dplyr::select(age_Focal, expectation, variance, kin, method, number_q05, number_q95, skew)
nes_sis_df <- newSISDF %>% reshape2::melt(id = c("age_Focal", "expectation", "variance", "kin", "method", "skew")) %>%
  transmute(age_Focal = age_Focal,
            expectation = expectation,
            sd = sqrt(variance),
            skew = skew,
            kin = kin,
            method = method,
            quantile = variable,
            quantile_value = value)

nes_sis_df1 <- nes_sis_df %>% reshape2::melt(id = c("age_Focal", "quantile", "quantile_value", "kin", "method", "skew")) %>%
  transmute(age_Focal = age_Focal,
            quantile = quantile,
            quantile_value = quantile_value,
            kin = kin,
            method = method,
            moment = variable,
            skew = skew,
            value = value)
## Figure 4 in manuscript

ylim.ys <- c(0, 2)   # in this example, precipitation
ylim.ys_sk <- c(0, 10)
b_ys <- diff(ylim.ys)/diff(ylim.ys_sk)
a_ys <- ylim.ys[1] - b_ys*ylim.ys_sk[1]

pys <- nes_sis_df1 %>% filter(age_Focal < 100, kin == "Younger Sisters") %>%
  ggplot(aes(x = age_Focal, color = method)) +
  facet_wrap(~kin) +
  geom_point(aes(y = value, color = method, shape = moment), size = 1, alpha = 1)  +
  geom_line(aes(y = quantile_value, linetype = quantile, color = method), alpha = 1) +
  geom_point(data = nes_sis_df1 %>% dplyr::filter(kin == "Younger Sisters", skew <= 10, skew >= 0),
             aes(y = a_ys + b_ys*skew), shape = "o", size = 1)  +
  scale_y_continuous(sec.axis = sec_axis(~ (. - a_ys)/b_ys, name = "")) +
  scale_color_manual(labels = c("Microsimulation","PMF model"), values = c("black","orange")) +
  scale_linetype_discrete(labels = c("5%","95%")) +
  scale_shape_manual(labels = c("Expectation","Std. Dev."), values = c(15, 4)) +
  theme_bw() + theme(legend.position = "top") + xlab("Age of Focal") +
  ylab("Expectation, St. Dev., Quantiles.") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14)) +
  theme(strip.text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 1,order = 1),
         shape = guide_legend(order = 3),
         linetype = guide_legend(order = 2)) + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
pys
nes_sis_df1 %>% filter(age_Focal < 100, kin == "Younger Sisters")

ylim.os <- c(0, 2)   # in this example, precipitation
ylim.os_sk <- c(0, 50)
b_os <- diff(ylim.os)/diff(ylim.os_sk)
a_os <- ylim.os[1] - b_os*ylim.os_sk[1]

pos <- nes_sis_df1 %>% filter(age_Focal<100, kin == "Older Sisters") %>%
  ggplot(aes(x = age_Focal, color = method)) +
  facet_wrap(~kin) +
  geom_point(aes(y = value, color = method, shape = moment), size = 1, alpha = 1)  +
  geom_line(aes(y = quantile_value, linetype = quantile, color = method), alpha = 1) +
  geom_point(data = nes_sis_df1 %>% dplyr::filter(kin == "Older Sisters", skew <= 50, skew >= 0),
             aes(y = a_os + b_os*skew), shape = "o", size = 1)  +
  scale_y_continuous(sec.axis = sec_axis(~ (. - a_os)/b_os, name = "Skew")) +
  scale_color_manual(labels = c("Microsimulation","PMF model"), values = c("black","orange")) +
  scale_linetype_discrete(labels = c("5%","95%")) +
  scale_shape_manual(labels = c("Expectation","Std. Dev."), values = c(15, 4)) +
  theme_bw() + theme(legend.position = "top") + xlab("Age of Focal") +
  ylab("") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14)) +
  theme(strip.text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 1,order = 1),
         shape = guide_legend(order = 3),
         linetype = guide_legend(order = 2)) + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
pos
ggsave(paste0(fig_out,"/ACC_sis_Q.png"), quan_sis_plot, width = 10, height = 6)
quan_sis_plot <- ggpubr::ggarrange(
  pys, pos, # list of plots
  labels = NULL, # labels
  common.legend = TRUE, # COMMON LEGEND
  legend = "top", # legend position
  align = "h", # Align them both, horizontal and vertical
  hjust= 0, vjust= 0,
  nrow = 1 # number of rows
)

### Prob plots]


full_sib_comp %>%
  dplyr::select(age_Focal, prob, number, kin, method) %>%
  reshape2::melt(id = c("age_Focal","kin","method")) %>%
  transmute()


quan_sis_plot <- full_sib_comp %>% filter(number < 5) %>% mutate(number = factor(number)) %>%
  dplyr::select(age_Focal, prob, number, kin, method) %>%
  filter(age_Focal<100) %>%
  ggplot() +
  facet_wrap(~kin) +
  geom_line(aes(x = age_Focal, y = prob, color = method, linetype = number ), size = 1)  +
  scale_color_manual(labels = c("Microsimulation","PMF model"), values = c("black","orange")) +
  theme_bw() + theme(legend.position = "top") + xlab("Age of Focal") +
  ylab("Probability Focal experiences # kin") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14)) +
  theme(strip.text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 1,order = 1),
         shape = guide_legend(order = 3),
         linetype = guide_legend(order = 2))


ggsave(paste0(fig_out,"/Prob_acc_sis_Q.png"), quan_sis_plot, width = 10, height = 6)
