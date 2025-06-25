

library(ggplot2)
library(dplyr)
library(reshape2)

fig_out <- here::here("Manuscript_Figures_code/Figures")

### load in matrix data
f_mat <- readr::read_rds("data/f_mat.Rds")
u_mat <- readr::read_rds("data/u_mat.Rds")

########### We plot the accumulated kin of Focal for each age of Focal

f_int <- which(mothers_age(u_mat,f_mat)>0) %>% length() ## fertility interval
nm <- which(mothers_age(u_mat,f_mat)>0)[1] ## min age of fertile period (n underbar)
nmm <- which(mothers_age(u_mat,f_mat)>0)[length(which(mothers_age(u_mat,f_mat)>0))] ## max age of fertile period (n bar)

########################## Mother #######################################
mother_list <- list()
pp_list <- list()
for(age_focal in 0:99){
  mother_conv <- list()
  for(age_mother in min((age_focal + nm-1 ),100) : min(age_focal + nmm+1, 100) ){

    mother <- unconditional_mother_pdf(age_focal, age_mother, u_mat, f_mat, 7)
    mother <- mother$prob[2]
    mother_conv[[(1+length(mother_conv))]] <- mother}
  pp <- sum(unlist(lapply(1:length(mother_conv), function(x){mother_conv[[x]]})))
  pp_list[[(1+length(pp_list))]] <- pp
  mother_df <- data.frame(number = c(0, 1),
                          prob = c(1-pp, pp),
                          age_focal = age_focal)
  mother_list[[(1+length(mother_list))]] <- mother_df
}
full_mf <- do.call("rbind", mother_list)
full_mf$prob <- ifelse(full_mf$prob < 0, 0, full_mf$prob)
full_mf$X <- full_mf$number*full_mf$prob
full_mf$X2 <- full_mf$number^2*full_mf$prob
full_mf$X3 <- full_mf$number^3*full_mf$prob

full_mf %>% group_by(age_focal) %>%
  summarise(check = sum(prob)) %>%
  ungroup()

## Save the data
#saveRDS(full_mf, "data/MUM_pdf_cum.Rds")
## Read in data
full_mf <- readRDS(paste0("data", "/MUM_pdf_cum.Rds"))
df_sim_mum <- readRDS(paste0("data", "/MUM_sim_py2.Rds"))

## Put ABM model into data frame
sim_mum1 <- df_sim_mum %>%
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
## Put ABM model into data frame
sim_mum2 <- sim_mum1 %>%
  dplyr::transmute(age_Focal = Age_Focal,
                   expectation = expectation,
                   variance = variance,
                   skew = skew,
                   number = number,
                   prob = prob,
                   kin = "Mother",
                   method = "Microsimulation",
                   inv_cum_p_95 = ifelse(cum_p >= 0.95 , number, 20 ),
                   inv_cum_p_05 = ifelse(cum_p >= 0.05 , number, 20 )) %>%
  group_by(age_Focal) %>%
  summarise(number_q95 = min(inv_cum_p_95),
            number_q05 = min(inv_cum_p_05),
            age_Focal = age_Focal,
            expectation = expectation,
            variance = variance,
            skew = skew,
            kin = "Mother",
            method = "Microsimulation",
            number = number,
            prob = prob) %>%
  ungroup()


## Put matrix model into data frame
pdf_mum_df <- full_mf %>%
  dplyr::group_by(age_focal) %>%
  dplyr::summarise(expectation = sum(X),
                   s2 = sum(X2),
                   variance = s2 - expectation^2,
                   s2 = sum(X2),
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
                   kin = "Mother",
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
            kin = "Mother",
            method = "PMF model",
            number = number,
            prob = prob) %>%
  ungroup()


pdf_mum_df$variance
### Daughters ###########
acc_dau_list <- list()
for(age_focal in 1:99){
  dau_conv <- list()
  for(age_dau in 0:(-1+age_focal)){
    DA <- unconditional_Focals_age_Dau(age_focal, age_dau, u_mat, f_mat, 7)
    dau_conv[[(1+length(dau_conv))]] <- DA
  }
  #pp <- fourier_fft_nth(length(dau_conv), dau_conv) ## Note that direct inverse fft slightly faster than convolution of list
  pp <- convoluion_nth(length(dau_conv), dau_conv)
  DA_df <- data.frame(number = seq(0,6),
                      prob = pp,
                      age_focal = age_focal)
  acc_dau_list[[(1+length(acc_dau_list))]] <- DA_df
}
acc_dau_list <- do.call("rbind", acc_dau_list)
## Check pdfs sum to 1 (over all age of Focal)
acc_dau_list %>% dplyr::group_by(age_focal) %>%
  dplyr::summarise(check = sum(prob)) %>%
  dplyr::ungroup()
acc_dau_list$X <- acc_dau_list$number*acc_dau_list$prob
acc_dau_list$X2 <- acc_dau_list$number^2*acc_dau_list$prob
acc_dau_list$X3 <- acc_dau_list$number^3*acc_dau_list$prob
## save data
#saveRDS(acc_dau_list, "data/DAU_pmf_accum.Rds")


df_sim_dau <- readRDS(paste0("data", "/DAU_sim_py2.Rds"))
acc_dau_list <- readRDS(paste0("data", "/DAU_pmf_accum.Rds"))
## ABM data as frame
df_sim_dau1 <- df_sim_dau %>% dplyr::filter(Age_Focal >=0 ) %>%
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

df_sim_dau2 <- df_sim_dau1 %>%
  dplyr::transmute(age_Focal = Age_Focal,
                   expectation = expectation,
                   variance = variance,
                   skew = skew,
                   kin = "Daughters",
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
            kin = "Daughters",
            method = "Microsimulation",
            number = number,
            prob = prob) %>%
  ungroup()
df_sim_dau2
## Matrix data as frame
pdf_dau_df1 <- acc_dau_list %>%
  dplyr::group_by(age_focal) %>%
  dplyr::summarise(expectation = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   variance = s2 - expectation^2,
                   m3 = s3 + 2*expectation^3 - 3*expectation*s2,
                   skew = m3/((variance^(0.5))^(3)),
                   prob = prob,
                   cum_p = cumsum(prob),
                   number = number) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_Focal = age_focal,
                   expectation = expectation,
                   variance = variance,
                   skew = skew,
                   kin = "Daughters",
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
            kin = "Daughters",
            method = "PMF model",
            number = number,
            prob = prob) %>%
  ungroup()

# Combine daughters and mother

head(pdf_mum_df)
head(sim_mum2)
head(pdf_dau_df1)
head(df_sim_dau2)

full_md_df <- rbind(pdf_dau_df1, df_sim_dau2 , pdf_mum_df , sim_mum2)
full_md_df$kin1 <- factor(full_md_df$kin, levels = c("Daughters","Mother"))

newMDDF <- full_md_df %>% dplyr::select(age_Focal, expectation, variance, kin, method, number_q05, number_q95, skew)
newMDDF <- newMDDF %>% group_by(age_Focal, kin, method) %>%
  summarise(expectation = expectation, variance = variance, number_q05 = number_q05, number_q95 = number_q95, skew = skew) %>%
  ungroup()
newMDDF
nes_md_df <- newMDDF %>% melt(id = c("age_Focal", "expectation", "variance", "kin", "method", "skew")) %>%
  transmute(age_Focal = age_Focal,
            expectation = expectation,
            sd = sqrt(variance),
            skew = skew,
            kin = kin,
            method = method,
            quantile = variable,
            quantile_value = value)

nes_md_df1 <- nes_md_df %>% melt(id = c("age_Focal", "quantile", "quantile_value", "kin", "method", "skew")) %>%
  transmute(age_Focal = age_Focal,
            quantile = quantile,
            quantile_value = quantile_value,
            kin = kin,
            method = method,
            moment = variable,
            skew = skew,
            value = value)
### Figure 7 in manuscript
head(nes_md_df1)

a
b
a2
b2
ylim.prim2 <- c(0, 1)   # in this example, precipitation
ylim.sec2 <- c(-50, 70)
b2 <- diff(ylim.prim2)/diff(ylim.sec2)
a2 <- ylim.prim2[1] - b2*ylim.sec2[1]

p2 <- nes_md_df1 %>% filter(age_Focal<100, kin == "Mother") %>%
  ggplot(aes(x = age_Focal, color = method)) +
  facet_wrap(~kin,  scales = "free") +
  geom_point(aes(y = value, color = method, shape = moment), size = 1, alpha = 1)  +
  geom_line(aes(y = quantile_value, linetype = quantile, color = method), alpha = 1) +
  geom_point(data = nes_md_df1 %>% dplyr::filter(kin == "Mother", skew <= 70, skew >= -50),
             aes(y = a2 + b2*skew), shape = "o", size = 1) +
  scale_y_continuous(sec.axis = sec_axis(~ (. - a2)/b2, name = "Skew")) +
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
         linetype = guide_legend(order = 2))+ theme(plot.margin = margin(0, 0, 0, 0, "cm"))

p2
ylim.prim <- c(0, 3)   # in this example, precipitation
ylim.sec <- c(0, 40)
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

p1 <- nes_md_df1 %>% filter(age_Focal<100, kin == "Daughters") %>%
  ggplot(aes(x = age_Focal, color = method)) +
  facet_wrap(~kin,  scales = "free") +
  geom_point(aes(y = value, color = method, shape = moment), size = 1, alpha = 1)  +
  geom_line(aes(y = quantile_value, linetype = quantile, color = method), alpha = 1) +
  geom_point(data = nes_md_df1 %>% dplyr::filter(kin == "Daughters", !is.na(skew), skew <= 40, skew >= 0),
             aes(y = a + b*skew), shape = "o", size = 1) +
  scale_y_continuous(sec.axis = sec_axis(~ (. - a)/b, name = "")) +
  scale_color_manual(labels = c("Microsimulation","PMF model"), values = c("black","orange")) +
  scale_linetype_discrete(labels = c("5%","95%")) +
  scale_shape_manual(labels = c("Expectation","Std. Dev."), values = c(15, 4)) +
  theme_bw() + theme(legend.position = "top") + xlab("Age of Focal") +
  ylab("Expectation, St. Dev. Quantiles.") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14)) +
  theme(strip.text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 1,order = 1),
         shape = guide_legend(order = 3),
         linetype = guide_legend(order = 2))+ theme(plot.margin = margin(0, 0, 0, 0, "cm"))


p1
cowplot::plot_grid(plotlist = list(p1, p2))


library(tidyverse)
library(ggpubr)

# first plot

fig_out <- "Revisions_Figures/Figures"
# Create grid
quan_md_plot <- ggpubr::ggarrange(
  p1, p2, # list of plots
  labels = NULL, # labels
  common.legend = TRUE, # COMMON LEGEND
  legend = "top", # legend position
  align = "h", # Align them both, horizontal and vertical
  hjust= 0, vjust= 0,
  nrow = 1 # number of rows
)
quan_md_plot
ggsave(paste0(fig_out,"/ACC_mumdau_Q.png"), quan_md_plot, width = 10, height = 6)


nes_md_df1 %>% dplyr::filter(kin == "Daughters", !is.na(skew), skew <= 40, skew >= 0) %>%
  dplyr::select(skew) %>% max()


### Prob plots

full_md_df$prob[full_md_df$kin == "Mother" & full_md_df$number == 0] <- NA

quan_md_plot <- full_md_df %>% filter(number < 5, age_Focal>0) %>% mutate(number = factor(number)) %>%
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
quan_md_plot
ggsave(paste0(fig_out,"/Prob_mumdau_Q.png"), quan_md_plot, width = 10, height = 6)
